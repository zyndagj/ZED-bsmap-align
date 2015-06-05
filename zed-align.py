#!/usr/bin/env python

from math import ceil
import os
import sys
import argparse
import multiprocessing
import subprocess as sp
import re
#from pprint import pprint
from array import array
from yaml import load, dump

contexts = ('CG','CHG','CHH')

def main():
	fCheck = fileCheck()
	parser = argparse.ArgumentParser(description="Wrapper for Bisulfite Methylation Alignment.")
	parser.add_argument('-R', metavar='FASTA', help='Reference for alignment', required=True, type=fCheck.fasta)
	parser.add_argument('-r1', metavar='FASTQ', help='Single or first fastq from pair', required=True, type=fCheck.fastq)
	parser.add_argument('-r2', metavar='FASTQ', help='Second read', type=fCheck.fastq)
	parser.add_argument('-N', '--name', metavar='STR', help='Name for run')
	parser.add_argument('-U', '--uniq', action='store_true', help="Only use unique alignments")
	parser.add_argument('-q', help="Fastq Quality Encoding (Default: %(default)s)", default=33, type=int)
	parser.add_argument('-C', metavar='Chrom', help="Chromosome to use for checking bisulfite conversion rate")
	parser.add_argument('-S', dest='tileSize', metavar='N', type=int, help="Window size", default=100)
	parser.add_argument('-d', metavar='N', type=int, help="Minimum coverage in tile for methylation to be printed", default=1)
	args = parser.parse_args()
	if not args.name:
		args.name = os.path.splitext(args.r1)[0]
	config = {'bsmap':{}, 'methratio':{}, 'tiles':{}}
	#bsmap section
	config['bsmap']['-a'] = {'value':args.r1, 'description':'R1 input'}
	config['bsmap']['-z'] = {'value':str(args.q), 'description':'Fastq quality encoding'}
	config['bsmap']['-p'] = {'value':str(int(multiprocessing.cpu_count()*1.5)), 'description':'Number of threads'}
	config['bsmap']['-q'] = {'value':'20', 'description':"Quality threshold for trimming 3' ends of reads"}
	#config['bsmap']['-V'] = {'value':'1', 'description':'Print major messages'}
	config['bsmap']['-d'] = {'value':args.R, 'description':'Reference'}
	#config['bsmap']['-o'] = {'value':args.name+".sam", 'description':'Output BAM'}
	#methratio section
	#config['methratio']['-q'] = {'value':'', 'description':'Quiet'}
	config['methratio']['-z'] = {'value':'', 'description':'Report locations with zero methylation'}
	config['methratio']['-r'] = {'value':'', 'description':'Remove duplicate reads'}
	config['methratio']['-d'] = {'value':args.R, 'description':'Reference'}
	config['methratio']['-o'] = {'value':args.name+"_methratio.txt", 'description':'Output methylation ratio file'}
	# Paired specific arguments
	if args.r2:
		config['bsmap']['-b'] = {'value':args.r2, 'description':'R2 input'}
		config['methratio']['-p'] = {'value':'', 'description':'Require propper pairings'}
	if args.uniq:
		config['bsmap']['-r'] = {'value':'0', 'description':'No non-unique hits reported'}
		config['methratio']['-u'] = {'value':'', 'description':'Only use unique alignments'}
	else:
		config['bsmap']['-r'] = {'value':'2', 'description':'non-unique hits reported'}
		config['bsmap']['-w'] = {'value':'20', 'description':'Only 20 equal best hits reported'}
	# Tile Section
	config['tiles']['size'] = {'value':args.tileSize, 'description':'Size of tiles for summarizing methylation'}
	config['tiles']['minCoverage'] = {'value':args.d, 'description':'Minimum Coverage'}
	# Check for Dependencies
	for d in ('bsmap','samtools','methratio.py','bedGraphToBigWig'):
		if not which(d):
			sys.exit("Please add %s to your path\n"%(d))
	# Parse FAI
	fai = args.R+'.fai'
	if not os.path.exists(fai):
		sys.exit("Please make a fai for your reference\n")
	faiDict = ParseFai(fai)
	# run BSMAP
	runBSMAP(config, args.name, args.r2)
	# run methratio.py
	runRatio(config)
	if args.C:
		calcConversion(config, args.C, faiDict)
	makeTile(config, args.name, faiDict)
	makeBigWig(config,fai)
	dump(config, open(args.name+'.yaml','w'), default_flow_style=False)

def calcConversion(config, chrom, faiDict):
	if not chrom in faiDict:
		chromStr = '\n - '.join(faiDict.keys())
		sys.exit("Chromosome: %s not in reference. Please choose a chromosome from:\n - %s"%(chrom, chromStr))
	ratioFile = config['methratio']['-o']['value']
	p = sp.Popen(["grep", chrom, ratioFile], stdout=sp.PIPE).stdout
	cSum = 0
	ctSum = 0
	for line in p:
		tmp = line.split('\t')
		cSum += int(tmp[6])
		ctSum += int(tmp[7])
	percent = round((1.0-float(cSum)/float(ctSum))*100.0, 2)
	config['conversion'] = {}
	config['conversion']['Chromosome'] = {'value':chrom, 'description':'Chromosome to calculate conversion efficiency from. No methylation should be expected on this chromosome.'}
	config['conversion']['C'] = {'value':cSum, 'description':'Number of methylated cytosines'}
	config['conversion']['CT'] = {'value':ctSum, 'description':'Number of un/methylated cytosines'}
	config['conversion']['percent'] = {'value':percent, 'description':'Conversion rate: (1-C/CT)*100'}
	p.close()

def runRatio(config):
	ratioCMD = makeCMD('methratio.py', config, 'methratio')+[config['bsmap_stats']['output']['value']]
	ratioOUT = sp.check_output(ratioCMD, stderr=sp.STDOUT)
	statLine = ratioOUT.split('\n')[-2]
	m = re.match(r".+total\s([0-9]+)\s.+,\s([0-9]+)\s.+age:\s(\w+\.\w+) fold", statLine)
	mappings, covered, coverage = m.groups()
	config['methratio_stats'] = {}
	config['methratio_stats']['mappings'] = {'value':mappings, 'description':'Number of valid mappings'}
	config['methratio_stats']['covered'] = {'value':covered, 'description':'Number of cytosines covered'}
	config['methratio_stats']['coverage'] = {'value':coverage, 'description':'Average coverage fold'}

def runBSMAP(config, name, r2):
	bsmapCMD = makeCMD('bsmap', config, 'bsmap')
	bsP = sp.Popen(bsmapCMD, stderr=sp.PIPE, stdout=sp.PIPE)
	samP = sp.Popen(['samtools','view','-bS@','5','-'], stdin=bsP.stdout, stdout=open(name+'.bam','wb'), stderr=sp.PIPE)
	bsP.stdout.close()
	bsOUT = bsP.stderr.read()
	samP.wait()
	if r2:
		total, aligned, unique, mult = map(int, re.findall(r'pairs:\s+([0-9]+)', bsOUT))
		unit='pairs'
	else:
		total, aligned, unique, mult = map(int, re.findall(r'reads:\s+([0-9]+)', bsOUT))
		unit='reads'
	config['bsmap_stats'] = {}
	config['bsmap_stats']['output'] = {'value':name+".bam", 'description':'Output BAM'}
	config['bsmap_stats']['input'] = {'value':total, 'description':'Total number of %s in input'%(unit)}
	config['bsmap_stats']['aligned'] = {'value':aligned, 'description':'Total number of %s aligned'%(unit)}
	config['bsmap_stats']['unique'] = {'value':unique, 'description':'Total number of %s uniquely aligned'%(unit)}
	config['bsmap_stats']['mult'] = {'value':mult, 'description':'Total number of %s with multiple alignments'%(unit)}

def makeCMD(baseBin, config, section):
	outCMD = [baseBin]
	cSec = config[section]
	for key in cSec.keys():
		outCMD.append(key)
		v = cSec[key]['value']
		if v: outCMD.append(v)
	return outCMD

def ParseFai(inFile):
	'''
	Parses a fa.fai into a python dictionary
	Paramteters
	================================
	inFile	FILE	fai file
	'''
	return dict(map(lambda y: (y[0], int(y[1])), map(lambda y: y.split('\t'), open(inFile,'r').readlines())))

class fileCheck:
	def check(self, file, exts):
		ext = os.path.splitext(file)[1][1:]
		fName = os.path.split(file)[1]
		if not ext in exts:
			raise argparse.ArgumentTypeError("%s not a %s"%(fName, exts[0]))
		if not os.path.exists(file):
			raise argparse.ArgumentTypeError("%s does not exist"%(file))
	def fastq(self, file):
		self.check(file, ['fastq','fq'])
		return file
	def fasta(self, file):
		self.check(file, ['fasta','fa'])
		return file

def makeBigWig(config,fai):
	bedgraphs = config['tiles']['output']['bedgraphs']['value']
	pool = []
	bws = []
	for bg in bedgraphs:
		bw = os.path.splitext(bg)[0]+'.bw'
		bws.append(bw)
		pool.append(sp.Popen(['bedGraphToBigWig',bg,fai,bw]))
	for p in pool:
		p.wait()
	config['bigwigs'] = {'value':bws,'description':'Bigwig versions of bedgraph files for jbrowse to load'}

def makeTile(config, name, faiDict):
# Make sure to do something with the coverage variable
	bgNames = map(lambda x: name+'_'+x+'.bedgraph', contexts)
	config['tiles']['output'] = {\
'bedgraphs':{'value':bgNames, 'description':'Mehtylation ratios for each methylation motif {CG, CHG, CHH} in bedgraph format.'},\
'tab':{'value':name+'.tab', 'description':'Tab delimited file of methylation ratios and coverage for each tile.'}}
	buffer = 100000
	bGs = map(lambda x: open(x, 'w', buffer), bgNames)
	tab = open(name+'.tab', 'w', buffer)
	# Write header
	headStr = '\t'.join(['Chr','Start','End']+[ c+'_'+t for c in contexts for t in ('ratio','C','CT')])
	tab.write(headStr+'\n')
	# Get parameters
	tileSize = config['tiles']['size']['value']
	ratioFile = config['methratio']['-o']['value']
	sortedChroms = sorted(faiDict.keys())
	# start writing by chromosome
	for chrom in sortedChroms:
		offset = int(ceil(faiDict[chrom]/float(tileSize)))
		C = array('H', [0]*(offset*3))
		CT = array('H', [0]*(offset*3))
		p = sp.Popen(["grep", '^'+chrom, ratioFile], stdout=sp.PIPE).stdout
		for line in p:
			chr, pos, cIndex, c, ct = formatLine(line)
			C[offset*cIndex+pos/tileSize] += c
			CT[offset*cIndex+pos/tileSize] += ct
		p.close()
		zCheck = [False, False, False]
		for posIndex in xrange(offset):
			start = posIndex*tileSize
			end = min(start+tileSize,faiDict[chrom])
			tabStr = '%s\t%i\t%i'%(chrom,start,end)
			for cIndex in range(3):
				loc = offset*cIndex+posIndex
				if C[loc]:
					if zCheck[cIndex]:
						outStr = '%i\t0\n'%(start,)
						zCheck[cIndex] = False
						bGs[cIndex].write(outStr)
					ratio = float(C[loc])/float(CT[loc])
					outStr = '%s\t%i\t%i\t%.2f\n'%(chrom,start,end,ratio)
					tabStr += '\t%.2f\t%i\t%i'%(ratio, C[loc], CT[loc])
					bGs[cIndex].write(outStr)
				else:
					if not zCheck[cIndex]:
						outStr = '%s\t%i\t'%(chrom,start)
						zCheck[cIndex] = True
						bGs[cIndex].write(outStr)
					tabStr += '\t0\t%i\t%i'%(C[loc], CT[loc])
			tab.write(tabStr+'\n')
		for cIndex in range(3):
			if zCheck[cIndex]:
				outStr = '%i\t0\n'%(end,)
				bGs[cIndex].write(outStr)
	# Close files
	for bg in bGs:
		bg.close()
	tab.close()

def formatLine(line):
	tmp = line.split('\t')
	chr = tmp[0]
	pos = int(tmp[1])-1
	cIndex = contexts.index(tmp[3])
	c = int(tmp[6])
	ct = int(tmp[7])
	return (chr, pos, cIndex, c, ct)

def which(program):
	def is_exe(fpath):
		return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
	fpath, fname = os.path.split(program)
	if fpath:
		if is_exe(program):
			return program
	else:
		for path in os.environ["PATH"].split(os.pathsep):
			path = path.strip('"')
			exe_file = os.path.join(path, program)
			if is_exe(exe_file):
				return exe_file
	return None

if __name__ == "__main__":
	main()
