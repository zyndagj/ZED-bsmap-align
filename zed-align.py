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
	fCheck = fileCheck() #class for checking parameters
	parser = argparse.ArgumentParser(description="Wrapper for Bisulfite Methylation Alignment.")
	parser.add_argument('-R', metavar='FASTA', help='Reference for alignment', required=True, type=fCheck.fasta)
	parser.add_argument('-r1', metavar='FASTQ', help='Single or first fastq from pair', required=True, type=fCheck.fastq)
	parser.add_argument('-r2', metavar='FASTQ', help='Second read', type=fCheck.fastq)
	parser.add_argument('-O', metavar='STR', help='Output directory (Default: %(default)s)', default='.', type=str)
	parser.add_argument('-N', '--name', metavar='STR', help='Name for run')
	parser.add_argument('-U', '--uniq', action='store_true', help="Only use unique alignments")
	parser.add_argument('-q', help="Fastq Quality Encoding (Default: %(default)s)", default=33, type=int)
	parser.add_argument('-C', metavar='Chrom', help="Chromosome to use for checking bisulfite conversion rate")
	parser.add_argument('-S', dest='tileSize', metavar='N', type=int, help="Window size (Default: %(default)s)", default=100)
	parser.add_argument('-d', metavar='N', type=int, help="Minimum coverage in tile for methylation to be printed (Default: %(default)s - all)", default=1)
	parser.add_argument('--CG', metavar='N', type=int, help="Minimum sites per tile (Default: %(default)s)", default=3)
	parser.add_argument('--CHG', metavar='N', type=int, help="Minimum sites per tile (Default: %(default)s)", default=3)
	parser.add_argument('--CHH', metavar='N', type=int, help="Minimum sites per tile (Default: %(default)s)", default=6)
	args = parser.parse_args()
	######################################################
	# Path Section
	######################################################
	if not args.name:
		args.name = os.path.splitext(args.r1)[0]
	if not os.path.exists(args.O): os.makedirs(args.O)
	outPrefix = os.path.join(args.O, args.name)
	######################################################
	# Arguments Section
	######################################################
	config = {'bsmap':{}, 'methratio':{}, 'tiles':{}}
	#-----------------------------------------------------
	# Arguments for running BSMAP
	#-----------------------------------------------------
	config['bsmap']['-a'] = {'value':args.r1, 'description':'R1 input'}
	config['bsmap']['-z'] = {'value':str(args.q), 'description':'Fastq quality encoding'}
	config['bsmap']['-p'] = {'value':str(int(multiprocessing.cpu_count()*1.5)), 'description':'Number of threads'}
	config['bsmap']['-q'] = {'value':'20', 'description':"Quality threshold for trimming 3' ends of reads"}
	config['bsmap']['-d'] = {'value':args.R, 'description':'Reference'}
	config['bsmap']['-S'] = {'value':'77345', 'description':'Hardcoded random seed for mapping reproducibility'}
	config['bsmap']['-w'] = {'value':'10000', 'description':'Number of candidate seeds to align against'}
	#config['bsmap']['-V'] = {'value':'1', 'description':'Print major messages'}
	#config['bsmap']['-o'] = {'value':args.name+".sam", 'description':'Output BAM'} # default SAM stdout is piped to samtools
	#-----------------------------------------------------
	# Arguments for methratio.py
	#-----------------------------------------------------
	#config['methratio']['-q'] = {'value':'', 'description':'Quiet'}
	config['methratio']['-z'] = {'value':'', 'description':'Report locations with zero methylation'}
	config['methratio']['-r'] = {'value':'', 'description':'Remove duplicate reads'}
	config['methratio']['-d'] = {'value':args.R, 'description':'Reference'}
	config['methratio']['-o'] = {'value':outPrefix+"_methratio.txt", 'description':'Output methylation ratio file'}
	#-----------------------------------------------------
	# Paired specific arguments
	#-----------------------------------------------------
	if args.r2:
		config['bsmap']['-b'] = {'value':args.r2, 'description':'R2 input'}
		config['methratio']['-p'] = {'value':'', 'description':'Require propper pairings'}
	if args.uniq:
		config['bsmap']['-r'] = {'value':'0', 'description':'No non-unique hits reported'}
		config['methratio']['-u'] = {'value':'', 'description':'Only use unique alignments'}
	else:
		config['bsmap']['-r'] = {'value':'2', 'description':'non-unique hits reported'}
		config['bsmap']['-w'] = {'value':'20', 'description':'Only 20 equal best hits reported'}
	#-----------------------------------------------------
	# Tile Section
	#-----------------------------------------------------
	config['tiles']['size'] = {'value':args.tileSize, 'description':'Size of tiles for summarizing methylation'}
	config['tiles']['minCoverage'] = {'value':args.d, 'description':'Minimum Coverage'}
	config['tiles']['CG'] = {'value':args.CG, 'description':'Minimum number of sites per tile'}
	config['tiles']['CHG'] = {'value':args.CHG, 'description':'Minimum number of sites per tile'}
	config['tiles']['CHH'] = {'value':args.CHH, 'description':'Minimum number of sites per tile'}
	######################################################
	# Check for Dependencies
	######################################################
	for d in ('bsmap','samtools','methratio.py','bedGraphToBigWig'):
		if not which(d):
			sys.exit("Please add %s to your path\n"%(d))
	# Parse FAI
	fai = args.R+'.fai'
	if not os.path.exists(fai):
		os.system("samtools faidx %s"%(args.R))
	######################################################
	# Run workflow
	######################################################
	faiDict = ParseFai(fai)
	#-----------------------------------------------------
	# run BSMAP
	#-----------------------------------------------------
	runBSMAP(config, outPrefix, args.r2)
	#-----------------------------------------------------
	# run methratio.py and calculate conversion rate
	#-----------------------------------------------------
	runRatio(config)
	if args.C:
		calcConversion(config, args.C, faiDict)
	#-----------------------------------------------------
	# Make Tiles and Bedgraphs
	#-----------------------------------------------------
	makeTile(config, outPrefix, faiDict)
	#-----------------------------------------------------
	# Make bigWig
	#-----------------------------------------------------
	makeBigWig(config,fai)
	#-----------------------------------------------------
	# Write YAML
	#-----------------------------------------------------
	dump(config, open(outPrefix+'.yaml','w'), default_flow_style=False, width=1000)

def calcConversion(config, chrom, faiDict):
	if not chrom in faiDict:
		chromStr = '\n - '.join(faiDict.keys())
		sys.exit("Chromosome: %s not in reference. Please choose a chromosome from:\n - %s"%(chrom, chromStr))
	ratioFile = config['methratio']['-o']['value']
	p = sp.Popen(["grep", "^%s\s"%chrom, ratioFile], stdout=sp.PIPE).stdout
	cSum = 0
	ctSum = 0
	for line in p:
		tmp = line.split('\t')
		cSum += int(tmp[6])
		ctSum += int(tmp[7])
	percent = round((1.0-float(cSum)/(float(ctSum)+1.0))*100.0, 2)
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

def runBSMAP(config, outPrefix, r2):
	bsmapCMD = makeCMD('bsmap', config, 'bsmap')
	bsP = sp.Popen(bsmapCMD, stderr=sp.PIPE, stdout=sp.PIPE)
	samP = sp.Popen(['samtools','view','-bS@','5','-'], stdin=bsP.stdout, stdout=open(outPrefix+'.bam','wb'), stderr=sp.PIPE)
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
	config['bsmap_stats']['output'] = {'value':outPrefix+".bam", 'description':'Output BAM'}
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

def makeTile(config, outPrefix, faiDict):
# Make sure to do something with the coverage variable
	bgNames = map(lambda x: outPrefix+'_'+x+'.bedgraph', contexts)
	config['tiles']['output'] = {\
'bedgraphs':{'value':bgNames, 'description':'Mehtylation ratios for each methylation motif {CG, CHG, CHH} in bedgraph format.'},\
'tab':{'value':outPrefix+'.tab', 'description':'Tab delimited file of methylation ratios and coverage for each tile.'}}
	buffer = 100000
	bGs = map(lambda x: open(x, 'w', buffer), bgNames)
	tab = open(outPrefix+'.tab', 'w', buffer)
	# Write header
	headStr = '\t'.join(['Chr','Start','End']+[ c+'_'+t for c in contexts for t in ('ratio','C','CT')]) ## old out format
	#headStr = '\t'.join(['Chr','Start','End']+[ c+'_'+t for c in contexts for t in ('ratio','C','CT','#sites')]) ## new out format
	tab.write(headStr+'\n')
	#######################################
	# Get parameters
	#######################################
	tileSize = config['tiles']['size']['value']
	ratioFile = config['methratio']['-o']['value']
	nSitesT = map(lambda y: config['tiles'][y]['value'], contexts)
	sortedChroms = sorted(faiDict.keys())
	#######################################
	# start writing by chromosome
	#######################################
	for chrom in sortedChroms:
		#----------------------------------
		# Create data arrays
		#----------------------------------
		offset = int(ceil(faiDict[chrom]/float(tileSize))) # number of tiles
		C, CT, nSites = makeDataArrays(offset)
		#----------------------------------
		# Read Chrom and populate arrays
		#----------------------------------
		p = sp.Popen(["grep", "^%s\s"%chrom, ratioFile], stdout=sp.PIPE).stdout
		for line in p:
			chr, pos, cIndex, c, ct = formatLine(line)
			index = offset*cIndex+pos/tileSize
			C[index] += c
			CT[index] += ct
			nSites[index] += 1
		p.close()
		# zCheck is true if loc-1 had zero methylation
		zCheck = [False, False, False]
		for posIndex in xrange(offset): # tile index
			start = posIndex*tileSize
			end = min(start+tileSize, faiDict[chrom])
			tabStr = '%s\t%i\t%i'%(chrom,start,end)
			for cIndex in range(3):
				loc = offset*cIndex+posIndex # data index
				tabStr += makeTabStr(C[loc], CT[loc], nSites[loc])
				#-------------------------
				# Generate BG
				#-------------------------
				if C[loc]: # if methylated
					if nSites[loc] < nSitesT[cIndex]:
						if not zCheck[cIndex]:
							bgStr = '%s\t%i\t'%(chrom,start)
							zCheck[cIndex] = True
							bGs[cIndex].write(bgStr)
					else:
						if zCheck[cIndex]: # if previous was 0
							bgStr = '%i\t0\n'%(start,)
							zCheck[cIndex] = False
							bGs[cIndex].write(bgStr)
						ratio = float(C[loc])/float(CT[loc])
						bgStr = '%s\t%i\t%i\t%.2f\n'%(chrom,start,end,ratio)
						bGs[cIndex].write(bgStr)
				else:
					if not zCheck[cIndex]:
						bgStr = '%s\t%i\t'%(chrom,start)
						zCheck[cIndex] = True
						bGs[cIndex].write(bgStr)
				#-------------------------
			tab.write(tabStr+'\n')
		#---------------------------------
		# Write out orphaned zeros
		#---------------------------------
		for cIndex in range(3):
			if zCheck[cIndex]:
				bgStr = '%i\t0\n'%(end,)
				bGs[cIndex].write(bgStr)
	######################################
	# Close files
	######################################
	for bg in bGs:
		bg.close()
	tab.close()

def makeTabStr(C, CT, nSites):
	'''
	Generates a tab-separated string for the .tab file.
	'''
	if C:
		ratio = float(C)/float(CT)
		return '\t%.2f\t%i\t%i\t%i'%(ratio, C, CT, nSites)
	return '\t0\t%i\t%i\t%i'%(C, CT, nSites)

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

def makeDataArrays(offset):
	'''
	Function for creating arrays that keep track of data from
	methratio.py output.

	>>> makeDataArrays(1)
	(array('H', [0, 0, 0]), array('H', [0, 0, 0]), array('H', [0, 0, 0]))
	'''
	C = array('H', [0]*(offset*3))
	CT = array('H', [0]*(offset*3))
	nSites = array('H', [0]*(offset*3)) # max is tile size
	return (C, CT, nSites)

if __name__ == "__main__":
	main()
