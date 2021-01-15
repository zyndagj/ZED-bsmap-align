#! /usr/bin/env python
import numpy as np
import sys
import argparse
from os.path import exists, splitext
from subprocess import Popen, PIPE
from multiprocessing import Process, Pipe

# Jawon tiles
#Start	End	Chrom	F CpG metC	F CpG C	R CpG metC	R CpG C
#200	299	1	0	8	0	0
#1000	1099	1	0	18	0	17
#1100	1199	1	0	0	0	3
#1700	1799	1	0	2	0	0
#1800	1899	1	0	19	0	0

contexts = ('CG','CHG','CHH')

def main():
	parser = argparse.ArgumentParser(description="Combine methratio files from bsmap into windows of fixed sizes.")
	parser.add_argument('files',metavar='FILE', type=str, nargs='+', help='Methratio files to combine and window')
	parser.add_argument('-R', dest='reference', metavar='REF',type=str, required=True, help="Reference genome that is indexed as a fai")
	parser.add_argument('-S', dest='tileSize', metavar='N', type=int, help="Window size", default=200)
	args = parser.parse_args()
	if exists(args.reference):
		ext = splitext(args.reference)[1]
		if ext == '.fa' or ext == '.fasta':
			if exists(args.reference+'.fai'):
				inFai = args.reference+'.fai'
			else:
				sys.exit("Please index reference using samtools faidx.")
		else:
			sys.exit("Reference is not a fasta.")
	else:
		sys.exit("Reference file doesn't exist.")
	chromDict, chroms = parseFAI(inFai)
	print "Parsing:"
	for f in args.files:
		print " - %s" % (f)
	samples = map(lambda y: y.split('_')[0], args.files)
	makeTile(args.tileSize, args.files, chromDict, chroms, samples, 2)

def makeTile(tileSize, files, chromDict, chroms, samples, covThresh=1):
	p = []
	pConns = []
	for f in files:
		pConn, cConn = Pipe()
		pConns.append(pConn)
		p.append(Process(target=ratioWorker, args=(f, chromDict, chroms, tileSize, cConn)))
		p[-1].start()
	outFiles = []
	for c in contexts:
		outFiles.append(open(str(tileSize)+'_'+c+'.tab', 'w'))
	writeHeader(outFiles, samples)
	for chrom in chroms:
		## Added to make naming consistent
		if chrom == 'chloroplast':
			outChr = "chrPt"
		elif chrom == 'mitochondrion':
			outChr = "chrMt"
		else:
			outChr = "chr%s"%(chrom)
		##################################
		for start in xrange(0, chromDict[chrom], tileSize):
			end = start+tileSize
			row = [[], [], []] #[CG, CHG, CHH]
			for i in xrange(len(files)):
				rObj = pConns[i].recv()
				for i in xrange(3):
					row[i].append(rObj[i])
			writeRows(outFiles, outChr, start, end, row)
	for of in outFiles:
		of.close()

def writeHeader(outFiles, samples):
	samList = []
	for s in samples:
		samList += [s+'_ratio', s+'_C', s+'_CT']
	outStr = '\t'.join(['chr','start','end']+samList)+'\n'
	for of in outFiles:
		of.write(outStr)

def writeRows(outFiles, chrom, start, end, row):
	for i in range(3):
		outStr = '%s\t%i\t%i'%((chrom, start, end))
		if not np.all(map(np.isnan, row[i])):
			for d in row[i]:
				outStr += '\t'+'\t'.join(map(str, d))
			outFiles[i].write(outStr+'\n')

def parseFAI(inFai):
	chromDict = {}
	chroms = []
	IF = open(inFai,'r')
	for line in IF:
		tmp = line.split('\t')
		chromDict[tmp[0]] = int(tmp[1])
		chroms.append(tmp[0])
	IF.close()
	chroms.sort()
	return (chromDict,chroms)

def ratioWorker(inFile, chromDict, chroms, tileSize, oPipe):
	#chr	pos	strand	context	ratio	eff_CT_count	C_count	CT_count	rev_G_count	rev_GA_count	CI_lower	CI_upper
	#1	212	+	CHH	0.000	1.00	0	1	0	0	0.000	0.793
	IP = Popen(['zcat',inFile], stdout=PIPE)
	IP.stdout.readline()
	fLine = formatLine(IP.stdout.readline())
	empty = [(np.nan,np.nan,np.nan),(np.nan,np.nan,np.nan),(np.nan,np.nan,np.nan)]
	for chrom in chroms:
		for start in xrange(0, chromDict[chrom], tileSize):
			end = start+tileSize
			c = [0,0,0]
			ct = [0,0,0]
			if not fLine:
				oPipe.send(empty)
			else:
				while fLine[1] < end and chrom == fLine[0]:
					con = fLine[2]
					c[con] += fLine[3]
					ct[con] += fLine[4]
					fLine = formatLine(IP.stdout.readline())
					if not fLine: break
				outA = []
				for i in xrange(3):
					if ct[i] == 0:
						outA.append((np.nan,np.nan,np.nan))
					else:
						ratio = np.round(np.float(c[i])/np.float(ct[i]), 3)
						outA.append((ratio, c[i], ct[i]))
				oPipe.send(outA)
	IP.stdout.close()

def formatLine(line):
	if not line: return None
	tmp = line.rstrip('\n').split('\t')
	chr = tmp[0]
	pos = int(tmp[1])-1
	cIndex = contexts.index(tmp[3])
	c = int(tmp[6])
	ct = int(tmp[7])
	return (chr, pos, cIndex, c, ct)

if __name__ == "__main__":
	main()
