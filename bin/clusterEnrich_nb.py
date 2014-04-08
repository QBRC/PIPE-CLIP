#!/usr/bin/python
# programmer : bbc
# usage:

import sys
import re
import random
import string
import math
from pysam import *
from pybedtools import BedTool
import argparse as ap
from collections import Counter


def prepare_argparser():
	description = "Find enriched clusters"
	epilog = "For command line options of each command, type %(prog)s COMMAND -h"
	argparser = ap.ArgumentParser(description=description,epilog=epilog)
	argparser.add_argument("-i","--input",dest="infile",type=str,required=True,help="input bam file")
	argparser.add_argument("-n",dest="readCount",help="Use reads in cluster to filter. Default=0",type = int, default = 0)
	argparser.add_argument("-f",dest="fdr",help="Use FDR to filter. Default=0.01. Exclusive to -n",type = float,default = 0.01)
	argparser.add_argument("-s",dest="species",type=str,default="hg19G")
	argparser.add_argument("-b",dest="binsize",type=int,default=25)
	return(argparser)

def freqRank(readCount,rev=False):#readCount is a dictionary
	key = sorted(readCount.keys(),reverse=rev)
	r_rank = {}
	rank = 0
	for i in key:
		rank +=  readCount[i]
		r_rank[i] = rank
		#print i,readCount[i],rank
	return r_rank

def BH(pvalue,pRank,N):
	a = math.log(N/float(pRank))
	q = a + pvalue
	qva = max(-1 * q,0)
	return qva


def nullNB(binCoverage,genomeSize,bin_size):
	#calculate the number of 0,1,2-read bins
	genome_bin_count = int(genomeSize/bin_size)
	read0_bin = genome_bin_count - len(binCoverage)
	reads = []
	testSet = []
	for item in binCoverage:
		reads.append(int(item[6]))
		if int(item[6])>2:
			testSet.append(item)
	reads_freq = dict(Counter(reads))
	read1_bin = reads_freq[1]
	read2_bin = reads_freq[2]
	N = float(read0_bin+read1_bin+read2_bin)
	p0 = read0_bin / N
	p1 = read1_bin / N
	p2 = read2_bin / N
	print "#",read0_bin,read1_bin,read2_bin
	r = (p1*p1)/(2*p0*p2-p1*p1)
	p = p1 / (r*p0)
	del reads_freq[1]
	del reads_freq[2]
	return (r,p,reads_freq,testSet)
	
def nullNB2(binCoverage):
	#calculate the number of 0,1,2-read bins
	reads = []
	testSet = []
	for item in binCoverage:
		reads.append(int(item[6]))
		if int(item[6])>3:
			testSet.append(item.fields)
	reads_freq = dict(Counter(reads))
	read1_bin = reads_freq[1]
	read2_bin = reads_freq[2]
	read3_bin = reads_freq[3]
	print reads_freq
	N = float(read3_bin+read1_bin+read2_bin)
	p3 = read3_bin / N
	p1 = read1_bin / N
	p2 = read2_bin / N
	
	r = (3*p1*p3-4*p2*p2)/(2*p2*p2-3*p1*p3)
	p = ((r+2)*p3) / (3*p2)
	del reads_freq[1]
	del reads_freq[2]
	del reads_freq[3]

	return (r,p,reads_freq,testSet)


def main():
	argparser = prepare_argparser()
	args = argparser.parse_args()
	GenomeAssembly = {'hg19T':501000000,'hg19G':3100000000}
	try:
		bin_coverage = BedTool(args.infile)
	except IOError,message:
		print >> sys.stderr,"Cannot open BAM file.",message
		sys.exit(1)
	
	(null_r,null_p,readsFreq,testSet) = nullNB(bin_coverage,GenomeAssembly[args.species],args.binsize)
	#(null_r,null_p,readsFreq,testSet) = nullNB2(bin_coverage)
	print "#",null_r,null_p
	#null_r=0.02818
	#null_p=0.85657
	filterByFDR(testSet,args.fdr,null_r,null_p,readsFreq)
	print "#chr\tstart\tend\tcluster_name\treads_in_cluster\tstrand\t-log(q)"	


if __name__=="__main__":
	main()
