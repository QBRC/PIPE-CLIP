#!/usr/bin/python
# programmer : bbc
# usage:


import sys
import re
import random
import string
import pysam
from pysam import *
import ghmm
import argparse as ap
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


#############Begin functions###########

def prepare_argparser():
	description = "Filter SAM file by tags"
	epilog = "For command line options of each command, type %(prog)s COMMAND -h"
	argparser = ap.ArgumentParser(description=description, epilog = epilog)
	#IO arguments
	argparser.add_argument("-i","--input",dest = "infile", type = str, required = True, help = "input bam file")
	argparser.add_argument("-o","--output",dest = "outfile", type = str,required = True, help = "output file, default is stdout")
	
	#CIGAR tag arguments
	argparser.add_argument("-m",dest = "matchLength", type = int, help = "Nucleotide number that are mapped to the genome")
	argparser.add_argument("-s",dest = "soft-clipped", type = int, help = "Minimum number of total base marked as 'S'(soft-clipped) in CIGAR tag")
	argparser.add_argument("-H",dest = "hard-clipped", type = int, help = "Maximum number of total base marked as 'H'(hard-clipped) in CIGAR tag")
	argparser.add_argument("-I",dest = "insertion", type = int, help = "Minimum number of total base marked as 'I'(insertion) in CIGAR tag")
	argparser.add_argument("-d",dest = "deletion", type = int, help = "Minimum number of total base marked as 'D'(deletion) in CIGAR tag")
	
	#Other tag arguments
	argparser.add_argument("-n",dest = "mismatch", type = int, default = 2, help = "NM tag number, default=2")
	
	#remove duplication arguments
	argparser.add_argument("-r",dest = "rm_loc", type = int, default = 1, help = "Remove redundant reads with the same mapping starting location ")
	
	return(argparser)

def countMatchNumber(b):
	myList = b
	m = 0
	for i in myList:
		if i[0]==0:
			m += i[1]
	return (m)
def countMatchLength(b):
	myList = b
	m = 0
	for i in myList:
		if i[0]<=1:
			m += i[1]
	return (m)

def countMismatch(b):
	myTag = b.tags
	for i in myTag:
		if i[0]=="NM":
			return (i[1])

def rmdup(duplist):
	mn = 0 #match number
	matchlist = []
	mq = 0 #match quality
	for i in duplist:
		#print i.qname,i.mapq
		m = countMatchNumber(i.cigar)
		q = i.mapq
		if m>mn  or (m==mn and q>mq):
			matchlist = [i]
			mn = m
			mq = q
		elif m == mn and q==mq:
			matchlist.append(i)
	if len(matchlist)>1: #need to select one randomly
		index =  random.randint(0,len(matchlist)-1)
		return matchlist[index]
	else:
		return matchlist[0]


def main():
	
	argparser = prepare_argparser()
	args = argparser.parse_args()
	#print args
	try:
		inputfile = pysam.Samfile(args.infile,"rb")
	except IOError,message:
		print >> sys.stderr, "cannot open SAM file",message
		sys.exit(1)
	#print "finish to read file"
	outname = args.outfile + ".filter.bam"
	outputfile = pysam.Samfile(outname,'wb',template=inputfile)
	coverageName = args.outfile + ".filter.bam.coverage"
	coveragefile = open(coverageName,"wa")
	lenList = []
	filterFile = []
	barcount = []
	barcount.append(inputfile.mapped)
	if args.mismatch or args.matchLength: #At least some CIGAR fileter if applied
		#print "length filter"
		for item in inputfile:
			b= item.cigar
			#print b
			flag = 0
			if b:
				if args.matchLength:
					if countMatchLength(b)>=args.matchLength:
						flag = 1
				if args.mismatch:
					if countMismatch(item)<=args.mismatch:
						flag = flag and 1
					else:
						flag = flag and 0
				if flag: #and b[0][0]==0:
					lenList.append(countMatchLength(b))
					filterFile.append(item)
	else:
		filterFile = inputfile
	barcount.append(len(lenList))
	if args.rm_loc > 0: #user require to remove duplicates
		former = []
		counter = 0
		reductantList = []
		checkNewStart = 1
		lenList = []  # re-make lengList for histogram
		for item in filterFile:
			if checkNewStart:
				former = item
				reductantList = [item]
				checkNewStart = 0
				counter += 1
				continue
			if item.tid == former.tid and item.pos == former.pos: #and (item.is_reverse and former.is_reverse):
				reductantList.append(item)
				counter += 1
			else: #not same with existing starts
				if len(reductantList)<=1: #no reductance
					#print former
					outputfile.write(former)
					lenList.append(countMatchLength(former.cigar))
					#continue
				else:
					remain = rmdup(reductantList)
					#print remain
					lenList.append(countMatchLength(remain.cigar))
					outputfile.write(remain)
				former = item
				reductantList=[item]
				counter += 1
			if counter >= len(filterFile): #reach the last line of file
				if len(reductantList)<=1: #no reductance
					#print item
					lenList.append(countMatchLength(item.cigar))
					outputfile.write(item)
				else:
					remain = rmdup(reductantList)
					#print remain
					lenList.append(countMatchLength(remain.cigar))
					outputfile.write(remain)

	else:
		if len(lenList)==0:
			for item in filterFile:
				lenList.append(countMatchLength(item.cigar))
				outputfile.write(item)
		else:
			for item in filterFile:
				outputfile.write(item)
	outputfile.close()
	barcount.append(len(lenList))
	print >> coveragefile, sum(lenList)
#generate remained reads mapped length distribution
	dis_name = "length_distribution.pdf"	
	pp = PdfPages(dis_name)
	fig = plt.figure(frameon = True)
	ax = fig.add_subplot(111)
	ax.hist(lenList,50,facecolor='green',alpha=0.5)
	ax.set_xlabel('Matched length')
	ax.set_ylabel('Frequency')
	ax.set_xlim(args.matchLength-1,max(lenList)+1)
	plt.savefig(pp,format='pdf')
	pp.close()

#generate barchar of reads number remained after each filter
	sta_name = "filter_statistics.pdf"
	bp = PdfPages(sta_name)
	bfig = plt.figure(frameon = True)
	bax = bfig.add_subplot(111)
	bax.bar([0.9,1.9,2.9],barcount,0.2,color='green')
	bax.set_ylabel('Frequency')
	bax.set_xlim(0.5,3.5)
	bax.set_xticks([1,2,3])
	bax.set_xticklabels(['Mapped','Length&mismatch','rmdup'])
	for i in range(3):
		bax.annotate(barcount[i],(i+1,barcount[i]+0.02*max(barcount)),va='bottom',ha='center')
	plt.savefig(bp,format='pdf')
	bp.close()

if __name__=="__main__":
	main()
