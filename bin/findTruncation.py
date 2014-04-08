#!/usr/bin/python
# programmer : bbc
# usage:

import sys
import re
import random
import string
import math
import pysam
from pybedtools import BedTool
import argparse as ap

def main():
	try:
		mapped_bam = pysam.Samfile(sys.argv[1],"rb")
	except IOError,message:
		print >> sys.stderr,"Cannot open BAM file.",message
		sys.exit(1)
	
	print "#chr\tstart\tend\tcluster_name\toffset\tstrand\ttype"
	for item in mapped_bam:
		end = item.pos+1
		if item.is_reverse:
			strand = '-'
		else:
			strand = '+'
		if item.tid>=0:
			chr = mapped_bam.getrname(item.tid)
		else:
			continue
		st = chr +"\t"+str(item.pos) + "\t"+ str(end) +"\t" + item.qname +"\t0\t" + strand+"\ttruncation"
		print st

if __name__=="__main__":
	main()
