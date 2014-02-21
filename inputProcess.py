#!/usr/bin/python
# Programmer : beibei.chen@utsouthwestern.edu
# Usage: Convert SAM to BAM, and sort it
# Input: SAM/BAM
# Output: sorted BAM
# Last modified: 20 Feb. 2014


import sys
import re
import string
import gzip
import pysam
from pysam import *
import argparse as ap
import mimetypes


def is_BAM(filename):
	try:
		infile  = gzip.open(filename)
	except IOError: # not a gzip file
		print "This is a IO error"
		infile = open(filename)
	magic = infile.read(3)
	if magic == "BAM":
		return True
	else:
		return False

def main():
	if is_BAM(sys.argv[1]): #is binary, BAM file	
		outname = sys.argv[2] + ".sorted"
		pysam.sort(sys.argv[1],outname)
		pysam.index(outname+".bam")
	else:#SAM file
		try:
			inputfile = pysam.Samfile(sys.argv[1])
			
		except IOError,message:
			print >> sys.stderr, "cannot open SAM file",message
			sys.exit(1)
		bamout = sys.argv[2] + ".bam"
		outname = sys.argv[2] + ".sorted"
		outputfile = pysam.Samfile(bamout,'wb',template=inputfile)
		for item in inputfile.fetch():
			outputfile.write(item)
		pysam.sort(bamout,outname)
		pysam.index(outname+".bam")
	

if __name__=="__main__":
	main()
