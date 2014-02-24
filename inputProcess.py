#!/usr/bin/python
# Programmer : beibei.chen@utsouthwestern.edu
# Usage: Convert SAM to BAM, and sort it
# Input: SAM/BAM
# Output: sorted BAM
# Last mofidication: 19 Dec. 2013


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
		print "This is not a BAM file"
		infile = open(filename)
	magic = infile.read(3)
	if magic == "BAM":
		#print >> sys.stderr,"This is a BAM"
		return True
	else:
		return False


def main():
	try:#Test file integrity
		pysam.view(outname+".bam","-H","-o"+outname+".header")
	except:
		print >> sys.stder, "Cannot read binary header, please check BAM file.)"
		sys.exit(1) #Stop whole analysis
	#If the file is OK, find out if it is BAM or SAM
	if is_BAM(sys.argv[1]): #is binary, BAM file	
		try:#Test file integrity
			pysam.view(outname+".bam","-H","-o"+outname+".header")
		except:
			print >> sys.stder, "Cannot read binary header, please check BAM file.)"
		outname = sys.argv[2] + ".sorted"

	else:#SAM file,converte to BAM first
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
		print >> sys.stderr,"started to sort"
	
	#Sort BAM and get index
	pysam.sort(bamout,outname)
	pysam.index(outname+".bam")
	

if __name__=="__main__":
	main()
