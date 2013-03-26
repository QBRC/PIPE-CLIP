#!/usr/bin/python
# programmer : bbc
# usage:


import sys
import re
import random
import string
import gzip
import pysam
from pysam import *
import ghmm
import argparse as ap
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import mimetypes

def is_binary(filename):
	fin = open(filename,'rb')
	try:
		CHUNKSIZE = 1024
		while 1:
			chunk = fin.read(CHUNKSIZE)
			if '\0' in chunk: #found null byte
				return True
			if len(chunk) < CHUNKSIZE:
				break
	finally:
		fin.close()
	return False

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


def is_bin(filename):
	pass

def main():
	#if is_BAM(sys.argv[1]):
	#	print "bam"
	#else:
	#	print "sam"
	
	if is_binary(sys.argv[1]): #is binary, BAM file	
		#print "This is a bam"
		outname = sys.argv[2] + ".sorted"
		pysam.sort(sys.argv[1],outname)
	else:#SAM file
		#print "This is a sam"
		try:
			inputfile = pysam.Samfile(sys.argv[1])
			
		except IOError,message:
			print >> sys.stderr, "cannot open SAM file",message
			sys.exit(1)
	#print "finish to read file"
		bamout = sys.argv[2] + ".bam"
		outname = sys.argv[2] + ".sorted"
		outputfile = pysam.Samfile(bamout,'wb',template=inputfile)
		for item in inputfile.fetch():
			outputfile.write(item)
		pysam.sort(bamout,outname)
	

if __name__=="__main__":
	main()
