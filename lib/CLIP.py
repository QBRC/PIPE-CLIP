'''Define CLIP class'''
import sys
import gzip
import logging
import pysam
import random


class CLIP:
	def __init__(self,fileaddr):
		self.filepath = fileaddr
		self.originalBAM = None
		self.filteredAlignment = []
		self.clusters = {} #Dictionary of bed instance
		self.mutations = {} #Dictionary of bed instance
		self.wig = strand
		self.coverage = 0 #"reads coverage of this sample"
		self.bamheader = None
	
	def __str__(self):
		pass

	def testInput(self):
	'''Test the input file format, modify self.filepath and return bool
		Status: True:file is ready to use; False: wrong file, program stop
	'''
		#test if file has header
	try:
			self.header = pysam.view("-H",self.filepath)
		except:
			print >> sys.stderr, "Input file does not have header, please check your file. Program quit"
			return (False,"None")
		#Header test passed, test if it is BAM
		try:
			infile = gzip.open(self.filepath)
			infile.readline()
		except:#cannot read line, should be sam
			print >> sys.stderr,"Input is SAM, converting to BAM..."
			bamout = ".".join(self.filepath.split(".")[0:-1])+"."+"bam"
			pysam.view("-bS",self.filepath,bamout)
			self.filepath = bamout
		#Now the infile is BAM,check if it is sorted
		if Utils.is_sorted(self.header):
			pysam.index(self.filepath)
			return True
		else:#sort the BAM
			print >> sys.stderr, "Input is not sorted, sorting file..."
			bamsort = ".".join(self.filepath.split(".")[0:-1])+"."+"sort"
			pysam.sort(self.filepath,bamsort)
			pysam.index(bamsort)
			self.filepath = bamsort+".bam" # change input file path
			return True
			





