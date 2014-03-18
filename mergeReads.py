#!/usr/bin/python
# Programmer : beibei.chen@utsouthwestern.edu
# Usage: generate clusters with reads count for cluster enrichment analysis
# Input: BAM
# Output: BED
# Last modified: 17 Mar. 2014


import sys
import string
import pysam
from pybedtools import BedTool

class BED:
		def __init__(self,chr,start,stop,name,score,strand):
			self.chr = chr
			self.start = start
			self.stop = stop
			self.name = name
			self.score = score
			self.strand = strand
	
		def __str__(self):
			st = "\t".join([self.chr,str(self.start),str(self.stop),self.name,str(self.score),self.strand])
			return st

		def merge(self,read):
				self.stop = read.stop
				self.score += 1

		def overlap(self,read):
				if self.chr == read.chr or self.strand == read.strand:
						if self.start <= read.stop and self.stop >=read.start:
								return True
				else:
						return False


class mergeReadsRunner:
	def __init__(self,mapped_bam,producedFile):
		self.mapped_bam =	mapped_bam
		self.output = producedFile

	def run(self):
		cluster_counter = 1
		for index,read in enumerate(self.mapped_bam):
			#change bam to bed
			if read.tid>=0:
					chr = self.mapped_bam.getrname(read.tid)
					if read.is_reverse:
							strand = "-"
					else:
							strand = "+"
					read_bed = BED(chr,read.pos,read.pos+read.alen,"cluster",1,strand)
					if index == 0:
							self.this_group = read_bed 
							self.this_group.name = "cluster"+str(cluster_counter)
							continue
					else:#compare recent read to recent group
							if self.this_group.overlap(read_bed):
									self.this_group.merge(read_bed)
							else:#print the recent group and make a new one
									print >> self.output, self.this_group
									cluster_counter += 1
									self.this_group = read_bed
									self.this_group.name = "cluster"+str(cluster_counter)
			else: #the read is unmapped, which should not be the case here
					continue
		#print the last group
		print >> self.output, self.this_group

def mergeReadsMain(bamFilePath,producedFilePath): #Eric, please use this main funtion when you modify
	try:
		mapped_bam = pysam.Samfile(bamFilePath) #change sys.argv[1] to bamFilePath
	except IOError,message:
		print >> sys.stderr,"Cannot open BAM file.",message
		sys.exit(1)
	
	producedFile = open(producedFilePath,"w+") #change sys.argv[2] to produceFilePath
	myRunner = mergeReadsRunner(mapped_bam,producedFile)
	myRunner.run()

def mergeReadsMainNoArgs():
	try:
		mapped_bam = pysam.Samfile(sys.argv[1]) #change sys.argv[1] to bamFilePath
	except IOError,message:
		print >> sys.stderr,"Cannot open BAM file.",message
		sys.exit(1)
	
	producedFile = open(sys.argv[2],"w+") #change sys.argv[2] to produceFilePath
	myRunner = mergeReadsRunner(mapped_bam,producedFile)
	myRunner.run()
	

	"""
	mapped_bed = mapped_bam.bam_to_bed()
	bed_merge = mapped_bed.merge(s=True, n=True)
	
	count = 1
	for item in bed_merge:
		name = "cluster_"+str(count)
		count += 1
		print "%s\t%s\t%s\t%s\t%s\t%s" % (item[0],item[1],item[2],name,item[3],item[4])
	"""
if __name__=="__main__":
	mergeReadsMainNoArgs()
