#!/usr/bin/python
# Programmer : beibei.chen@utsouthwestern.edu
# Usage:Report the truncation site of iCLIP BAM file 
# Input: BAM
# Output: BED
# Last modified: 19 Dec 2013

import sys
import re
import string
import pysam

class findTruncation:
  def __init__(self,bamFile):
    self.mapped_bam = bamFile

  def run(self):
    print "#chr\tstart\tend\tcluster_name\toffset\tstrand\ttype"
	  for item in self.mapped_bam:
		  end = item.pos+1
		  if item.is_reverse:
			  strand = '-'
		  else:
			  strand = '+'
		  if item.tid>=0:
			  chr = self.mapped_bam.getrname(item.tid)
		  else:
			  continue
		  st = chr +"\t"+str(item.pos) + "\t"+ str(end) +"\t" + item.qname +"\t0\t" + strand+"\ttruncation"
		  print st

def findTruncationMain():
	try:
		mapped_bam = pysam.Samfile(sys.argv[1],"rb")
	except IOError,message:
		print >> sys.stderr,"Cannot open BAM file.",message
		sys.exit(1)

  findTruncationRunner = findTruncation(mapped_bam)
  fundTruncationRunner.run()
	"""
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
  """
if __name__=="__main__":
	findTruncationMain()
