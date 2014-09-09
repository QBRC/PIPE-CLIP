#!/usr/bin/python
# Programmer : beibei.chen@utsouthwestern.edu
# Usage:Report the truncation site of iCLIP BAM file 
# Input: BAM
# Output: BED
# Last modified: 9 Sep 2014

import sys
import re
import string
import pysam

class findTruncation:
  def __init__(self,bamFile,outputFile):
    self.mapped_bam = bamFile
    self.outputFile = outputFile

  def run(self):
    self.outputFile.write("#chr\tstart\tend\tcluster_name\toffset\tstrand\ttype");
    for item in self.mapped_bam:
      if item.is_reverse:
        strand = '-'
        start = item.pos + item.alen
      else:
        strand = '+'
        start = item.pos
      end = start + 1
      if item.tid>=0:
        chr = self.mapped_bam.getrname(item.tid)
      else:
        continue
      st = chr +"\t"+str(item.pos) + "\t"+ str(end) +"\t" + item.qname +"\t0\t" + strand+"\ttruncation"
      self.outputFile.write(st)

def findTruncationMain(mapped_bam_path,outputPath):
  try:
    mapped_bam = pysam.Samfile(mapped_bam_path,"rb")
  except IOError,message:
    print >> sys.stderr,"Cannot open BAM file.",message
    sys.exit(1)
  outputFile = open(outputPath,"w+")
  findTruncationRunner = findTruncation(mapped_bam,outputFile)
  findTruncationRunner.run()

def findTruncationMainNoArgs():
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
  findTruncationMainNoArgs()
