#!/usr/bin/python
# Programmer : beibei.chen@utsouthwestern.edu
# Usage: Convert SAM to BAM, and sort it
# Input: SAM/BAM
# Output: sorted BAM
# Last mofidication: 25 Feb. 2014


import sys
import gzip
import pysam
from pysam import *

def is_BAM(inputFilePath):
  try:
    infile  = gzip.open(inputFilePath)
    print "done!"
  except IOError: # not a gzip file
    infile = open(inputFilePath)
    print "error"
  magic = infile.read(3)
  if magic == "BAM":
    return True
  else:
    return False
def good_header(inputFilePath):
  outputFileRoot = "apple"
  try:#Test file integrity
    #pysam.view(inputFilePath,"-H","-o"+outputFileRoot+".header",False)
    out = pysam.view(inputFilePath)
    print out
    return True
  except ValueError as e:
    print str(e)
    return False
if __name__=="__main__":
  print good_header(sys.argv[1])
