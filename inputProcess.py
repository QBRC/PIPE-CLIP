#!/usr/bin/python
# Programmer : beibei.chen@utsouthwestern.edu
# Usage: Convert SAM to BAM, and sort it
# Input: SAM/BAM
# Output: sorted BAM
# Last mofidication: 5 Mar. 2014


import sys
import re
import string
import gzip
import pysam
from pysam import *
import argparse as ap
import mimetypes

class inputProcessRunner:
  def __init__(self,inputFilePath,outputRoot):
    self.inputFilePath = inputFilePath
    self.outputFileRoot = outputRoot

  def is_BAM(self):
    try:
      self.infile  = gzip.open(self.inputFilePath)
    except IOError: # not a gzip file
      self.infile = open(self.inputFilePath)
    magic = self.infile.read(3)
    if magic == "BAM":
      return True
    else:
      return False


  def good_header(self):
    try:#Test file integrity
      header = pysam.view("-H",self.inputFilePath)
      content = pysam.view(self.inputFilePath)
      outFile = open(self.outputFileRoot+".header","w+")
      outFile.write(''.join(header))
      outFile.write(''.join(content))
      outFile.close()
      return True
    except Exception as e:
      print str(e)
      print >> sys.stderr, "Cannot read binary header, please check BAM file.)"
      return False

  def run(self):
    if self.good_header():
      if self.is_BAM():
        outname = self.outputFileRoot + ".sorted"
        bamout = self.inputFilePath
      else: #This is a SAM file, convert into BAM first
        try:
          inputfile = pysam.Samfile(self.inputFilePath)
        except IOError,message:
          print >> sys.stderr, "cannot open SAM file",message
          sys.exit(1)
        bamout = self.outputFileRoot + ".bam"
        outname = self.outputFileRoot + ".sorted"
        outputfile = pysam.Samfile(bamout,'wb',template=inputfile)
        for item in inputfile.fetch():
          outputfile.write(item)
      #Sort BAM and get index, write to file
      pysam.sort(bamout,outname)
      pysam.index(outname+".bam")
    else: #There is something wrong with the file itself
      print >> sys.stderr,"File corrupted, please check your file."
      sys.exit(1)

def inputProcessMain(inputFilePath,outputRoot):
  print inputFilePath
  ainputProcessRunner = inputProcessRunner(inputFilePath,outputRoot)
  ainputProcessRunner.run()
  
def inputProcessMainNoArgs():
  myinputProcessRunner = inputProcessRunner(sys.argv[1],sys.argv[2])
  myinputProcessRunner.run()
    
if __name__=="__main__":
  inputProcessMainNoArgs()
