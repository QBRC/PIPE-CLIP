#!/usr/bin/python
# Programmer : beibei.chen@utsouthwestern.edu
# Usage:pre-process BAM file by user specified paramters
# Input:Sorted BAM file (index file in the same folder) and parameters
# Output: Filtered BAM file and two statistic figures
# Last modification: 20 Feb 2014


import sys
import re
import random
import string
import pysam
from pysam import *
import argparse as ap
import matplotlib
matplotlib.use("PDF")
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


#############Begin functions###########

def prepare_argparser():
  description = "Filter SAM file by tags"
  epilog = "For command line options of each command, type %(prog)s COMMAND -h"
  argparser = ap.ArgumentParser(description=description, epilog = epilog)
  
  #IO arguments
  argparser.add_argument("-i","--input",dest = "infile", type = str, required = True, help = "input bam file")
  argparser.add_argument("-o","--output",dest = "outfile", type = str,required = True, help = "output file, default is stdout")
  argparser.add_argument("-t","--type",dest="clipType",type = int, required = True, help = "CLIP type, 0:HITS-CLIP, 1:PAR-CLIP(4SG),2:PAR_CLIP(6SG),3:iCLIP")
  
  #CIGAR tag arguments
  argparser.add_argument("-m",dest = "matchLength", type = int, help = "Nucleotide number that are mapped to the genome")
  
  #Other tag arguments
  argparser.add_argument("-n",dest = "mismatch", type = int, default = 2, help = "NM tag number, default=2")
  
  #remove duplication arguments
  argparser.add_argument("-r",dest = "rm_loc", type = int, default = 1, help = "Remove redundant reads with the same mapping starting location ")
  
  return(argparser)
class SAMFILTERRunner:
  def __init__(self,inputfile,outputfile,coveragefile,matchLength,mismatch,rm_loc,clipType):
    self.inputfile = inputfile
    self.outputfile = outputfile
    self.coveragefile = coveragefile
    self.matchLength = matchLength
    self.mismatch = mismatch
    self.rm_loc = rm_loc
    self.clipType = clipType

  def countMatchNumber(b):
    myList = b
    m = 0
    for i in myList:
      if i[0]==1:
        m += i[1]
    return (m)

  def countMatchLength(b):
    myList = b
    m = 0
    for i in myList:
      if i[0]<=1:
        m += i[1]
    return (m)

  def countMismatch(b):
    myTag = b.tags
    for i in myTag:
      if i[0]=="NM":
        return (i[1])

  def PCRdupRm(rlist,method):
    if method == 1:
      return rmdup(rlist)
    elif method ==2:
      return rmdup_seq(rlist)

  def rmdup(duplist): #According to the same start
    mn = 0 #match number
    matchlist = []
    mq = 0 #match quality
    for i in duplist:
      m = countMatchNumber(i.cigar)
      q = i.mapq
      if m>mn  or (m==mn and q>mq):
        matchlist = [i]
        mn = m
        mq = q
      elif m == mn and q==mq:
        matchlist.append(i)
    if len(matchlist)>1: #need to select one randomly
      index =  random.randint(0,len(matchlist)-1)
      return [matchlist[index]] #return one SAM entry
    else:
      return [matchlist[0]]


  def rmdup_seq(duplist):
    '''
    Consider the same start along with sequence. Reads with same start but different sequences are not considered to be PCR duplicates.
    This function simply splits the same-start list into sub-list and use rmdup to find the representative sequence.
    '''
    seqbase = {}
    seq_remain = []
    for item in duplist:
      if seqbase.has_key(item.seq):
        seqbase[item.seq].append(item)
      else:
        seqbase[item.seq] = [item]
    for seq in seqbase.values():
      seq_remain.append(rmdup(seq))
    
  def countMatchLength(b):
    myList = b
    m = 0
    for i in myList:
      if i[0]<=1:
        m += i[1]
    return (m)

  def countMismatch(b):
    myTag = b.tags
    for i in myTag:
      if i[0]=="NM":
        return (i[1])

  def PCRdupRm(rlist,method):
    if method == 1:
      return [rmdup(rlist)]
    elif method ==2:
      return rmdup_seq(rlist)

  def rmdup(duplist): #According to the same start
    mn = 0 #match number
    matchlist = []
    mq = 0 #match quality
    for i in duplist:
      m = countMatchNumber(i.cigar)
      q = i.mapq
      if m>mn  or (m==mn and q>mq):
        matchlist = [i]
        mn = m
        mq = q
      elif m == mn and q==mq:
        matchlist.append(i)
    if len(matchlist)>1: #need to select one randomly
      index =  random.randint(0,len(matchlist)-1)
      return matchlist[index] #return one SAM entry
    else:
      return matchlist[0]


  def rmdup_seq(duplist):
    '''
    Consider the same start along with sequence. Reads with same start but different sequences are not considered to be PCR duplicates.
    This function simply splits the same-start list into sub-list and use rmdup to find the representative sequence.
    '''
    seqbase = {}
    seq_remain = []
    for item in duplist:
      if seqbase.has_key(item.seq):
        seqbase[item.seq].append(item)
      else:
        seqbase[item.seq] = [item]
    for seq in seqbase.values():
      seq_remain.append(rmdup(seq))
    return seq_remain # A list of SAM entries

  def run():
    inputfile = self.inputfile
    outputfle = self.outputfile
    coveragefile = self.coveragefile
    matchLength = self.matchLength
    mismatch = selfmismatch
    rm_loc = self.rm_loc
    clipType = self.clipType

    lenList = []
    filterFile = []
    barcount = []
    barcount.append(inputfile.mapped)
    if mismatch or matchLength: #At least some CIGAR fileter if applied
      for item in inputfile:
        b= item.cigar
        flag = 0
        if b:
          if matchLength:
            if countMatchLength(b)>=matchLength:
              flag = 1
          if mismatch:
            if countMismatch(item)<=mismatch:
              flag = flag and 1
            else:
              flag = flag and 0
          if flag:
            lenList.append(countMatchLength(b))
            filterFile.append(item)
    else:
      filterFile = inputfile
    barcount.append(len(lenList)) # For the bar chart
    if clipType==3: # if iCLIP, do not remove duplicates
      rm_loc = 0
    if rm_loc > 0: #user require to remove duplicates
      former = []
      counter = 0
      reductantList = []
      checkNewStart = 1
      lenList = []  # re-make lengList for histogram
      for item in filterFile:
        if checkNewStart:
          former = item
          reductantList = [item]
          checkNewStart = 0
          counter += 1
          continue
        if item.tid == former.tid and item.pos == former.pos and (item.is_reverse and former.is_reverse):
          reductantList.append(item)
          counter += 1
        else: #not same with existing starts
          if len(reductantList)==1: #no reductance
            outputfile.write(former)
            lenList.append(countMatchLength(former.cigar))
          else:#choose removal method
            remain = PCRdupRm(reductantList,rm_loc)
            for r in remain:
              lenList.append(countMatchLength(r.cigar))
              outputfile.write(r)
            
          former = item
          reductantList=[item]
          counter += 1
        if counter >= len(filterFile): #reach the last line of file
          if len(reductantList)<=1: #no reductance
            lenList.append(countMatchLength(item.cigar))
            outputfile.write(item)
          else:
            remain = PCRdupRm(reductantList,rm_loc)
            for r in remain:
              lenList.append(countMatchLength(r.cigar))
              outputfile.write(r)
    else:
      if len(lenList)==0:
        for item in filterFile:
          lenList.append(countMatchLength(item.cigar))
          outputfile.write(item)
      else:
        for item in filterFile:
          outputfile.write(item)
    outputfile.close()
    pysam.index(outname)
    barcount.append(len(lenList))
    print >> coveragefile, sum(lenList)


#generate remained reads mapped length distribution
  name1 = "Length_Distribution.pdf"
  pp = PdfPages(name1)
  fig = plt.figure(frameon = True)
  ax = fig.add_subplot(111)
  ax.hist(lenList,50,facecolor='green',alpha=0.5)
  ax.set_xlabel('Matched length')
  ax.set_ylabel('Frequency')
  ax.set_xlim(matchLength-1,max(lenList)+1)
  plt.savefig(pp,format='pdf')
  pp.close()

#generate barchar of reads number remained after each filter
  name2 = "Filter_Statistics.pdf"
  bp = PdfPages(name2)
  bfig = plt.figure(frameon = True)
  bax = bfig.add_subplot(111)
  bax.bar([0.9,1.9,2.9],barcount,0.2,color='green')
  bax.set_ylabel('Frequency')
  bax.set_xlim(0.5,3.5)
  bax.set_xticks([1,2,3])
  bax.set_xticklabels(['Mapped','Length&mismatch','rmdup'])
  for i in range(3):
    bax.annotate(barcount[i],(i+1,barcount[i]+0.02*max(barcount)),va='bottom',ha='center')
  plt.savefig(bp,format='pdf')
  bp.close()



def SAMFILTERRMain(infile,outfile,matchLength,mismatch,rm_loc,clipType):
  try:
    inputfile = pysam.Samfile(args.infile,"rb")
  except IOError,message:
    print >> sys.stderr, "cannot open SAM file",message
    sys.exit(1)
  outname = args.outfile+".bam"
  outputfile = pysam.Samfile(outname,'wb',template=inputfile)
  coverageName = args.outfile + ".coverage"
  coveragefile = open(coverageName,"wa")
  SAMFILTERRunner = SAMFILTERRunner(inputfile,outputfile,coveragefile,matchLength,mismatch,rm_loc,clipType)
  SAMFILTERRunner.run()

if __name__=="__main__":
  SAMFILTERMain()
