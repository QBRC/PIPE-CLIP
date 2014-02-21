#!/usr/bin/python
# Programmer : beibei.chen@utsouthwestern.edu
# Usage: remove PCR duplicates by comparing the sequences and return fastq file with barcode removed.
# Input: fastq file
# Output: fastq file
# Last modified: Feb 19 Dec. 2014 - Eric Roos <eric.roos@utsouthwestern.edu>

import sys
import re
import string

class fastq:
	def __init__(self,x):
		self.id = x[0]
		self.seq = x[1]
		self.name = x[2]
		self.quality = x[3]
	
  #Adding a more direct constructor, instead of the paramater array version
  #this should be a standard
  def __init__(self,id,seq,name,quality):
		self.id = id
		self.seq = seq
		self.name = name
		self.quality = quality 

  def __str__(self):
		st = self.id+self.seq+self.name+self.quality
		return st

class fqList:
	def __init__(self):
		self.data = []
		self.unique = {}

	def readFq(self,fh,bl):
		st = []
		buf =  fh.readlines()
		for i in range(len(buf)):
			if buf[i][0]=="@" and re.search(" ",buf[i]):
				st.append(buf[i])
				st.append(buf[i+1])
				st.append(buf[i+2])
				st.append(buf[i+3].rstrip())
				self.addFq(st)
				if not self.unique.has_key(buf[i+1]):
					self.addUniqFq(st,bl)
				st = []

	def addFq(self,s):
		newFq = fastq(s)
		self.data.append(newFq)

	def addUniqFq(self,s,b):
		seq_key = s[1]
		offset = b
		newFq = fastq([s[0],s[1][offset:],s[2],s[3][offset:]])
		self.unique[seq_key] = newFq

class barcodeRemover:
  def __init__(self,infile,barLen):
    self.infile = infile
    self.barLen = barLen

  def run(self):
	  myfq = fqList()
	  myfq.readFq(self.infile,self.barLen)
	  for item in myfq.unique.values():
		  if len(item.seq)>=15:
			  print item

def barCodeRemovalMain():
	try:
		infile = open(sys.argv[1],"r+")
	except IOError,message:
		print >> sys.stderr, "cannot open file",message
		sys.exit(1)	
  try:
    barLen = int(sys.argv[2])
  except:
    barLen = 5
  barcodeRemover = barcodeRemover(infile,barLen)
  barcodeRemover.run()	

if __name__=="__main__":
	barCodeRemovalMain()
