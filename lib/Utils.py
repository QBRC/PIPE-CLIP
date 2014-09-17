#!/usr/bin/python
# programmer: beibei.chen@utsouthwestern.edu
# Usage: definition of utilities to handel pysam classes

import sys
import Alignment
import subprocess
import logging

def is_sorted(header):
	'''Get a BAM header, and check if this file is sorted or not
	   Return: True: sorted, False: unsorted
	   Header variable is a list
	'''
	for row in header:
		buf = row.rstrip().split("\t")
		if buf[0]=="@HD":
			if buf[2] in ["SO:unsorted","SO:unknown"] :
				return False
			elif buf[2] in ["SO:coordinate","SO:queryname"]:
				return True
	#If no HD header contains SO info, return False
	return False 

def readQuaFilter(read,mlen,mis):
	matchlen = 0
	mismatch = 0
	try:
		for i in read.cigar:
			if i[0]<=1:
				matchlen += i[1]
		for j in read.tags:
			if j[0]=="NM":
				mismatch = j[1]
				break
	except:
		logging.debug("Problematic read %s" % (read))
	#print >> sys.stderr,"mLen and mis for read",read.qname,"are",matchlen,mismatch
	if matchlen>=mlen and mismatch<=mis:
		return (True,matchlen,mismatch)
	else:
		return (False,0,0)

def rmdupKey_Start(read):
	#print >>sys.stderr,"Remove duplicate by alignment start"
	k = str(read.tid)+":"+str(read.pos)+":"+str(read.is_reverse)
	#print >>sys.stderr,k
	return k


def rmdupKey_Seq(read):
	k = str(read.tid)+":"+str(read.pos)+":"+str(read.is_reverse)+":"+read.seq
	return k

def filterMutations(mlist,feature,keep=True):
	newList = []
	for i in mlist:
		if i.type==feature:
			if keep:
				newList.append(i)
		else:
			if not keep:
				newList.append(i)
	return newList

def annotation(fileaddr,species):
	'''Call Homer annotation script'''
	outfile = open(fileaddr+".anno","w")
	args = ['annotatePeaks.pl',fileaddr,species]
	p = subprocess.Popen(args,stdout=outfile)
	stdout_value = p.communicate()[0]
	outfile.close()

def bisort(alist,b):
	'''List is a list of bed instances, b is also an instance'''
	if isinstance(alist[0],Alignment.BED) and isinstance(b,Alignment.BED):
		pass 
