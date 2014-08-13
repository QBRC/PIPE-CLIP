#!/usr/bin/python
# programmer: beibei.chen@utsouthwestern.edu
# Usage: definition of utilities to handel pysam classes

import sys

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
	for i in read.cigar:
		if i[0]<=1:
			matchlen += i[1]
	for j in read.tags:
		if j[0]=="NM":
			mismatch = j[1]
			break
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

def parseMD(self,b):
	myList = b
	for j in myList:
		if j[0]=="MD":
			md = copy.deepcopy(j[1])
			num = md.replace('A','^').replace('G','^').replace('T','^').replace('C','^').replace('^^','^').split('^')
			for c in range(num.count("")):
				num.remove("")
		buf = [num[0]]
		counter = 1
		afterAlpha = 0
		for i in j[1]:#walk thought MD string to record mutation and deletion
			if i.isalpha() or i == '^':
				buf.append(i)
				afterAlpha = 1
			else:
				if afterAlpha and counter <= len(num)-1:
					buf.append(num[counter])
					afterAlpha = 0
					counter += 1
		break
	return buf
