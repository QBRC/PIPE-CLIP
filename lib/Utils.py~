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

def readQuafilter(read,mlen,mis):
	matchlen = 0
	mismatch = 0
	for i in read.cigar:
		if i[0]<=1:
			matchlen += 1
	for j in read.tags:
		if j[0]=="NM":
			mismatch = j[1]
			break
	if matchlen>=mlen and mismatch<=mis:
		return True
	else:
		return False

