#!/usr/bin/python
# programmer : bbc
# usage:
# output in BED format, the score is the offset of the mutation from 5' end

import sys
import re
import random
import string
import copy
import pysam
from pysam import *
import argparse as ap
import Alignment


class MutationBed(Alignment.BED):
	def __init__(self,chr,start,stop,name,score,strand,type):
		BED(chr,start,stop,name,score,strand,type)
		self.type = type #insertion,deletion,type of substitution

def countMatchNumber(b):
	myList = b
	m = 0
	for i in myList:
		if i[0]==0:
			m += i[1]
	return (m)

def countInsertionNumber(b):
	myList = b
	m = 0
	for i in myList:
		if i[0]==1:
			m += i[1]
	return (m)

def countDeletionNumber(b):
	myList = b
	m = 0
	for i in myList:
		if i[0]==2:
			m += i[1]
	return (m)

def countMismatch(b):
	myList = b
	for i in myList:
		if i[0]=="NM":
			return (i[1])

def survey(entry):
	mismatchNumber = countMismatch(entry.tags)
	insertionNumber = countInsertionNumber(entry.cigar)
	deletionNumber = countDeletionNumber(entry.cigar)
	substitutionNumber = mismatchNumber-insertionNumber-deletionNumber
	#print insertionNumber,deletionNumber
	return ([insertionNumber,deletionNumber,substitutionNumber])


def SBeforeFirstM(ci):
	for index in range(len(ci)):
		if ci[index][0] == 0:
			if index == 0:
				return 0
			else:
				s = 0
				for i in range(0,index):
					if ci[i][0] == 4:
						s += ci[i][1]
				return s
def parseCIGAR(ci): #calculate match segment length list
	matchSeg = []
	preIndex = 0
	for index in range(len(ci)):
		if ci[index][0] == 0: #match tag
			preIndex = index
			matchSeg.append(ci[index][1])
			if index > 0: #M is not the first tag
				s = 0
				for i in range(preIndex-1,index): #check all the tags before M
					if ci[i][0]==4:
						s += ci[i][1]
				matchSeg[-1] += s 
			#seach for the insertion behind (before the next match)
			for insertionIndex in range(index+1,len(ci)):
				inse = 0
				if ci[insertionIndex][0]==0:
					break
				elif ci[insertionIndex][0]==1:#insertion exists
					inse += ci[insertionIndex][1]
			matchSeg[-1]+=inse
	return matchSeg	

def parseMD(b):
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


def insertionLocation(entry,num):#sam entry and insertion total number
	insertionLoc = {} #key:ref_offset; value:list of seq_offset
	ref_offset = 0
	p
	seq_offset = [0]
	preInsertion = 0 #used to check if it is the last insertion tag
	for i in entry.cigar:
		if i[0]==0 or i[0]==4:#match or soft clip
			ref_offset += i[1]
			seq_offset[-1] += i[1]
		elif i[0]== 2: #deletion counts on ref but not on read seq
			ref_offset+=i[1]
		elif i[0]==1:#record this insertion tag and prepare for the next one
			preInsertion += i[1]
			for ins in range(1,i[1]):
				newloc = seq_offset[-1]+1
				seq_offset.append(newloc)
			st = []
			for item in seq_offset:
				st.append(item)
			insertionLoc[ref_offset]=st
			if preInsertion == num: #this is the last insertion tag
				return insertionLoc
			else:#insertion only count on read seq
				seq_index = seq_offset[0]+i[1]

def countInsertionBefore(seqLoc,insertLocList):
	if len(insertLocList)==0:
		return 0
	else:
		ins = 0
		for i in insertLocList:
			if i>seqLoc:
				return ins
			else:
				ins += 1
		return ins

def  mutationLocation(entry,insertLoc):#return mutation location in 
	mutations = []
	match = entry
	myList = match.tags
	S_count = SBeforeFirstM(match.cigar)#hard clip will cut the seq, doesn't count
	mis = countMismatch(myList)
	mdlist = parseMD(myList)
	mdtag = ''.join(mdlist)
	counter = 0
	if not mdtag.isdigit(): #insertions may exist
		st_seq = 0
		st_genome = 0
		offset = 0
		pre = ':'
		for ch in mdlist:#i[1] is the MD tag
			if ch.isdigit():
				st_seq += int(ch)
				st_genome += int(ch)
				pre =  ch
			elif ch == "^":
				st_genome += 1
				pre = ch
			elif ch.isalpha():
				if not pre == '^':
					origin = ch
					index = st_seq+S_count+offset
					insertionBefore = countInsertionBefore(index,insertLoc)
					loc = st_genome+match.pos+offset#-insertionBefore #0-based 
					index += insertionBefore # add on 9 Oct
					#print index,len(match.seq)
					mu = match.seq[index]
					offset = index-S_count+1
					st_seq = 0
					st_genome = 0
					if match.is_reverse:
						chr = '-'
						t = RC([origin,mu])
						origin = t[0]
						mu = t[1]
					else:
						chr = '+'
					mutation = [str(loc),str(loc+1),match.qname,str(index-S_count),chr,origin+"->"+mu]
					yield mutation
				else:
					loc = st_genome+match.pos+offset-1 #0-based 
					if match.is_reverse:
						strand = '-'
						ch = RC([ch])[0]
					else:
						strand = '+'
					index1 = loc - match.pos 
					insertionBefore = countInsertionBefore(index1,insertLoc)
					index1 += insertionBefore #added 9 Oct
					mutation = [str(loc),str(loc+1),match.qname,str(index1),strand,"Deletion->"+ch]
					yield mutation
					pre = ch

	return

def RC(strList):
	rc = []	
	for item in strList:
		st = ""
		for chIndex in range(len(item)):
			rcIndex = len(item)-1
			if item[rcIndex].upper()== "A":
				st += 'T'
			elif item[rcIndex].upper()=="C":
				st += 'G'
			elif item[rcIndex].upper()=="T":
				st += 'A'
			elif item[rcIndex].upper()=="G":
				st += 'C'
			else:
				st += 'N'
		rc.append(st)
	return(rc)

def getMutaions(infile,read):
	#tmp = pysam.Samfile("demo.bam",'wb',template=infile)
	mutationList = []
	b= item.tags
	if countMismatch(b)>0: #and countMismatch(b)<2 and countMatchNumber(item.cigar)>=20:
		#tmp.write(item)
		sur = survey(item)
		insertion = sur[0]
		deletion = sur[1]
		substi = sur[2]
		insertionSeqLoc = []
		if insertion > 0:
			insertionDic = insertionLocation(item,insertion)
			for k in insertionDic.keys():
				for loc_index in range(len(insertionDic[k])):
					insertionSeqLoc.append(insertionDic[k][loc_index])
					mu = item.seq[insertionDic[k][loc_index]]
					loc = k+loc_index+item.pos
					if item.tid >=0:
						chr = infile.getrname(item.tid)
					if item.is_reverse:
						strand = '-'
						mu = RC([mu])[0]
					else:
						strand = "+"
					mutationList.append(MutationBed(chr,loc,loc+1,item.qname,insertionDic[k][loc_index]),strand,"Insertion->"+mu)
			insertionSeqLoc.sort()
		if deletion + substi > 0:
			for mu in mutationLocation(item,insertionSeqLoc):
				if item.tid>=0:
					chr  = infile.getrname(item.tid)
					newMu = MutationBed(chr,mu[0],mu[1],mu[2],mu[3],mu[4],mu[5])
					mutationList.append(newMu)
		return mutationList


def getTruncation(infile,read):
	'''One read only have one truncation, but in order to be consistent with muations, return a list'''
	if read.is_reverse:
		strand = "-"
		start = read.pos + read.alen
	else:
		strand = "+"
		start = read.pos
	stop = start + 1
	if read.tid>=0:
		chr = infile.getrname(read.tid)
	newTr = MutationBed(chr,start,stop,item.qname,1,strand,"truncation")
	return [newTr] 

