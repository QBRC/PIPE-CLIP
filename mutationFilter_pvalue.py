#!/usr/bin/python
# programmer : bbc
# usage:


import sys
import re
import random
import string
import pysam
from pysam import *
import ghmm
import argparse as ap
from pybedtools import BedTool
import copy
import rpy2.robjects as robject
from rpy2.robjects.packages import importr
from rpy2.robjects import FloatVector
import clusterEnrich as ce
import math

stats = importr('stats')


def prepare_argparser():
	description = "Looking for reliable mutations"
	epilog = "For command line options of each command, type %(prog)s COMMAND -h"
	argparser = ap.ArgumentParser(description = description, epilog = epilog)
	argparser.add_argument("-a","--bam",dest = "saminput",type = str, required = True, help = "input mapping result BAM file")
	argparser.add_argument("-b","--bed",dest = "bedinput",type = str, required = True, help = "input mutation BED file")
	argparser.add_argument("-o","--output",dest = "output",type = str, required = True, help = "output file prefix")
	argparser.add_argument("-p","--par",dest = "par",type = int, help = "specify the PAR-CLIP", default = 0)
	argparser.add_argument("-f",dest = "fdr",type = float, help = "FDR cutoff", default = 0.001)
	argparser.add_argument("-c",dest="coveragefile",type = str,required = True, help = "The sum of mapped reads length")

	return(argparser)
		
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

def mutationUniq(mufile):
	muDic = {}
	muList = []
	for m in mufile:
		name  = m[0]+"\t"+m[1]+"\t"+m[2]+"\t"+m[5]+"\t"+m[6]
		if muDic.has_key(name):
			muDic[name]+=1
		else:
			muDic[name]=1
	for k in muDic.keys():
		buf = k.split("\t")
		buf.insert(3,muDic[k])
		muList.append(buf)#chr,start,stop,offset,mutationType
	return muList

def KMvalue(mapfile,mufile):
	km = []
	mu_merge = mutationUniq(mufile)
	count  = 1
	for item in mu_merge:
		st = []
		mutation = item[-1][-1] #record the mutation type
		strand = item[-2] #record the strand, 10/11 added 
		M = item[-3]
		for pileupColumn in mapfile.pileup(item[0],int(item[1]),int(item[2])):
			if pileupColumn.pos == int(item[1]): #find the mutation site
				K = 0#pileupColumn.n edited 1023
				for pileupRead in pileupColumn.pileups:
					if pileupRead.alignment.is_reverse:
						if strand == "-":
							K += 1
					else: #pileup alignment is on plus strand
						if strand == "+": #changed - into +
							K += 1
		for i in item:
			st.append(str(i))
		mu_name="mutation_"+str(count)
		count += 1
		st.insert(3,mu_name)
		st.append(str(K))
		st.append(str(M))
		km.append(st)
	return km


def uniq(b): #b is a list
	uniqElements = []
	for i in b:
		if uniqElements.count(i)==0:
			uniqElements.append(i)
	uniqElements.sort()
	return uniqElements


def muEvaluate(mapfile,mufile,cover,threshold):
	original_KM = KMvalue(mapfile,mufile)
	R = robject.r
	reliableList = []
	MList = []
	for item in original_KM:
		MList.append(int(item[-1])*int(item[4]))
	P = (sum(MList)*1.0)/cover
	cdf = []
	for entry in original_KM:
		p = R.pbinom(int(entry[-1])-1,int(entry[-2]),P,False,True)[0]	
		cdf.append((entry,p))
	
	cdf_adjust = ce.BH(cdf)
	for i in cdf_adjust:
		if i[-1] >= math.log(threshold)*(-1):
			i[0][4]=str(i[-1])
			reliableList.append(i[0])
	return reliableList
	


def main():
	argparser = prepare_argparser()
	args = argparser.parse_args()
	
	try:
		bamfile = pysam.Samfile(args.saminput,"rb")
	except IOError,message:
		print >> sys.stderr, "cannot open mapping BAM file",message
		sys.exit(1)

	try:
		mutationfile = BedTool(args.bedinput)
	except IOError,message:
		print >> sys.stderr, "cannot open mutation BED file",message
		sys.exit(1)
	
	try:
		coverageFile = open(args.coveragefile,"r")
	except IOError,message:
		print >> sys.stderr, "cannot open coverage file",message
		sys.exit(1)
	coverage = int(coverageFile.readline().rstrip())
	if args.par > 0: #input is a par,no need to split the file
		filename = args.output+".bed"
		outputfile = open(filename,"wa")
		print >> outputfile,"#chr\tstart\tend\tname\t-log(q)\tstrand\ttype\tk\tm"
		for reliable_mu in muEvaluate(bamfile,mutationfile,coverage,args.fdr):
			print >>outputfile,'\t'.join(reliable_mu)
	else: #splitfile to insertion, deletion, substitution
		insertion = []
		deletion = []
		substitution = []
		for item in mutationfile:
			if item[-1].count("Deletion")>0:
				deletion.append(item)
			elif item[-1].count("Insertion")>0:
				insertion.append(item)
			else:
				substitution.append(item)
		del_name = args.output+"_deletion.bed"
		ins_name = args.output+"_insertion.bed"
		sub_name = args.output+"_substitution.bed"

		outfile_del = open(del_name,"wa")
		outfile_ins = open(ins_name,"wa")
		outfile_sub = open(sub_name,"wa")
		print >> outfile_ins,"#chr\tstart\tend\tname\t-log(q)\tstrand\ttype\tk\tm"
		print >> outfile_del,"#chr\tstart\tend\tname\t-log(q)\tstrand\ttype\tk\tm"
		print >> outfile_sub,"#chr\tstart\tend\tname\t-log(q)\tstrand\ttype\tk\tm"
		for reliable_mu in muEvaluate(bamfile,insertion,coverage,args.fdr):
			print >> outfile_ins,'\t'.join(reliable_mu)
		for reliable_mu in muEvaluate(bamfile,deletion,coverage,args.fdr):
			print >> outfile_del,'\t'.join(reliable_mu)
		for reliable_mu in muEvaluate(bamfile,substitution,coverage,args.fdr):
			print >> outfile_sub,'\t'.join(reliable_mu)
	

if __name__=="__main__":
	main()
