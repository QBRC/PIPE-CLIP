#!/usr/bin/python
# programmer : bbc
# usage:


import sys
import re
import random
import string
import pysam
from pysam import *
import argparse as ap
from pybedtools import BedTool
import copy
import rpy2.robjects as robject
from rpy2.robjects.packages import importr
from rpy2.robjects import FloatVector
import clusterEnrich_nb as ce
import math
from collections import Counter

stats = importr('stats')

def prepare_argparser():
	description = "Looking for reliable mutations"
	epilog = "For command line options of each command, type %(prog)s COMMAND -h"
	argparser = ap.ArgumentParser(description = description, epilog = epilog)
	argparser.add_argument("-a","--bam",dest = "saminput",type = str, required = True, help = "input mapping result BAM file")
	argparser.add_argument("-b","--bed",dest = "bedinput",type = str, required = True, help = "input mutation BED file")
	argparser.add_argument("-o","--output",dest = "output",type = str, required = True, help = "output file prefix")
	argparser.add_argument("-p","--par",dest = "par",type = int, help = "CLIP type, 0 for HITS-CLIP, 1 for PAR-CLIP and iCLIP", default = 0)
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
	'''
	Get the mutation locations. Strand sensitive.
	
	'''
	muDic = {}
	muList = []
	for m in mufile:
		name  = m[0]+"\t"+m[1]+"\t"+m[2]+"\t"+m[5]+"\t"+m[6] #chr,start,stop,strand,mutationType
		if muDic.has_key(name):
			muDic[name]+=1
		else:
			muDic[name]=1
	for k in muDic.keys():
		#print k,muDic[k]
		buf = k.split("\t")
		buf.insert(3,muDic[k])
		muList.append(buf)#chr,start,stop,copyNumber,strand,mutationType
	return muList

def KMvalue(mapfile,mufile):
	'''
	Calculate K(coverage),M(mutation count) value for each mutation location
	Output in BED-like format. Score column(4th column) is the copy number of that mutation, which is the value of M.(strand specific)
	'''
	km = []
	km_pair = {}#Dic of count tuples of (k,m),key:"KkMm"
	mu_merge = mutationUniq(mufile)
	count  = 1
	for item in mu_merge:
		#print item
		st = []
		mutation = item[-1][-1] #record the mutation type
		strand = item[-2] #record the strand, 10/11 added 
		M = item[-3]
		K = 0
		for pileupColumn in mapfile.pileup(item[0],int(item[1]),int(item[2])):
			if pileupColumn.pos == int(item[1]): #find the mutation site
				K = 0 #pileupColumn.n #edited 1023
				for pileupRead in pileupColumn.pileups:
					if pileupRead.alignment.is_reverse:
						if strand == "-":
							K += 1
					else: #pileup alignment is on plus strand
						if strand == "+": #changed - into +
							K += 1
		if K>=M:
			for i in item:
				st.append(str(i))
			mu_name="mutation_"+str(count)
			count += 1
			st.insert(3,mu_name)
			st.append(str(K))
			st.append(str(M))
			#print st #chr, start, stop, name, copyNumber, strand, type, K, M
			pair_name = str(K)+"_"+str(M)
			if km_pair.has_key(pair_name):
				km_pair[pair_name] += 1
			else:
				km_pair[pair_name] = 1
			km.append(st)
	#print km #test line
	return (km,km_pair)


def uniq(b): #b is a list
	uniqElements = []
	for i in b:
		if uniqElements.count(i)==0:
			uniqElements.append(i)
	uniqElements.sort()
	return uniqElements


def muEvaluate(mapfile,mufile,cover,threshold):
	(original_KM,KM_test) = KMvalue(mapfile,mufile)
	R = robject.r
	reliableList = []
	P = len(mufile)/(cover*1.0)
	km_p = {}#store km and corresponding p value
	pvalues = []
	for k in KM_test:
		parameters = k.split("_")
		p = R.pbinom(int(parameters[1])-1,int(parameters[0]),P,False,True)[0]	
		pvalues.append(p)
		km_p[k]=p
	pCount = dict(Counter(pvalues))
	#print pCount
	pRank = ce.freqRank(pCount,True)
	total_test = len(mufile)
	pqDic={}
	for i in pRank.keys():
		try:
			p_rank = pRank[i]
			q_value = ce.BH(i,p_rank,total_test)
			pqDic[i]=q_value
		except:
			print >> sys.stderr,"Cannot find p value in dictionary"
			continue
	for record in original_KM:
		name = str(record[-2])+"_"+str(record[-1])
		record_p = km_p[name]
		record_q = pqDic[record_p]
		if record_q >= math.log(threshold)*(-1):
			record[4]=str(record_q)
			reliableList.append(record)
	return reliableList
	

def muEvaluate_check(kmfile,cover,threshold):
	original_KM = open(kmfile,'r')
	R = robject.r
	reliableList = []
	MList = []
	for row in original_KM:
		item = row.rstrip().split("\t")
		#print item
		MList.append(int(item[-1]))
	P = (sum(MList)*1.0)/cover
	#print "P is",P
	cdf = []
	pvalues = []
	original_KM.seek(0)
	for row in original_KM:
		entry = row.rstrip().split("\t")
		p = R.pbinom(int(entry[-1])-1,int(entry[-2]),P,False,True)[0]	
		pvalues.append(p)
		cdf.append((entry,p))
	
	p_adjust = ce.BH(pvalues)
	for i in cdf_adjust:
		if i[-1] >= 0:#math.log(threshold)*(-1):
			i[0][4]=str(i[-1])
			reliableList.append(i[0])
	return reliableList

def main():
	#for item in muEvaluate_check(sys.argv[1],115065019,0.01):
	#	print item
	
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
