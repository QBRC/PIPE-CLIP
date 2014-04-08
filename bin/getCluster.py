#!/usr/bin/python
# programmer : bbc
# usage:

import sys
import re
import random
import string

def main():
	try:
		infile = open(sys.argv[1],"r+")
	except IOError,message:
		print >> sys.stderr, "cannot open file",message
		sys.exit(1)
	
	try:
		nbfile = open(sys.argv[2],"r+")
	except IOError,message:
		print >> sys.stderr, "cannot open file",message
		sys.exit(1)
	nbDic ={}
	for item in nbfile:
		buf = item.rstrip().split("\t")
		nb_key = "_".join(buf[0:2])
		if not nbDic.has_key(nb_key):
			nbDic[nb_key]=buf[-1]
	for record in infile:
		bufr = record.rstrip().split("\t")
		r_key = bufr[4]+"_"+str(int(bufr[2])-int(bufr[1]))
		if nbDic.has_key(r_key):
			bufr.append(nbDic[r_key])
			print "\t".join(bufr)
if __name__=="__main__":
	main()
