#!/usr/bin/python
# programmer : bbc
# usage:

import sys
import re
import random
import string
import math
from pysam import *
from pybedtools import BedTool



def main():
	try:
		mapped_bam = BedTool(sys.argv[1])
	except IOError,message:
		print >> sys.stderr,"Cannot open BAM file.",message
		sys.exit(1)
	
	mapped_bed = mapped_bam.bam_to_bed()
	bed_merge = mapped_bed.merge(s=True, n=True)
	
	#print "#chr\tstart\tend\tcluster_name\treads_in_cluster\tstrand"
	count = 1
	for item in bed_merge:
		name = "cluster_"+str(count)
		count += 1
		print "%s\t%s\t%s\t%s\t%s\t%s" % (item[0],item[1],item[2],name,item[3],item[4])

if __name__=="__main__":
	main()
