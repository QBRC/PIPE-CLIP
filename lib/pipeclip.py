#Main pipeline connects all the scripts together
#Original Programmer: beibei.chen@utsouthwestern.edu
#Refactored by: eric.roos@utsouthwestern.edu
#Usage: python pipeclip.py input.sam output_prefix match_length mismatch_number pcr_rm fdr_cluster clip_type fdr_mutation species
#Required packages: pysam, ghmm, pybedtools
#Last modification: 3 March 2014

#from lib import *

import sys
import argparse
import logging
from subprocess import call
import os
import CLIP
import Alignment
import Utils
import Enrich
import OptValidator

def prepare_argparser():
	description = "Find mutations"
	epilog = "For command line options of each command, type %(prog)s COMMAND -h"
	argparser = argparse.ArgumentParser(description=description, epilog = epilog)
	argparser.add_argument("-i","--input",dest = "infile", type = str, required = True, help = "input bam file")
	argparser.add_argument("-o","--output",dest = "outfile", type = str,required = True, help = "output file, default is stdout")
	argparser.add_argument("-l","--matchLength",dest = "matchLength", type = int ,required = True, help = "shorted matched segment length")
	argparser.add_argument("-m","--mismatch",dest = "mismatch", type = int,required = True, help = "maximum mismatch number")
	argparser.add_argument("-c","--clipType",dest = "clipType", type = int,required = True, help = "CLIP type (0)HITS-CLIP; (1)PAR-4SU; (2)PAR-6SG; (3)iCLIP", choices=[0,1,2,3])
	argparser.add_argument("-r","--rmdup",dest = "dupRemove", type = int,required = True, help = "Remove PCR duplicate (0)No removal; (1)Remove by read start; (2)Remove by sequence; ", choices=[0,1,2])
	argparser.add_argument("-M","--fdrMutation",dest = "fdrMutation", type = float,required = True, help = "FDR for reliable mutations")
	argparser.add_argument("-C","--fdrCluster",dest = "fdrCluster", type = float,required = True, help = "FDR for enriched clusters")
  #argparser.add_argument("-s","--species",dest = "species", type = str,required = True, help = "Species [\"mm10\",\"hg19\"]",choices=["mm10","hg19"])
	return(argparser)

def runPipeClip(infile,outputPrefix,matchLength,mismatch,rmdup,fdrEnrichedCluster,clipType,fdrReliableMutation,species):
	myClip = CLIP.CLIP(infile)
	logging.info("Start to run")
	if myClip.testInput():#check input
		logging.info("Input file OK,start to run PIPE-CLIP")
		if myClip.readfile():
			myClip.filter(matchLength,mismatch,clipType,rmdup)
			#myClip.printFilteredReads()
			#myClip.printClusters()
			#myClip.printMutations()
			if len(myClip.clusters)>0:
				logging.info("Get enriched clusters")
				Enrich.clusterEnrich(myClip,fdrEnrichedCluster)
				logging.debug("Found %d enriched clusters" % myClip.sigClusterCount)
				myClip.printReliableList(myClip.clusters)
			else:
				logging.error("There is no clusters found. Please check input.Exit program.")
				sys.exit(1)
			
#BC#			if len(myClip.mutations.keys())>0:
#BC#				logging.info("Get reliable mutations")
#BC#				Enrich.mutationEnrich(myClip,fdrReliableMutation)
#BC#				#myClip.printEnrichedItem(myClip.sigMutations)
#BC#			else:
#BC#				logging.warning("There is no mutation found in this BAM file.")
#BC#			if myClip.sigClusterCount > 0 and myClip.sigMutationCount>0:
#BC#				logging.info("Get cross-linking sites")
#BC#				myClip.getCrosslinking()
#BC#				if (len(myClip.crosslinking.keys())>0):
#BC#					myClip.printCrosslinkingSites()
#BC#					myClip.printList(myClip.crosslinkingMutations)
#BC#				else:
#BC#					logging.warning("There is no crosslinking found. May be caused by no reliable mutations in enriched clusters. Print out enriched clusters instead.")
#BC#					myClip.printEnrichClusters()
#BC#			else:
#BC#				if myClip.sigClusterCount <= 0:
#BC#					logging.error("There is no enriched clusters for this sample, please check your input file. Exit.")
#BC#					sys.exit(2)
#BC#				elif myClip.sigMutationCount <=0:
#BC#					logging.warning("There is no reliable mutations found. PIPE-CLIP will provide enriched clusters as crosslinking candidates.")
#BC#					myClip.printEnrichClusters()
	else:
		print >> sys.stderr, "File corruputed, program exit."
		sys.exit(0)
	
	

if __name__=="__main__":
	arg_parser = prepare_argparser()
	args = arg_parser.parse_args()
	OptValidator.opt_validate()
	infile = args.infile                  # Input SAM/BAM file
	outputPrefix = args.outfile           # Output prefix
	matchLength = args.matchLength        # Shorted matched segment length
	mismatch = args.mismatch              # Maximum mismatch number
	rmcode = args.dupRemove
	fdrEnrichedCluster = args.fdrCluster  # FDR for enriched clusters
	clipType =args.clipType               # CLIP type (0)HITS-CLIP; (1)PAR-4SU; (2)PAR-6SG; (3)iCLIP
	fdrReliableMutation = args.fdrMutation# FDR for reliable mutations
#BC#  #species = args.species                # Species ["mm10","hg19"]
	runPipeClip(infile,outputPrefix,matchLength,mismatch,rmcode,fdrEnrichedCluster,clipType,fdrReliableMutation,None)
