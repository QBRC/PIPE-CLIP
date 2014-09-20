#Main pipeline connects all the scripts together
#Original Programmer: beibei.chen@utsouthwestern.edu
#Refactored by: eric.roos@utsouthwestern.edu
#Usage: python pipeclip.py input.sam output_prefix match_length mismatch_number pcr_rm fdr_cluster clip_type fdr_mutation species
#Required packages: pysam, ghmm, pybedtools
#Last modification: 18 Sep 2014


import sys
import argparse
import logging
import os
from lib import *

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
	argparser.add_argument("-s","--species",dest = "species", type = str, help = "Species [\"mm10\",\"hg19\"]",choices=["mm10","hg19"])
	return(argparser)

def runPipeClip(infile,outputPrefix,matchLength,mismatch,rmdup,fdrEnrichedCluster,clipType,fdrReliableMutation,species):
	myClip = CLIP.CLIP(infile,outputPrefix)
	logging.info("Start to run")
	if myClip.testInput():#check input
		logging.info("Input file OK,start to run PIPE-CLIP")
		logging.info("Species info %s" % species)
		if myClip.readfile():
			myClip.filter(matchLength,mismatch,clipType,rmdup)
			#myClip.printMutations()
			if len(myClip.clusters)>0:
				logging.info("Get enriched clusters")
				status = Enrich.clusterEnrich(myClip,fdrEnrichedCluster)
				if status:
					logging.info("Found %d enriched clusters" % myClip.sigClusterCount)
					myClip.printEnrichedClusters()
				else:
					logging.error("There is no enriched cluster found. Exit program")
					sys.exit(1)
			else:
				logging.error("There is no clusters found. Please check input.Exit program.")
				sys.exit(1)
			
			if len(myClip.mutations.keys())>0:
				logging.info("Get reliable mutations")
				Enrich.mutationEnrich(myClip,fdrReliableMutation)
				myClip.printReliableMutations()
			else:
				logging.warning("There is no mutation found in this BAM file.")
			#Start to get crosslinking sites
			if myClip.sigClusterCount > 0 and myClip.sigMutationCount>0:
				logging.info("Get cross-linking sites")
				myClip.getCrosslinking()
				if (len(myClip.crosslinking.keys())>0):
					outfilelist = myClip.printCrosslinkingSites()
					myClip.printCrosslinkingMutations()
				else:
					logging.warning("There is no crosslinking found. May be caused by no reliable mutations in enriched clusters. Print out enriched clusters instead.")
					outfilelist = myClip.printEnrichClusters()
			else:
				if myClip.sigClusterCount <= 0:
					logging.error("There is no enriched clusters for this sample, please check your input file. Exit.")
					sys.exit(2)
				elif myClip.sigMutationCount <=0:
					logging.warning("There is no reliable mutations found. PIPE-CLIP will provide enriched clusters as crosslinking candidates.")
					outfilelist = myClip.printEnrichClusters()
			#annotation if possible
		if species in ["mm10","mm9","hg19"]:
			for name in outfilelist:
				#logging.debug("Start to do annotation for %s" % name)
				Utils.annotation(name,species)
		#output a status log file
		logfile = open(outputPrefix+".pipeclip.summary.log","w")
		print >> logfile, "There are %d mapped reads in input BAM file. After filtering,%d reads left" % (myClip.originalMapped,myClip.filteredAlignment)
		print >> logfile, "%d out of %d clusters are enriched." % (self.sigClusterCount,len(self.clusters))
		print >> logfile, "%d out of %d mutations are reliable." % (self.sigMutationCount,self.mutationCount)
		print >> logfile, "%d crosslinking site candidates are found, with %d supporting reliable mutations." % (len(myClip.crosslinking.keys()),len(myClip.crosslinkingMutations))
		logfile.close()
		
	else:
		logging.error("File corruputed, program exit.")
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
	species = args.species                # Species ["mm10","hg19"]
	runPipeClip(infile,outputPrefix,matchLength,mismatch,rmcode,fdrEnrichedCluster,clipType,fdrReliableMutation,species)
