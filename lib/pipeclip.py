#Main pipeline connects all the scripts together
#Original Programmer: beibei.chen@utsouthwestern.edu
#Refactored by: eric.roos@utsouthwestern.edu
#Usage: python pipeclip.py input.sam output_prefix match_length mismatch_number pcr_rm fdr_cluster clip_type fdr_mutation species
#Required packages: pysam, ghmm, pybedtools
#Last modification: 3 March 2014

#from lib import *

import sys
import argparse
from subprocess import call
import os
import CLIP
import Alignment
import Utils

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

def runPipeClip(infile):#,outputPrefix,matchLength,mismatch,pcr,fdrEnrichedCluster,clipType,fdrReliableMutation,species):
  ########################## Check input #######################
  myClip = CLIP.CLIP(infile)
  if myClip.testInput():
	  pass
  else:
	  print >> sys.stderr, "File corruputed, program exit."
	  sys.exit(0)
	
	#########################Process the input file #############
	#1: check the match-length and mismatch count
	#2: check the PCR duplication removal options
	#3: After get the results from rmdup, merge them to clusters 
	#   and recored the mutations at the same type, also count coverage
	if myClip.readfile():
		myClip.filter(args.matchLength,args.mismatch,args.clipType,args.dupRemove)

if __name__=="__main__":
  arg_parser = prepare_argparser()
  args = arg_parser.parse_args()

  infile = args.infile                  # Input SAM/BAM file
#BC#  outputPrefix = args.outfile           # Output prefix
#BC#  matchLength = args.matchLength        # Shorted matched segment length
#BC#  mismatch = args.mismatch              # Maximum mismatch number
#BC#  pcr = args.clipType                    # PCR removal: (0)no removal; (1)same-start removal; (2)same-seq removal  
#BC#  fdrEnrichedCluster = args.fdrCluster  # FDR for enriched clusters
#BC#  clipType =args.clipType               # CLIP type (0)HITS-CLIP; (1)PAR-4SU; (2)PAR-6SG; (3)iCLIP
#BC#  fdrReliableMutation = args.fdrMutation# FDR for reliable mutations
#BC#  #species = args.species                # Species ["mm10","hg19"]
  runPipeClip(infile)#,outputPrefix,matchLength,mismatch,pcr,fdrEnrichedCluster,clipType,fdrReliableMutation,None)
