#Main pipeline connects all the scripts together
#Original Programmer: beibei.chen@utsouthwestern.edu
#Refactored by: eric.roos@utsouthwestern.edu
#Usage: python pipeclip.py input.sam output_prefix match_length mismatch_number pcr_rm fdr_cluster clip_type fdr_mutation species
#Required packages: pysam, ghmm, pybedtools
#Last modification: 3 March 2014

import sys
import argparse
import annotatePeaks
import barcodeRemoval
import findMutation
import findTruncation
import getCluster
import getCrosslinking
import inputProcess
from mergeReads import *
import mutationFilter
import SAMFilter
from subprocess import call

def prepare_argparser():
  description = "Find mutations"
  epilog = "For command line options of each command, type %(prog)s COMMAND -h"
  argparser = argparse.ArgumentParser(description=description, epilog = epilog)
  argparser.add_argument("-i","--input",dest = "infile", type = str, required = True, help = "input bam file")
  argparser.add_argument("-o","--output",dest = "outfile", type = str,required = True, help = "output file, default is stdout")
  argparser.add_argument("-l","--matchLength",dest = "matchLength", type = int ,required = True, help = "shorted matched segment length")
  argparser.add_argument("-m","--mismatch",dest = "mismatch", type = int,required = True, help = "maximum mismatch number")
  argparser.add_argument("-c","--clipType",dest = "clipType", type = int,required = True, help = "CLIP type (0)HITS-CLIP; (1)PAR-4SU; (2)PAR-6SG; (3)iCLIP", choices=[0,1,2,3])
  argparser.add_argument("-M","--fdrMutation",dest = "fdrMutation", type = float,required = True, help = "FDR for reliable mutations")
  argparser.add_argument("-C","--fdrCluster",dest = "fdrCluster", type = float,required = True, help = "FDR for enriched clusters")
  argparser.add_argument("-s","--species",dest = "species", type = str,required = True, help = "Species [\"mm10\",\"hg19\"]",choices=["mm10","hg19"])
  return(argparser)

def runPipeClip(infile,outputPrefix,matchLength,mismatch,pcr,fdrEnrichedCluster,clipType,fdrReliableMutation,species):
  ########################## Process input #######################
  print "input process main"
  inputProcess.inputProcessMain(infile,outputPrefix)
  print "samfilter main"
  SAMFilter.SAMFILTERMain(outputPrefix+".sorted.bam",outputPrefix+".filter",matchLength,mismatch,pcr,clipType)

  ######################### Enrich clusters ######################
  print "mergeReadsMain"
  mergeReadsMain(outputPrefix+".filter.bam",outputPrefix+".filter.rehead.merge")

  ######################### Enrich clusters ######################
  print "R Analysis"
  call(["Rscript","ZTNB.R",outputPrefix+".filter.rehead.merge",str(fdrEnrichedCluster)])
  print "getCluster main"
  getCluster.getClusterMain(outputPrefix+".filter.rehead.merge",outputPrefix+".filter.rehead.merge.ztnb", outputPrefix+".filter.cluster.bed")

  ########################### Mutation #############################
  print "Mutation section"
  if clipType == "3":
    findTruncation.findTruncationMain(outputPrefix+".filter.bam",outputPrefix+".filter.mutation.bed")
  else:
    findMutation.findMutationMain(outputPrefix+".filter.bam",outputPrefix+".filter.mutation.bed",clipType)
  print "mutaiton filter"
  mutationFilter.mutationFilterMain(outputPrefix+".filter.bam",outputPrefix+".filter.mutation.bed",outputPrefix+".filter.reliable",clipType,fdrReliableMutation,outputPrefix+".filter.coverage")

  ######################### Merge and annotation ################
  #print "Merge and annotaiton"
  #if clipType == "0":
  #  getCrosslinking.getCrossLinkingMain(outputPrefix+".filter.cluster.bed",outputPrefix+".reliable_deletion.bed", outputPrefix+"crosslinking.deletion.bed")
  #  getCrosslinking.getCrossLinkingMain(outputPrefix+".filter.cluster.bed",outputPrefix+".reliable_insertion.bed", outputPrefix+"crosslinking.insertion.bed")
  #  getCrosslinking.getCrossLinkingMain(outputPrefix+".filter.cluster.bed",outputPrefix+".reliable_substition.bed", outputPrefix+"crosslinking.substition.bed")
  #  if species is not None:
  #    annotatePeaks.annotatePeak( outputPrefix+".crosslinking.deletion.bed", species, outputPrefix+".crosslinking.deletion.anno.txt")
  #    annotatePeaks.annotatePeak( outputPrefix+".crosslinking.insertion.bed", species,outputPrefix+".crosslinking.insertion.anno.txt")
  #    annotatePeaks.annotatePeak( outputPrefix+".crosslinking.substitution.bed", species, outputPrefix+".crosslinking.substitution.anno.txt")
  #  else:
  #    getCrosslinking.getCrossLinkingMain(outputPrefix+".filter.cluster.bed",outputPrefix+".filter.reliable.bed",outputPrefix+".crosslinking.bed")

  #if species is not None:
  #  annotatePeaks.annotatePeak( outputPrefix+".crosslinking.bed", species,outputPrefix+".crosslinking.anno.txt")

if __name__=="__main__":
  arg_parser = prepare_argparser()
  args = arg_parser.parse_args()

  infile = args.infile                  # Input SAM/BAM file
  outputPrefix = args.outfile           # Output prefix
  matchLength = args.matchLength        # Shorted matched segment length
  mismatch = args.mismatch              # Maximum mismatch number
  pcr = args.clipType                    # PCR removal: (0)no removal; (1)same-start removal; (2)same-seq removal  
  fdrEnrichedCluster = args.fdrCluster  # FDR for enriched clusters
  clipType =args.clipType               # CLIP type (0)HITS-CLIP; (1)PAR-4SU; (2)PAR-6SG; (3)iCLIP
  fdrReliableMutation = args.fdrMutation# FDR for reliable mutations
  species = args.species                # Species ["mm10","hg19"]

  runPipeClip(infile,outputPrefix,matchLength,mismatch,pcr,fdrEnrichedCluster,clipType,fdrReliableMutation,species)
