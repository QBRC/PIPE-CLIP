#Main pipeline connects all the scripts together
#Original Programmer: beibei.chen@utsouthwestern.edu
#Refactored by: eric.roos@utsouthwestern.edu
#Usage: python pipeclip.py input.sam output_prefix match_length mismatch_number pcr_rm fdr_cluster clip_type fdr_mutation species
#Required packages: pysam, ghmm, pybedtools
#Last modification: 3 March 2014

import sys
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
  argparser = ap.ArgumentParser(description=description, epilog = epilog)
  argparser.add_argument("-i","--input",dest = "infile", type = str, required = True, help = "input bam file")
  argparser.add_argument("-o","--output",dest = "outfile", type = str,required = True, help = "output file, default is stdout")
  argparser.add_argument("-l","--matchLength",dest = "matchLength", type = int ,required = True, help = "shorted matched segment length")
  argparser.add_argument("-m","--mismatch",dest = "mismatch", type = str,required = True, help = "maximum mismatch number")
  argparser.add_argument("-c","--clipType",dest = "clipType", type = str,required = True, help = "CLIP type (0)HITS-CLIP; (1)PAR-4SU; (2)PAR-6SG; (3)iCLIP")
  argparser.add_argument("-M","--fdrMutation",dest = "fdr", type = str,required = True, help = "FDR for reliable mutations")
  argparser.add_argument("-C","--fdrCluster",dest = "fdr", type = str,required = True, help = "FDR for enriched clusters")
  argparser.add_argument("-s","--species",dest = "species", type = str,required = True, help = "Species [\"mm10\",\"hg19\"]")
  return(argparser)

def runPipeClip(infile,outputPrefix,matchLength,mismatch,pcr,fdrEnrichedCluster,clipType,fdrReliableMutation,species):
  ########################## Process input #######################
  inputProcess.inputProcessMain(infile,outputPrefix)
  SAMFilter.SAMFILTERMain(outputPrefix+".sorted.bam",outputPrefix+".filter",matchLength,mismatch,pcr,clipType)

  ######################### Enrich clusters ######################
  mergeReadsMain(outputPrefix+".filter.bam",outputPrefix+".filter.merge")

  ######################### Enrich clusters ######################
  call(["./Rscript","ZTNB.R",outputPrefix+".filter.merge",fdrEnrichedCluster])
  getClusterMain(outputPrefix+".filter.merge",outputPrefix+".merge.ztnb", outputPrefix+".filter.cluster.bed")

  ########################### Mutation #############################
  if clipType == "3":
    findTruncationMain(outputPrefix+".filter.bam")
  else:
    findMutationMain(outputPrefix+".filter.bam",outputPrefix+".filter.merge.ztnb",clipType);
  mutationFilterMain(outputPrefix+".filter.bam",outputPrefix+".filter.mutation.bed",outputPrefix+".filter.reliable",clipType,fdrReliableMutation,outputPrefix+".filter.coverage")

  ######################### Merge and annotation ################
  if clipType == "0":
    getCrossLinkingMain(outputPrefix+".filter.cluster.bed",outputPrefix+".reliable_deletion.bed", outputPrefix+"crosslinking.deletion.bed")
    getCrossLinkingMain(outputPrefix+".filter.cluster.bed",outputPrefix+".reliable_insertion.bed", outputPrefix+"crosslinking.insertion.bed")
    getCrossLinkingMain(outputPrefix+".filter.cluster.bed",outputPrefix+".reliable_substition.bed", outputPrefix+"crosslinking.substition.bed")
    if species is not None:
      print "not yet implemented"
      #annotatePeaks.pl $2.crosslinking.deletion.bed $9 > $2.crosslinking.deletion.anno.txt
      #annotatePeaks.pl $2.crosslinking.insertion.bed $9 > $2.crosslinking.insertion.anno.txt
      #annotatePeaks.pl $2.crosslinking.substitution.bed $9  > $2.crosslinking.substitution.anno.txt
    else:
      getCrossLinkingMain(outputPrefix+".filter.cluster.bed",outputPrefix+".filter.reliable.bed",outputPrefix+".crosslinking.bed")

  if species is not None:
    print "Not yet implemented"
    #annotatePeaks.pl $2.crosslinking.bed $9  > $2.crosslinking.anno.txt
if __name__=="__main__":
  arg_parser = prepare_argparser()
  args = argparser.parse_args()

  infile = args.infile                  # Input SAM/BAM file
  outputPrefix = args.outfile           # Output prefix
  matchLength = args.matchLength        # Shorted matched segment length
  mismatch = args.mismatch              # Maximum mismatch number
  pcr = arg.clipType                    # PCR removal: (0)no removal; (1)same-start removal; (2)same-seq removal  
  fdrEnrichedCluster = args.fdrCluster  # FDR for enriched clusters
  clipType =args.clipType               # CLIP type (0)HITS-CLIP; (1)PAR-4SU; (2)PAR-6SG; (3)iCLIP
  fdrReliableMutation = args.fdrMutation# FDR for reliable mutations
  species = args.species                # Species ["mm10","hg19"]

  runPipeClip(infile,outputPrefix,matchLength,mismatch,pcr,fdrEnrichedCluster,clipType,fdrReliableMutation,species)
