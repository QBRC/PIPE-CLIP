#Main pipeline connects all the scripts together
#Original Programmer: beibei.chen@utsouthwestern.edu
#Refactored by: eric.roos@utsouthwestern.edu
#Usage: python pipeclip.py input.sam output_prefix match_length mismatch_number pcr_rm fdr_cluster clip_type fdr_mutation species
#Required packages: pysam, ghmm, pybedtools
#Last modification: 27 Feb 2014

import sys
import annotatePeaks
import barcodeRemoval
import findMutation
import findTruncation
import getCluster
import getCrosslinking
import inputProcess
import mergeReads
import mutationFilter
from subprocess import call

if __name__=="__main__":
  infile = sys.argv[1]              # Input SAM/BAM file
  outputPrefix = sys.argv[2]        # Output prefix
  matchLength = sys.argv[3]         # Shorted matched segment length
  mismatch = sys.argv[4]            # Maximum mismatch number
  pcr = sys.argv[5]                 # PCR removal: (0)no removal; (1)same-start removal; (2)same-seq removal  
  fdrEnrichedCluster = sys.argv[6]  # FDR for enriched clusters
  clipType =sys.argv[7]             # CLIP type (0)HITS-CLIP; (1)PAR-4SU; (2)PAR-6SG; (3)iCLIP
  fdrReliableMutation = sys.argv[8] # FDR for reliable mutations
  species = sys.argv[9]             # Species ["mm10","hg19"]

  ########################## Process input #######################
  inputProcessMain(infile,outputPrefix)
  SAMFILTERMain(outputPrefix+".sorted.bam",outputPrefix+".filter",matchLength,mismatch,pcr,clipType)

  ######################### Enrich clusters ######################
  mergeReadsMain(outputPrefix+".filter.rehead.bam")

  ######################### Enrich clusters ######################
  call(["./Rscript","ZTNB.R",outputPrefix+".filter.rehead.merge",fdrEnrichedCluster])
  getClusterMain(outputPrefix+".filter.rehead.merge",outputPrefix+".reheadmerge.ztnb", outputPrefix+".filter.cluster.bed")

  ########################### Mutation #############################
  if clipType == "3":
    findTruncationMain(outputPrefix+".filter.rehead.bam")
  else
    findMutationMain(outputPrefix+".filter.rehead.bam",outputPrefix+".filter.rehead.merge.ztnb",clipType);
  mutationFilterMain(outputPrefix+".filter.rehead.bam",outputPrefix+".filter.mutation.bed",outputPrefix+".filter.reliable",clipType,fdrReliableMutation,outputPrefix+".filter.coverage")

  ######################### Merge and annotation ################
  if clipType == "0":
    getCrossLinkingMain(outputPrefix+".filter.cluster.bed",outputPrefix+".reliable_deletion.bed", outputPrefix+"crosslinking.deletion.bed")
    getCrossLinkingMain(outputPrefix+".filter.cluster.bed",outputPrefix+".reliable_insertion.bed", outputPrefix+"crosslinking.insertion.bed")
    getCrossLinkingMain(outputPrefix+".filter.cluster.bed",outputPrefix+".reliable_substition.bed", outputPrefix+"crosslinking.substition.bed")
    if species is not None
        #annotatePeaks.pl $2.crosslinking.deletion.bed $9 > $2.crosslinking.deletion.anno.txt
        #annotatePeaks.pl $2.crosslinking.insertion.bed $9 > $2.crosslinking.insertion.anno.txt
        #annotatePeaks.pl $2.crosslinking.substitution.bed $9  > $2.crosslinking.substitution.anno.txt
    else
      getCrossLinkingMain(outputPrefix+".filter.cluster.bed",outputPrefix+".filter.reliable.bed",outputPrefix+".crosslinking.bed")
      if species is not None
        #annotatePeaks.pl $2.crosslinking.bed $9  > $2.crosslinking.anno.txt
