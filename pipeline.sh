#Require pysam,ghmm package, python2.7 or higher version
#Require bedtools and samtools
#10 required parameters
#$1: Input SAM file
#$2: Output prefix
#$3: Shortest matched segment length
#$4: Maximum mismatch number
#$5: Remove PCR duplicates or not
#$6: FDR to determine enriched cluster
#$7: Minimum reads number in a cluster
#$8: CLIP type (0:HITS-CLIP; 1:PAR-CLIP(4SU); 2:PAR-CLIP(6SG); 3:iCLIP)
#$9: FDR to determine reliable mutations
#$10: species, recently, should be "mm10" or "hg19"


############Process input################################
#Get matched entries and make them into sorted bam and index
#Use python program to pre-process input file
#########################################################
python inputProcess.py $1 $2  #output $2.sorted.bam
samtools index $2.sorted.bam
samtools view -H $2.bam > header
############Process input ends###########################


#############Filter BAM by customized parameters#########
#-i: input file
#-o: output file prefix, prefix.fileter.sam will be provided to user
#-m: Shortest matched segment length
#-n: Maximum number of mismatch
#########################################################
python SAMFilter_rmdupV2.py -i $2.sorted.bam -o $2.filter -m $3 -n $4 -r $5 -t $8
samtools reheader header $2.filter.bam > $2.filter.rehead.bam # should be sorted
samtools index $2.filter.rehead.bam
############Input filtering ends#########################


#############Clustering the mapped reads################
#-i: input file
#-f: FDR threshold
#-n: minimum reads number required in a cluster
########################################################
python clusterEnrich.py -i $2.filter.rehead.bam  -f $6 -n $7  > $2.filter.cluster.bed
#############Enriched cluster ends#######################


#############Looking for mutations #########################
# This segment gets the mutation file and finds out reliable mutations
# Programs will treat iCLIP separately
# Parameters in mutaionFilter_pvalue.py:
#-a: Input BAM file - in order to calculate (K,M)
#-b: Original mutation sites (can be repetitive)
#-o: Output prefix
#-p: CLIP-seq type, '3' for iCLIP
#-f: FDR threshold
#-c: Coverage file, sum of reads length after filtering, generated in filtering step
############################################################

if test "$8" = "3" #iCLIP
	then
		python findTruncation.py $2.filter.rehead.bam > $2.filter.mutation.bed
	else
		python findMutationV2.py -i $2.filter.rehead.bam -o $2.filter.mutation.bed -p $8
fi

python mutationFilter_pvalue.py -a $2.filter.rehead.bam -b $2.filter.mutation.bed -o $2.filter.reliable -p $8 -f $9 -c $2.filter.coverage
##############Mutation filtering ends########################

##############Intersect cluster with mutations ###############
# Final process: intersect clusters and mutations
# Treated HITS-CLIP separately
##############################################################
if test "$8" = "0" #HITS-CLIP, 3 output
	then
	#echo $8,"HITS"
		python finalMerge.py $2.filter.cluster.bed $2.filter.reliable_deletion.bed > $2.clr.deletion.bed
		python finalMerge.py $2.filter.cluster.bed $2.filter.reliable_insertion.bed > $2.clr.insertion.bed
		python finalMerge.py $2.filter.cluster.bed $2.filter.reliable_substitution.bed > $2.clr.substitution.bed
		annotatePeaks.pl $2.clr.deletion.bed ${10} > $2.clr.deletion.annotation
		annotatePeaks.pl $2.clr.insertion.bed ${10} > $2.clr.insertion.annotation
		annotatePeaks.pl $2.clr.substitution.bed ${10} > $2.clr.substitution.annotation
fi	

if test "$8" != "0" #PAR-CLIP and iCLIP, 1 output
	then
	#echo $8,"PAR"
		python finalMerge.py $2.filter.cluster.bed $2.filter.reliable.bed > $2.clr.bed
		annotatePeaks.pl $2.clr.bed ${10} > $2.clr.clr.annotation
fi
###############End of final merge###############################

###########End of pipeline######################################

