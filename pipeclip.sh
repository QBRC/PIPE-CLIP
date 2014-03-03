#Main pipeline connects all the scripts together
#Programmer: beibei.chen@utsouthwestern.edu
#Usage: sh pipeclip.sh input.sam output_prefix match_length mismatch_number pcr_rm fdr_cluster clip_type fdr_mutation species
#Required packages: pysam, ghmm, pybedtools, samtools
#Required parameters
#$1: Input SAM/BAM file
#$2: Output prefix
#$3: Shorted matched segment length
#$4: Maximum mismatch number
#$5: PCR removal: (0)no removal; (1)same-start removal; (2)same-seq removal  
#$6: FDR for enriched clusters
#$7: CLIP type (0)HITS-CLIP; (1)PAR-4SU; (2)PAR-6SG; (3)iCLIP
#$8: FDR for reliable mutations
#$9: Species ["mm10","hg19"]
#
#Last modification: 19 Dec 2013
#
##########################Process input#######################
python inputProcess.py $1 $2  #output $2.sorted.bam
samtools index $2.sorted.bam
samtools view -H $2.sorted.bam > $2.header
python SAMFilter.py -i $2.sorted.bam -o $2.filter -m $3 -n $4 -r $5 -t $7
samtools reheader $2.header $2.filter.bam > $2.filter.rehead.bam
samtools index $2.filter.rehead.bam
#
#
#########################Clustering filtered mapped reads######
python mergeReads.py $2.filter.rehead.bam > $2.filter.rehead.merge
#
#
#########################Enrich clusters######################
./Rscript ZTNB.R $2.filter.rehead.merge $6
python getCluster.py $2.filter.rehead.merge $2.filter.rehead.merge.ztnb > $2.filter.cluster.bed


#BC##########################Mutation#############################
if test "$7" = "3"  #iCLIP
	then
		python findTruncation.py $2.filter.rehead.bam > $2.filter.mutation.bed
	else
		python findMutation.py -i $2.filter.rehead.bam -o $2.filter.mutation.bed -p $7
fi
python mutationFilter.py -a $2.filter.rehead.bam -b $2.filter.mutation.bed -o $2.filter.reliable -p $7 -f $8 -c $2.filter.coverage
#
#
#########################Merge and annotation################
if test "$7" = "0" #HITS-CLIP
	then
	python getCrosslinking.py $2.filter.cluster.bed $2.filter.reliable_deletion.bed > $2.crosslinking.deletion.bed
	python getCrosslinking.py $2.filter.cluster.bed $2.filter.reliable_insertion.bed > $2.crosslinking.insertion.bed
	python getCrosslinking.py $2.filter.cluster.bed $2.filter.reliable_substitution.bed > $2.crosslinking.substitution.bed
#	echo "HITS finish"
	if test "$9" != "Null"
		then
			annotatePeaks.pl $2.crosslinking.deletion.bed $9 > $2.crosslinking.deletion.anno.txt
			annotatePeaks.pl $2.crosslinking.insertion.bed $9 > $2.crosslinking.insertion.anno.txt
			annotatePeaks.pl $2.crosslinking.substitution.bed $9  > $2.crosslinking.substitution.anno.txt
	fi
else
	python getCrosslinking.py $2.filter.cluster.bed $2.filter.reliable.bed > $2.crosslinking.bed
	if test "$9" != "Null"
	then
		annotatePeaks.pl $2.crosslinking.bed $9  > $2.crosslinking.anno.txt
	fi
fi
#
