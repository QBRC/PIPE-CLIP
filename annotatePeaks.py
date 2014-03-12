# Author: Min Soo Kim
# Usage: Annotate peaks for hg19 and mm10.
# Input: BED file
# Output: BED file with additional columns for annotation

import sys
import os
import pybedtools

def cleanBed(feature):
    return feature[0:3]

def returnFeatureWithGene(feature, region):
    geneNamePosition = 6
    region = [region]
    if int(feature[9]) == 0:
       return feature[0:3] + feature[geneNamePosition].split('/') + region

def returnRealPeaks(feature,region):
    region = [region]
    if int(feature[9]) != -1:
       return feature[0:3] + feature[6].split('/') + region

def returnUnannotatedPeaks(feature,region):
    region = [region]
    if int(feature[9]) == -1:
       return feature[0:3] + feature[6].split('/') + region

def returnLeftoverPeaks(feature):
    if int(feature[9]) != 0:
       return feature

def promoter(feature):
    if feature.strand == '+':
       feature.end = feature.start
       feature.start = feature.start - 1000
    elif feature.strand == '-':
       feature.start = feature.end
       feature.end = feature.end + 1000
    return feature

def upstream(feature):
    if feature.strand == '+':
       feature.end = feature.start - 1000
       feature.start = feature.start - 10000
    elif feature.strand == '-':
       feature.start = feature.end + 1000
       feature.end = feature.end + 10000
    return feature

def downstream(feature):
    if feature.strand == '+':
       feature.start = feature.end
       feature.end = feature.end + 5000
    elif feature.strand == '-':
       feature.start = feature.end + 1000
       feature.end = feature.end + 10000
    return feature

def annotatePeaks(peakFile,genome):
    pybedtools.BedTool(peakFile).each(cleanBed).moveto('peaks.bed')

    if (genome=='hg19'):
        genome = pybedtools.BedTool('hg19.RefSeq.bed')
    elif (genome=='mm10'):
        genome = pybedtools.BedTool('mm10.RefSeq.bed')
    else:
        sys.exit("Unknown or unavailable genome.")

    genebodyPeaks = pybedtools.BedTool('peaks.bed').closest(genome,d=True,t='first')
    genebodyPeaks.each(returnLeftoverPeaks).each(cleanBed).moveto('leftoverpeaks.bed')
    genebodyPeaks = genebodyPeaks.each(returnFeatureWithGene,'GeneBody')

    promotergenome = genome.each(promoter)
    promoterPeaks = pybedtools.BedTool('leftoverpeaks.bed').closest(promotergenome,d=True,t='first')
    promoterPeaks.each(returnLeftoverPeaks).each(cleanBed).moveto('leftoverpeaks.bed')
    promoterPeaks = promoterPeaks.each(returnFeatureWithGene,'Promoter')

    upstreamgenome = genome.each(upstream)
    upstreamPeaks = pybedtools.BedTool('leftoverpeaks.bed').closest(upstreamgenome,d=True,t='first') 
    upstreamPeaks.each(returnLeftoverPeaks).each(cleanBed).moveto('leftoverpeaks.bed')
    upstreamPeaks = upstreamPeaks.each(returnFeatureWithGene,'Upstream')

    downstreamgenome = genome.each(downstream)
    downstreamPeaks = pybedtools.BedTool('leftoverpeaks.bed').closest(downstreamgenome,d=True,t='first')
    downstreamPeaks.each(returnLeftoverPeaks).each(cleanBed).moveto('leftoverpeaks.bed')
    downstreamPeaks = downstreamPeaks.each(returnFeatureWithGene,'Downstream')

    intergenicPeaks = pybedtools.BedTool('leftoverpeaks.bed').closest(genome,d=True,t='first').each(returnRealPeaks,'Intergenic')
    unannotatedPeaks = pybedtools.BedTool('leftoverpeaks.bed').closest(genome,d=True,t='first').each(returnUnannotatedPeaks,'No Annotation')

#Concat BED. Write to file
    genebodyPeaks.moveto('genebodyPeaks.bed')
    promoterPeaks.moveto('promoterPeaks.bed')
    upstreamPeaks.moveto('upstreamPeaks.bed')
    downstreamPeaks.moveto('downstreamPeaks.bed')
    intergenicPeaks.moveto('intergenicPeaks.bed')
    unannotatedPeaks.moveto('unannotatedPeaks.bed')

    filenames = ['genebodyPeaks.bed','promoterPeaks.bed','upstreamPeaks.bed','downstreamPeaks.bed','intergenicPeaks.bed','unannotatedPeaks.bed']
    with open('Annotated.bed', 'w') as outfile:
       for fname in filenames:
          with open(fname) as infile:
             for line in infile:
                outfile.write(line)
    
    
    filenames = ['genebodyPeaks.bed','promoterPeaks.bed','upstreamPeaks.bed','downstreamPeaks.bed','intergenicPeaks.bed','unannotatedPeaks.bed','leftoverpeaks.bed','peaks.bed']
    for fname in filenames:
        os.remove(fname)
