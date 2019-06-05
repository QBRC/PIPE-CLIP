PIPE-CLIP
=========
This is version 1.1.0
Pipeline for CLIP-seq Analysis.
- Galaxy site: the galaxy service is discontinued in 2019.
- Publication: http://genomebiology.com/2014/15/1/R18
- Google Code site: https://code.google.com/p/pipe-clip/


Requirement:
-  Python 2.7; 
-  R 3.0 and above;
-  Perl 5 and above;
-  Python packages: pysam, pybedtools and ghmm;
-  R packages: MASS,VGAM and their dependencies.
-  Other packages: HOMER and annotation files


Installation tips:
- Make sure HOMER are in your PATH. You can test this by type "annotatePeaks.pl" from anywhere and you should get help information of this command.
 
How to use:

- After unzip the package, you cd into the program folder and run PIPE-CLIP by typing:
- python pipeclip.py -i input.bam -o output_prefix -c CLIP_type -l minimum_matchlength -m maximum_mismatchcount  -r Remove_PCR_duplicate -M FDR_for_mutations -C FDR_for_clusters -s species

- -i input BAM
- -o output prefix
- -c CLIP type,[0,1,2,3] (0)HITS-CLIP; (1)PAR-4SU; (2)PAR-6SG; (3)iCLIP
- -l minimum match length
- -m maximum mismatch count
- -r method to remove PCR duplicate,[0,1,2] (0)No removal; (1)Remove by read start; (2)Remove by sequence
- -M FDR to get significant mutations
- -C FDR to get enriched clusters
- -s species 
Here, species might be hg19.

Contact: Zhiqun.Xie@UTSouthwestern.edu
