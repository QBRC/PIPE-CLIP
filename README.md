# PIPE-CLIP

Pipeline for CLIP-seq Analysis.
This is a fork version (2.0.0) of [PIPE-CLIP](https://github.com/QBRC/PIPE-CLIP).

## Requirement:

- Python >=3.6;
- R >=3.0
- Perl >=5.0
- Python packages: `pysam`, `pybedtools` and `ghmm`;
- R packages: `MASS`, `VGAM` and their dependencies.
- Other packages: `HOMER` (annotatePeaks.pl) and annotation files
  - Make sure HOMER are in your PATH. You can test this by type "annotatePeaks.pl" from anywhere and you should get help information of this command.

## How to use:

- After unzip the package, you cd into the program folder and run PIPE-CLIP by typing:

```bash
python pipeclip.py -i input.bam -o output_prefix -c CLIP_type -l minimum_matchlength -m maximum_mismatchcount  -r Remove_PCR_duplicate -M FDR_for_mutations -C FDR_for_clusters -s species
```

- `-i` input BAM
- `-o` output prefix
- `-c` CLIP type,[0,1,2,3] (0)HITS-CLIP; (1)PAR-4SU; (2)PAR-6SG; (3)iCLIP
- `-l` minimum match length
- `-m` maximum mismatch count
- `-r` method to remove PCR duplicate,[0,1,2] (0)No removal; (1)Remove by read start; (2)Remove by sequence
- `-M` FDR to get significant mutations
- `-C` FDR to get enriched clusters
- `-s` species. (species might be hg19, mm10, mm9.)

## Footnote

- Contact: Zhiqun.Xie@UTSouthwestern.edu
- Publication: http://genomebiology.com/2014/15/1/R18
- Google Code site: https://code.google.com/p/pipe-clip/
