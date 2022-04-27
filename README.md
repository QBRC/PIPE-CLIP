# PIPE-CLIP

Pipeline for CLIP-seq Analysis.
This is a fork version (2.0.0) of [PIPE-CLIP](https://github.com/QBRC/PIPE-CLIP).

## Requirement:

- Python >=3.6;
- Python packages: `pysam`, `pybedtools` and `rpy2`. (Python packages will be installed automaticallly)

- R >=3.0;
- R packages: `MASS`, `VGAM` and their dependencies. (R packages will be installed automatically)

- Perl >=5.0
- Other packages: `HOMER` (annotatePeaks.pl) and annotation files
  - Make sure HOMER are in your PATH. You can test this by type "annotatePeaks.pl" from anywhere and you should get help information of this command.

## How to install:

```bash
pip install pipeclip
```

## How to use:

```bash
pipeclip -i input.bam -o output_prefix -c CLIP_type -l minimum_matchlength -m maximum_mismatchcount -r Remove_PCR_duplicate -M FDR_for_mutations -C FDR_for_clusters -s species
```

- `-i` input BAM
- `-t` control BAM
- `-o` output prefix
- `-c` CLIP type,[0,1,2,3] (0)HITS-CLIP; (1)PAR-4SU; (2)PAR-6SG; (3)iCLIP
- `-r` method to remove PCR duplicate,[0,1,2] (0)No removal; (1)Remove by read start; (2)Remove by sequence
- `-l` minimum match length
- `-m` maximum mismatch count
- `-M` FDR to get significant mutations
- `-C` FDR to get enriched clusters
- `-s` species. (species might be hg19, mm10, mm9.) Leave blank to skip annotation step.

## Footnote

> Cite:

- Chen, B., Yun, J., Kim, M.S. et al. PIPE-CLIP: a comprehensive online tool for CLIP-seq data analysis. Genome Biol 15, R18 (2014). https://doi.org/10.1186/gb-2014-15-1-r18
