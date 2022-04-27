#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2022 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2022-04-26 18:14

"""Fork from:

-  Main pipeline connects all the scripts together
-  Original Programmer: beibei.chen@utsouthwestern.edu
-  Refactored by: eric.roos@utsouthwestern.edu
-  Usage: python pipeclip.py input.sam output_prefix match_length mismatch_number pcr_rm fdr_cluster clip_type fdr_mutation species
-  Required packages: pysam, pybedtools, rpy2
"""


import pathlib
import sys

import click

from .lib import CLIP, Enrich, Utils

LOGGER = Utils.get_logger(__name__)


@click.command(help="PIPECLIP v2.0.7", no_args_is_help=True)
@click.option(
    "--infile", "-i", "infile", help="Input bam file.", required=True
)
@click.option(
    "--control", "-t", "control", help="Control bam file.", required=False
)
@click.option(
    "--output",
    "-o",
    "output_prefix",
    help="Output file, default is stdout.",
    required=True,
)
@click.option(
    "--matchLength",
    "-l",
    "matchLength",
    help="Minimum matched segment length.",
    type=int,
    required=True,
)
@click.option(
    "--mismatch",
    "-m",
    "mismatch",
    help="Maximum mismatch number.",
    type=int,
    required=True,
)
@click.option(
    "--clipType",
    "-c",
    "clipType",
    help="CLIP type (0)HITS-CLIP; (1)PAR-4SU; (2)PAR-6SG; (3)iCLIP",
    type=click.IntRange(0, 3),
    required=True,
)
@click.option(
    "--rmdup",
    "-r",
    "rmdup",
    help="Remove PCR duplicate (0)No removal; (1)Remove by read start; (2)Remove by sequence;",
    type=click.IntRange(0, 2),
    required=True,
)
@click.option(
    "--fdrMutation",
    "-M",
    "fdrReliableMutation",
    help="FDR for reliable mutations.",
    type=float,
    required=True,
)
@click.option(
    "--fdrCluster",
    "-C",
    "fdrEnrichedCluster",
    help="FDR for enriched clusters.",
    type=float,
    required=True,
)
@click.option(
    "--species",
    "-s",
    "species",
    help="Species (hg19, mm9, mm10...)",
    type=click.Choice(["mm9", "mm10", "hg19"]),
    required=False,
)
def runPipeClip(
    infile,
    control,
    output_prefix,
    matchLength,
    mismatch,
    rmdup,
    fdrEnrichedCluster,
    clipType,
    fdrReliableMutation,
    species,
):
    # make output directory
    pathlib.Path(output_prefix).parent.mkdir(parents=True, exist_ok=True)
    # run clip analysis
    myClip = CLIP.CLIP(infile, output_prefix)
    controlFlag = False
    if control != None:
        controlClip = CLIP.CLIP(control, output_prefix + "Control")
    LOGGER.info("Start to run")
    if myClip.testInput():  # check input
        LOGGER.info("Input file OK,start to run PIPE-CLIP")
        LOGGER.info("Species info %s" % species)
        if control != None:  # test control file
            if controlClip.testInput():
                LOGGER.info(
                    "Control file OK. Use control in mutation enrichment."
                )
                controlFlag = True
            else:
                LOGGER.info(
                    "Control file format error. Continue without control."
                )
        if myClip.readfile():
            myClip.filter2(matchLength, mismatch, clipType, rmdup)
            if controlFlag:
                LOGGER.info("Read in control file")
                controlClip.readfile()
                controlClip.filter2(matchLength, mismatch, clipType, rmdup)
            # myClip.printClusters()
            # myClip.printMutations()
            if myClip.clusterCount > 0:
                LOGGER.info("Get enriched clusters")
                status = Enrich.clusterEnrich_outsource(
                    myClip, fdrEnrichedCluster
                )
                if status:
                    LOGGER.info(
                        "Found %d enriched clusters" % myClip.sigClusterCount
                    )
                    myClip.printEnrichedClusters()
                else:
                    LOGGER.error(
                        "There is no enriched cluster found. Exit program"
                    )
                    sys.exit(1)
            else:
                LOGGER.error(
                    "There is no clusters found. Please check input.Exit program."
                )
                sys.exit(1)

            if myClip.mutationCount > 0:
                LOGGER.info("Get reliable mutations")
                if controlFlag:  # use control
                    Enrich.mutationEnrichWCtrl(
                        myClip, controlClip, fdrReliableMutation
                    )
                else:
                    Enrich.mutationEnrich(myClip, fdrReliableMutation)
                LOGGER.info(
                    "There are %d reliable mutations" % myClip.sigMutationCount
                )
                myClip.printReliableMutations()
            else:
                LOGGER.warning("There is no mutation found in this BAM file.")
            # Start to get crosslinking sites
            if myClip.sigClusterCount > 0 and myClip.sigMutationCount > 0:
                LOGGER.info("Get cross-linking sites")
                myClip.getCrosslinking()
                if len(list(myClip.crosslinking.keys())) > 0:
                    outfilelist = myClip.printCrosslinkingSites()
                    myClip.printCrosslinkingMutations()
                else:
                    LOGGER.warning(
                        "There is no crosslinking found. May be caused by no reliable mutations in enriched clusters. Print out enriched clusters instead."
                    )
                    outfilelist = myClip.printEnrichedClusters()
            else:
                if myClip.sigClusterCount <= 0:
                    LOGGER.error(
                        "There is no enriched clusters for this sample, please check your input file. Exit."
                    )
                    sys.exit(2)
                elif myClip.sigMutationCount <= 0:
                    LOGGER.warning(
                        "There is no reliable mutations found. PIPE-CLIP will provide enriched clusters as crosslinking candidates."
                    )
                    outfilelist = myClip.printEnrichedClusters()
            # annotation if possible
        if species in ["mm10", "mm9", "hg19"]:
            LOGGER.info("Started to annotate cross-linking sits using HOMER")
            for name in outfilelist:
                # LOGGER.debug("Start to do annotation for %s" % name)
                Utils.annotation(name, species)
        # output a status log file
        logfile = open(output_prefix + ".pipeclip.summary.log", "w")
        print("PIPE-CLIP run finished. Parameters are:", file=logfile)
        print(
            "Input BAM: %s \nOutput prefix: %s \nMinimum matched length: %d \nMaximum mismatch count: %d \nPCR duplicate removal code: %d \nFDR for enriched clusters: %f \nFDR for reliable mutations: %f"
            % (
                infile,
                output_prefix,
                matchLength,
                mismatch,
                rmdup,
                fdrEnrichedCluster,
                fdrReliableMutation,
            ),
            file=logfile,
        )
        print(
            "There are %d mapped reads in input BAM file. After filtering,%d reads left"
            % (myClip.originalMapped, myClip.filteredAlignment),
            file=logfile,
        )
        print(
            "%d out of %d clusters are enriched."
            % (myClip.sigClusterCount, len(myClip.clusters)),
            file=logfile,
        )
        print(
            "%d out of %d mutations are reliable."
            % (myClip.sigMutationCount, myClip.mutationCount),
            file=logfile,
        )
        print(
            "%d crosslinking site candidates are found, with %d supporting reliable mutations."
            % (
                len(list(myClip.crosslinking.keys())),
                len(myClip.crosslinkingMutations),
            ),
            file=logfile,
        )
        logfile.close()
        LOGGER.info(
            "PIPE-CLIP finished the job, please check your results. :)"
        )
    else:
        LOGGER.error("File corruputed, program exit.")
        sys.exit(0)


if __name__ == "__main__":
    runPipeClip()
