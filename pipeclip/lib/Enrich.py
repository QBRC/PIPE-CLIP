#!/usr/bin/python
# Programmer : beibei.chen@utsouthwestern.edu
# Usage: Get reliable mutations using binomial distribution
# Input: Filtered BAM, reads coverage (generated by SAMFilter.py), mutation file
# Output: BED
# Last modified: 18 Sep.2013


import datetime
import gc
import logging
import math
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pysam
import rpy2.robjects as robject
from pysam import *
from rpy2.robjects import FloatVector
from rpy2.robjects.packages import importr

from . import Alignment, Utils

gc.enable()

stats = importr("stats")


def freqRank(readCount, rev=False):
    key = sorted(list(readCount.keys()), reverse=rev)
    r_rank = {}
    rank = 0
    for i in key:
        rank += readCount[i]
        r_rank[i] = rank
    return r_rank


def BH(pvalue, pRank, N):
    a = N / float(pRank)
    q = a * pvalue
    qva = max(pvalue, q)
    return qva


# @profile
def KMvalue(posmapfile, negmapfile, mufile):
    """Calculate K(coverage) value for each mutation location Mutations are
    already unique."""
    km = []  # store mutations with updated k value
    km_pair = {}  # Dic of count tuples of (k,m),key:"K_M"
    count = 0
    logging.debug("make wig %s" % str(datetime.datetime.now()))
    poswig = Utils.makeWig(posmapfile)
    negwig = Utils.makeWig(negmapfile)
    start_time = datetime.datetime.now()
    logging.debug("finish making wig %s" % str(start_time))
    for item in mufile:
        count += 1
        if count % 5000 == 0:
            stop_time = datetime.datetime.now()
            logging.debug(
                "Counting K-M for %d mutation sites, using %s"
                % (count, str(stop_time - start_time))
            )
            start_time = stop_time
        st = []
        strand = item.strand
        M = item.score
        K = 0
        # logging.debug("Time begin to pileup is %s" % (str(datetime.datetime.now())))
        if strand == "+":
            try:
                K = poswig[item.chr][str(item.start)]
            except:
                continue

        elif strand == "-":
            try:
                K = negwig[item.chr][str(item.start)]
            except:
                continue
        if K >= M:
            item.updateK(K)
            # print item
            # logging.debug("K value for item %s is %d" % (item, K))
            pair_name = str(K) + "_" + str(M)
            if pair_name in km_pair:
                km_pair[pair_name] += 1
            else:
                km_pair[pair_name] = 1
            # km.append(item)
    gc.collect()
    return km_pair


def KMvalue_test(clip, mutations, chr, chrlen):
    """Calculate K(coverage) value for each mutation location Mutations are
    already unique."""
    km = []  # store mutations with updated k value
    count = 0
    # logging.debug("make wig %s" % str(datetime.datetime.now()))
    posBAM = pysam.Samfile(clip.posfilteredBAM, "rb")
    negBAM = pysam.Samfile(clip.negfilteredBAM, "rb")
    start_time = datetime.datetime.now()
    poswig = Utils.makeWigByChr(posBAM, chr)
    negwig = Utils.makeWigByChr(negBAM, chr)
    stop_time = datetime.datetime.now()
    # logging.debug("Finished making wig for %s using %s" % (chr,str(stop_time-start_time)))
    start_time = stop_time
    for item in mutations:
        count += 1
        if count % 100000 == 0:
            stop_time = datetime.datetime.now()
            logging.debug(
                "Counting K-M for %d mutation sites, using %s"
                % (count, str(stop_time - start_time))
            )
            start_time = stop_time
        st = []
        M = item.score
        K = 0
        strand = item.strand
        if strand == "+":
            try:
                K = poswig[item.start]
            except:
                # log.warning("Did not find mutation in poswig")
                continue
        elif strand == "-":
            try:
                K = negwig[item.start]
            except:
                continue
        if K >= M:
            item.updateK(K)
            pair_name = str(K) + "_" + str(M)
            if pair_name in clip.kmpair:
                clip.kmpair[pair_name] += 1
            else:
                clip.kmpair[pair_name] = 1
    posBAM.close()
    negBAM.close()


def uniq(b):  # b is a list
    uniqElements = []
    for i in b:
        if uniqElements.count(i) == 0:
            uniqElements.append(i)
    uniqElements.sort()
    return uniqElements


def mutationEnrich(clip, threshold=0.01):
    coverage = clip.coverage * 1.0
    totalMuCount = clip.mutationCount
    mutations = []
    total_test = 0
    for chr, chrlen in clip.refInfo:
        # logging.debug(chr)
        try:
            mufile = open(clip.outprefix + "." + chr + ".mutations.bed")
        except:
            logging.info(
                "Cannot open mutation file %s , move on."
                % (clip.outprefix + "." + chr + ".mutations.bed")
            )
            continue
        for record in mufile:
            total_test += 1
            info = record.rstrip().split("\t")
            new_mu = Alignment.MutationBed(
                info[0],
                int(info[1]),
                int(info[2]),
                info[3],
                int(info[4]),
                info[5],
                info[6],
            )
            mutations.append(new_mu)
            try:
                pass
                # os.remove(clip.outprefix+"."+chr+".mutations.bed")
            except:
                pass
        logging.debug(len(mutations))
        KMvalue_test(
            clip, mutations, chr, chrlen
        )  # check after doing KM, if clip.mutations changed
    try:
        os.remove(clip.posfilteredBAM)
        os.remove(clip.negfilteredBAM)
        os.remove(clip.posfilteredBAM + ".bai")
        os.remove(clip.negfilteredBAM + ".bai")
    except:
        pass
    del clip.posfilteredBAM
    del clip.negfilteredBAM
    gc.collect()  # logging.info("Finished K-M counting, starting fitting.")

    R = robject.r
    reliableList = []
    P = totalMuCount / coverage
    km_p = {}  # store km and corresponding p value
    pvalues = []
    for k in clip.kmpair:  # KM_test:
        parameters = k.split("_")
        p = R.pbinom(int(parameters[1]) - 1, int(parameters[0]), P, False)[0]
        pvalues.append(p)
        km_p[k] = p
    pCount = dict(Counter(pvalues))
    pRank = freqRank(pCount, True)
    pqDic = {}
    for i in list(pRank.keys()):
        try:
            p_rank = pRank[i]
            q_value = BH(i, p_rank, total_test)
            pqDic[i] = q_value
        except:
            print("Cannot find p value in dictionary", file=sys.stderr)
            continue
    count = 0
    for mu in mutations:
        name = str(mu.kvalue) + "_" + str(mu.score)
        try:
            mu.pvalue = km_p[name]
        except:
            # logging.debug("Problem with %s" % mu)
            continue
        mu.qvalue = pqDic[mu.pvalue]
        if mu.qvalue <= threshold:
            count += 1
            new_mutationName = "Mutation_" + str(count)
            mu.name = new_mutationName
            mu.sig = True
            clip.sigMutationCount += 1
            clip.addSigToDic(clip.sigMutations, mu)

    mutations = None
    gc.collect()


def mutationEnrichWCtrl(clip, ctrlclip, threshold=0.01):
    coverage = ctrlclip.coverage * 1.0
    totalMuCount = ctrlclip.mutationCount
    mutations = []
    total_test = 0
    for chr, chrlen in clip.refInfo:
        try:
            mufile = open(clip.outprefix + "." + chr + ".mutations.bed")
        except:
            logging.info(
                "Cannot open mutation file %s , move on."
                % (clip.outprefix + "." + chr + ".mutations.bed")
            )
            continue
        for record in mufile:
            total_test += 1
            info = record.rstrip().split("\t")
            new_mu = Alignment.MutationBed(
                info[0],
                int(info[1]),
                int(info[2]),
                info[3],
                int(info[4]),
                info[5],
                info[6],
            )
            mutations.append(new_mu)
            try:
                os.remove(clip.outprefix + "." + chr + ".mutations.bed")
                os.remove(controlclip.outprefix + "." + chr + ".mutations.bed")
            except:
                pass
        KM_test = KMvalue_test(
            clip, mutations, chr, chrlen
        )  # check after doing KM, if clip.mutations changed
    try:
        os.remove(clip.posfilteredBAM)
        os.remove(clip.negfilteredBAM)
        os.remove(clip.posfilteredBAM + ".bai")
        os.remove(clip.negfilteredBAM + ".bai")
    except:
        pass
    del clip.posfilteredBAM
    del clip.negfilteredBAM
    gc.collect()  # logging.info("Finished K-M counting, starting fitting.")

    R = robject.r
    reliableList = []
    P = totalMuCount / coverage
    km_p = {}  # store km and corresponding p value
    pvalues = []
    for k in clip.kmpair:  # KM_test:
        parameters = k.split("_")
        p = R.pbinom(int(parameters[1]) - 1, int(parameters[0]), P, False)[0]
        pvalues.append(p)
        km_p[k] = p
    pCount = dict(Counter(pvalues))
    pRank = freqRank(pCount, True)
    pqDic = {}
    for i in list(pRank.keys()):
        try:
            p_rank = pRank[i]
            q_value = BH(i, p_rank, total_test)
            pqDic[i] = q_value
        except:
            print("Cannot find p value in dictionary", file=sys.stderr)
            continue
    count = 0
    for mu in mutations:
        name = str(mu.kvalue) + "_" + str(mu.score)
        try:
            mu.pvalue = km_p[name]
        except:
            # logging.debug("Problem with %s" % mu)
            continue
        mu.qvalue = pqDic[mu.pvalue]
        if mu.qvalue <= threshold:
            count += 1
            new_mutationName = "Mutation_" + str(count)
            mu.name = new_mutationName
            mu.sig = True
            clip.sigMutationCount += 1
            clip.addSigToDic(clip.sigMutations, mu)

    mutations = None
    gc.collect()


def clusterEnrich(clip, threshold=0.01):
    cluster_filename = (
        clip.outprefix + ".clusters.bed"
    )  # clip.filepath.split("/")[-1].split(".")[0]
    # Call R code and get result
    epsilon = [0.01, 0.15, 0.1]
    step = [0.1, 0.08, 0.05]
    for index in range(len(epsilon)):
        e = epsilon[index]
        s = step[index]
        r_script = Path(__file__).parent.resolve() / "ZTNB_tryCatch.R"
        r_args = [
            "Rscript",
            r_script,
            cluster_filename,
            str(threshold),
            str(e),
            str(s),
        ]
        gc.collect()
        p = subprocess.Popen(r_args)
        stdout_value = p.communicate()[0]
        try:
            r_output_log = open(cluster_filename + ".pipeclip.ztnblog", "r")
            # logging.debug("Log file opened")
            flag = r_output_log.read(1)
            if flag == "Y":  # converged
                break
            elif flag == "N":
                continue
        except:
            logging.info(
                "No log file was produced by R code, continue regression using other parameters anyway."
            )
            continue

    # check ztnb file
    try:
        enrich_parameter = open(cluster_filename + ".pipeclip.ztnb", "r")
    except IOError as message:
        logging.error("Cannot open ztnb result file")
        return False
    nbDic = {}
    for item in enrich_parameter:
        buf = item.rstrip().split("\t")
        if buf[0] != "#":
            nb_key = "_".join(buf[0:2])  # reads_length as key
            # logging.debug("NB_key %s" % nb_key)
            if nb_key not in nbDic:
                nbDic[nb_key] = (buf[2], buf[3])  # pvalue and qvalue
    # logging.info("There are %d read-length pairs" % (len(nbDic.keys())))
    if len(list(nbDic.keys())) == 0:
        logging.error("There are no read-length pairs found by ZTNB. Exit.")
        return False
    else:
        for i in range(len(clip.clusters)):
            r_key = (
                str(clip.clusters[i].score)
                + "_"
                + str(clip.clusters[i].stop - clip.clusters[i].start)
            )
            # logging.debug("Keys from clip.clusters,%s" % r_key)
            if r_key in nbDic:
                clip.clusters[i].pvalue = nbDic[r_key][0]
                clip.clusters[i].qvalue = nbDic[r_key][1]
                clip.clusters[i].sig = True
                clip.sigClusterCount += 1
        nbDic = None
        try:
            os.remove(cluster_filename)
        except:
            pass
        if clip.sigClusterCount == 0:
            return False
        else:
            return True


def clusterEnrich_outsource(clip, threshold=0.01):
    cluster_filename = (
        clip.outprefix + ".clusters.bed"
    )  # clip.filepath.split("/")[-1].split(".")[0]
    # Call R code and get result
    # epsilon = [0.01,0.15,0.1]
    # step = [0.1,0.08,0.05]

    sh_script = Path(__file__).parent.resolve() / "runR1.sh"
    sh_args = ["sh", sh_script, cluster_filename, str(threshold)]
    p = subprocess.Popen(sh_args)
    stdout_value = p.communicate()[0]

    # check ztnb file
    try:
        enrich_parameter = open(cluster_filename + ".pipeclip.ztnb", "r")
    except IOError as message:
        logging.error("Cannot open ztnb result file")
        return False
    nbDic = {}
    for item in enrich_parameter:
        buf = item.rstrip().split("\t")
        if buf[0] != "#":
            nb_key = "_".join(buf[0:2])  # reads_length as key
            # logging.debug("NB_key %s" % nb_key)
            if nb_key not in nbDic:
                nbDic[nb_key] = (buf[2], buf[3])  # pvalue and qvalue
    # logging.info("There are %d read-length pairs" % (len(nbDic.keys())))
    if len(list(nbDic.keys())) == 0:
        logging.error("There are no read-length pairs found by ZTNB. Exit.")
        return False
    else:
        for i in range(len(clip.clusters)):
            r_key = (
                str(clip.clusters[i].score)
                + "_"
                + str(clip.clusters[i].stop - clip.clusters[i].start)
            )
            # logging.debug("Keys from clip.clusters,%s" % r_key)
            if r_key in nbDic:
                clip.clusters[i].pvalue = nbDic[r_key][0]
                clip.clusters[i].qvalue = nbDic[r_key][1]
                clip.clusters[i].sig = True
                clip.sigClusterCount += 1
        nbDic = None
        try:
            os.remove(cluster_filename)
        except:
            pass
        if clip.sigClusterCount == 0:
            return False
        else:
            return True


def fisherTest(clusterp, mutationp):
    R = robject.r
    min_mp = min(mutationp)
    product = clusterp * min_mp
    if product == 0:
        fps = 0
    else:
        xsq = -2 * math.log(clusterp * min_mp)
        fp = R.pchisq(xsq, **{"df": 4, "lower.tail": False, "log.p": True})[0]
        fps = math.exp(fp)
    return fps
