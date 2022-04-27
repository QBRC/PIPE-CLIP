#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Define CLIP class."""

import datetime
import gc
import gzip
import random

import pysam

from . import Alignment, Mutation2, Utils

LOGGER = Utils.get_logger(__name__)

gc.enable()


class CLIP:
    def __init__(self, fileaddr, prefix):
        self.filepath = fileaddr
        self.originalBAM = None
        self.originalMapped = 0
        self.refInfo = []
        self.posfilteredBAM = None
        self.negfilteredBAM = None
        self.filteredAlignment = 0
        self.type = 0
        self.outprefix = prefix
        self.currentGroupKey = "None"
        self.currentGroup = []  # Used in rmdup
        self.previousQul = [0, 0, 0]  # for rmdup,[matchlen,mapq,mismatch]
        self.clusters = []
        self.clusters_plus = []
        self.clusters_minus = []
        self.clusterCount = 0
        self.currentCluster_plus = Alignment.BED("", 0, 0, "", 0, "+")
        self.currentCluster_minus = Alignment.BED("", 0, 0, "", 0, "-")
        # self.sigClusters = {}
        self.mutations = {}  # Dictionary of bed instance
        self.mutationCount = 0
        self.mutationList = []
        self.sigMutations = {}  # same as sigClusters
        self.kmpair = {}
        self.sigMutationCount = 0
        self.sigClusterCount = 0
        self.coverage = 0  # "reads coverage of this sample"
        self.bamheader = None
        self.crosslinking = {}
        self.crosslinkingMutations = []

    def __str__(self):
        pass

    def testInput(self):
        """Test the input file format, modify self.filepath and return bool
        Status: True:file is ready to use; False: wrong file, program stop
        """
        # test if file has header
        try:
            self.header = pysam.view("-H", self.filepath)
        except:
            try:
                self.header = pysam.view("-SH", self.filepath)
            except:
                LOGGER.error(
                    "Input file does not have header, please check your file. Program quit"
                )
                return (False, "None")
        # Header test passed, test if it is BAM
        try:
            infile = gzip.open(self.filepath)
            infile.readline(10)
        except:  # cannot read line, should be sam
            LOGGER.info("Input is SAM, converting to BAM...")
            bamout = ".".join(self.filepath.split(".")[0:-1]) + "." + "bam"
            infile = pysam.Samfile(self.filepath, "r", header=self.header)
            # print >> sys.stderr,pysam.view("-SH",infile)
            outfile = pysam.Samfile(bamout, "wb", template=infile)
            for i in infile.fetch():
                outfile.write(i)
            self.filepath = bamout
        # Now the infile is BAM,check if it is sorted
        if Utils.is_sorted(self.header):
            pysam.index(self.filepath)
            return True
        else:  # sort the BAM
            1 / 0
            LOGGER.info("Input is not sorted, sorting file...")
            bamsort = ".".join(self.filepath.split(".")[0:-1]) + "." + "sort.bam"
            pysam.sort("-o", bamsort, self.filepath)
            pysam.index(bamsort)
            self.filepath = bamsort  # change input file path
            self.header = pysam.view("-H", bamsort)
            LOGGER.info("Input file sorted")
            # if Utils.is_sorted(self.header):
            # 	print >> sys.stderr, "The file is sorted"
            return True

    def readfile(self):
        try:
            self.originalBAM = pysam.Samfile(self.filepath, "rb")
            self.refInfo = list(
                zip(self.originalBAM.references, self.originalBAM.lengths)
            )
            return True
        except IOError:
            LOGGER.error("Cannot open input file")
            return False

    # 	def printFilteredReads(self):
    # 		for i in self.filteredAlignment:
    # 			print i

    def printClusters(self):
        outfile = open(self.outprefix + ".clusters.bed", "w")
        for i in [self.clusters_plus, self.clusters_minus]:
            for j in i:
                print(j, file=outfile)
        outfile.close()

    def printMutations(self):
        outfile = open(self.outprefix + ".mutations.bed", "w")
        for i in list(self.mutations.values()):
            print(i, file=outfile)
        outfile.close()

    def printMutations_chr(self):
        fhs = {}  # file handle dic, key is chr, value is file handle
        for chr, chrlen in self.refInfo:
            outfile = open(self.outprefix + "." + chr + ".mutations.bed", "w")
            fhs[chr] = outfile
        for i in self.mutations.values():
            print(i, file=fhs[i.chr])
        # clean up
        for fh in fhs.values():
            fh.close()
        del fhs

    def printMutationSingleChr(self, chr):
        outfile = open(self.outprefix + "." + chr + ".mutations.bed", "w")
        for i in self.mutations.values():
            print(i, file=outfile)
        # clean up
        outfile.close()

    def printReliableMutations(self):
        outfile = open(self.outprefix + ".reliableMutations.pipeclip.bed", "w")
        header = "#chr\tstart\tstop\tmutation_name\tM_value\tstrand\ttype\tK_value\tp_value\tfdr"
        print(header, file=outfile)
        self.printEnrichedItem(self.sigMutations, outfile)

    def printEnrichedClusters(self):
        outfile = open(self.outprefix + ".enrichedClusters.pipeclip.bed", "w")
        header = "#chr\tstart\tstop\tcluster_name\tread_count\tstrand\tp_value\tfdr"
        print(header, file=outfile)
        self.printReliableList(self.clusters, outfile)
        return [self.outprefix + ".enrichedClusters.pipeclip.bed"]

    def printCrosslinkingMutations(self):
        outfile = open(self.outprefix + ".crosslinkingMutations.pipeclip.bed", "w")
        header = "#chr\tstart\tstop\tmutation_name\tM_value\tstrand\ttype\tK_value\tp_value\tfdr"
        print(header, file=outfile)
        self.printReliableList(self.crosslinkingMutations, outfile)

    def printReliableList(self, mylist, fh):
        for i in mylist:
            if i.sig:
                st = i.__str__()
                st += "\t" + str(i.pvalue) + "\t" + str(i.qvalue)
                print(st, file=fh)

    def printList(self, mylist, fh):
        for i in mylist:
            print(i, file=fh)

    def printEnrichedItem(self, dic, fh):
        for k in list(dic.keys()):
            for i in dic[k]:
                st = i.__str__()
                st += "\t" + str(i.pvalue) + "\t" + str(i.qvalue)
                print(st, file=fh)

    def printCrosslinkingSites(self):
        output = self.outprefix
        header = "#chr\tstart\tstop\tcluster_name\treads_count\tstrand\tcluster_fdr\tcrosslinking_fisherP\tmutation_pos\tmutation_name"
        if self.type == 0:  # HITS-CLIP, three output
            output_del = open(output + "_deletion_crosslinking.pipeclip.txt", "w")
            output_sub = open(output + "_substitution_crosslinking.pipeclip.txt", "w")
            output_ins = open(output + "_insertion_crosslinking.pipeclip.txt", "w")
            print(header, file=output_del)
            print(header, file=output_sub)
            print(header, file=output_ins)
        else:
            output_name = open(output + "_crosslinking.pipeclip.txt", "w")
            print(header, file=output_name)
        for k in list(self.crosslinking.keys()):
            st = self.crosslinking[k].__str__()
            st += "\t" + "\t".join(
                [
                    str(self.crosslinking[k].qvalue),
                    str(self.crosslinking[k].fisherP),
                ]
            )
            st += "\t" + ",".join(self.crosslinking[k].mutationStarts)
            st += "\t" + ",".join(self.crosslinking[k].mutationNames)
            if self.type == 0:
                output_key = k.split("_")[-1]
                if output_key == "Deletion":
                    print(st, file=output_del)
                elif output_key == "Insertion":
                    print(st, file=output_ins)
                elif output_key == "Substitution":
                    print(st, file=output_sub)
            else:
                print(st, file=output_name)

        if self.type == 0:
            output_del.close()
            output_sub.close()
            output_ins.close()
            return [
                output + "_insertion_crosslinking.pipeclip.txt",
                output + "_deletion_crosslinking.pipeclip.txt",
                output + "_substitution_crosslinking.pipeclip.txt",
            ]
        else:
            output_name.close()
            return [output + "_crosslinking.pipeclip.txt"]

    def updatePreviousQul(self, n, q, m):
        self.previousQul[0] = n
        self.previousQul[1] = q
        self.previousQul[2] = m

    def updateCurrentGroup(self, read, mlen, mis):
        """Compare read to current duplication group parameters, determine to
        add to, to drop or to replace, make sure duplication group only has
        reads with best quality."""
        if (
            mlen >= self.previousQul[0]
            and read.mapq >= self.previousQul[1]
            and mis <= self.previousQul[2]
        ):  # consider to append or replace only when read has no worse quality
            if (
                mlen > self.previousQul[0]
                or read.mapq > self.previousQul[1]
                or mis < self.previousQul[2]
            ):  # read has better quality,replace
                self.currentGroup = [read]
                self.updatePreviousQul(mlen, read.mapq, mis)
            else:
                self.currentGroup.append(read)

    def iniDupGroupInfo(self, read, group_key, mlength, mismatch):
        self.currentGroupKey = group_key
        self.currentGroup = [read]
        self.updatePreviousQul(mlength, read.mapq, mismatch)

    def rmdup(self):
        """Return random one read of highest quality from list."""
        # print "get one from group"
        if len(self.currentGroup) == 1:
            # print self.currentGroup[0]
            return self.currentGroup[0]
        else:
            random.seed(1)
            index = random.randint(0, len(self.currentGroup) - 1)
            # print self.currentGroup[index]
            return self.currentGroup[index]

    def updateCluster(self, read):
        """Cluster new read to known clusters and update cluster reads
        count."""
        strandDic = {"True": "-", "False": "+"}
        clusterName = (
            "cluster"
            + "_"
            + str(len(self.clusters_plus) + len(self.clusters_minus) + 1)
        )
        newRead = Alignment.ClusterBed(
            self.originalBAM.getrname(read.tid),
            read.pos,
            read.pos + len(read.seq),
            clusterName,
            1,
            strandDic[str(read.is_reverse)],
        )
        if read.is_reverse:
            if self.currentCluster_minus.chr == "":  # Initiate cluster
                self.currentCluster_minus = newRead
                self.clusters_minus.append(self.currentCluster_minus)
                self.clusterCount += 1
            else:
                if self.currentCluster_minus.overlap(newRead):
                    self.currentCluster_minus.merge(newRead)
                    self.clusters_minus[-1] = self.currentCluster_minus
                else:  # New read is a new cluster
                    # self.clusters.append(self.currentCluster)
                    self.currentCluster_minus = newRead
                    self.clusters_minus.append(self.currentCluster_minus)
                    self.clusterCount += 1
        else:
            if self.currentCluster_plus.chr == "":  # Initiate cluster
                self.currentCluster_plus = newRead
                self.clusters_plus.append(self.currentCluster_plus)
                self.clusterCount += 1
            else:
                if self.currentCluster_plus.overlap(newRead):
                    self.currentCluster_plus.merge(newRead)
                    self.clusters_plus[-1] = self.currentCluster_plus
                else:  # New read is a new cluster
                    # self.clusters.append(self.currentCluster)
                    self.currentCluster_plus = newRead
                    self.clusters_plus.append(self.currentCluster_plus)
                    self.clusterCount += 1

    def updateMutation(self, read, mis):
        mutations = []
        if self.type == 3:  # iclip,find truncation
            mutations = Mutation2.getTruncations(self.originalBAM, read)
        else:
            mutations = Mutation2.getMutations(self.originalBAM, read)
            if self.type == 1:
                mutation_filter = Utils.filterMutations(mutations, "T->C", True)
                mutations = mutation_filter
            elif self.type == 2:
                mutation_filter = Utils.filterMutations(mutations, "G->A", True)
                mutations = mutation_filter
        # LOGGER.debug("read %s " % read)
        if len(mutations) > 0:
            for m in mutations:
                # print m
                self.mutationCount += 1
                m_key = "_".join([m.chr, str(m.start), m.strand, m.type])
                if m_key in self.mutations:
                    self.mutations[m_key].increaseScore()
                else:
                    self.mutations[m_key] = m

    def updateCLIPinfo(self, read, matchlen, miscount):
        """Update sample coverage info, clustering, mutation info."""
        # LOGGER.debug("read %s, cigar %s,mismatch %d" % (read.qname,read.cigar,miscount))
        # update sample coverage info
        self.coverage += matchlen
        # update cluster info
        self.updateCluster(read)
        # update mutation info
        if miscount > 0:
            self.updateMutation(read, miscount)

    def addSigToDic(self, dic, mu):
        """Add new mutation into the dictionary.Mutations should be sorted."""
        if mu.chr in dic:
            dic[mu.chr].append(mu)
        else:
            dic[mu.chr] = [mu]

    def getCrosslinking(self):
        """Merge enriched clusters and reliable mutations together Call
        Enrich.fisherTest() to calculate joint p vlaue."""
        for cluster in self.clusters:
            # LOGGER.debug("P value of cluster is %f" % cluster.pvalue)
            if cluster.sig and cluster.chr in self.sigMutations:
                for mutation in self.sigMutations[cluster.chr]:
                    # LOGGER.debug("Cluster loc %d,%d ; Mutation loc %d,%d, mutation type %s" % (cluster.start,cluster.stop,mutation.start,mutation.stop,mutation.type))
                    if cluster.overlap(mutation):
                        # LOGGER.debug("Overlapped")
                        if self.type == 0:  # HITS-CLIP
                            mutation_key = mutation.type.split("->")[0]
                            if mutation_key in ["A", "C", "G", "T"]:
                                mutation_key = "Substitution"
                            cross_key = cluster.name + "_" + mutation_key
                        else:
                            cross_key = cluster.name
                        if cross_key in self.crosslinking:
                            # LOGGER.debug("Existing mutation pvalue:",self.crosslinking[cluster.name].mutationP)
                            self.crosslinking[cross_key].addMutation(mutation)
                            self.crosslinkingMutations.append(mutation)
                        else:
                            self.crosslinking[cross_key] = Alignment.CrosslinkingBed(
                                cluster.chr,
                                cluster.start,
                                cluster.stop,
                                cluster.name,
                                cluster.score,
                                cluster.strand,
                                cluster.pvalue,
                                cluster.qvalue,
                                mutation.start,
                                mutation.name,
                                mutation.pvalue,
                            )
                            self.crosslinkingMutations.append(mutation)
        # start to calculate fisher test p value
        for k in self.crosslinking.keys():
            self.crosslinking[k].fishertest()

    def filter(self, matchLen, mismatch, cliptype, duprm):
        """Filter the input BAM file according to parameters.

        Make clusters and mutations at the same time
        """
        LOGGER.info("Start to filter alignment using parameters:")
        LOGGER.info("match length:%d" % (matchLen))
        LOGGER.info("mismatch count: %d" % (mismatch))
        LOGGER.info("CLIP type:%s" % (str(cliptype)))
        LOGGER.info("Rmdup code:%s" % (str(duprm)))
        LOGGER.info(
            "There are %d reads in origianl input file" % (self.originalBAM.mapped)
        )
        # print >> sys.stderr,zip(self.originalBAM.references,self.originalBAM.lengths)
        # sys_info = os.system("free -m")
        # LOGGER.debug(sys_info)
        # outBAM = pysam.Samfile(outprefix+".filtered.bam","wb",template=self.originalBAM)
        self.originalMapped = self.originalBAM.mapped
        outBAM_pos = pysam.Samfile(
            self.outprefix + ".pos.filtered.bam",
            "wb",
            template=self.originalBAM,
        )
        outBAM_neg = pysam.Samfile(
            self.outprefix + ".neg.filtered.bam",
            "wb",
            template=self.originalBAM,
        )
        self.type = cliptype
        if cliptype == 3:  # make sure there is no rmdup for iCLIP data
            duprm = 0
        count = 0
        start_time = datetime.datetime.now()
        for alignment in self.originalBAM:
            # print "Now processing",alignment.qname
            if not alignment.cigar:  # reads is unmapped
                continue
            count += 1
            if count % 500000 == 0:
                stop_time = datetime.datetime.now()
                LOGGER.debug(
                    "Processed %d reads in %s" % (count, str(stop_time - start_time))
                )
                start_time = stop_time
            flag, mlen, mis = Utils.readQuaFilter(alignment, matchLen, mismatch)
            if flag:
                # print "Qualified read"
                # print	alignment
                # print "current Gourp key",self.currentGroupKey
                if duprm > 0:
                    # get group key
                    if duprm == 1:
                        groupkey = Utils.rmdupKey_Start(alignment)
                    elif duprm == 2:
                        groupkey = Utils.rmdupKey_Seq(alignment)
                    # check current group
                    if (
                        groupkey == self.currentGroupKey
                    ):  # overlap with current group, update group
                        self.updateCurrentGroup(alignment, mlen, mis)
                    else:  # get read from current group and discard it, use current read to start a new current group
                        if self.currentGroupKey != "None":  # current group exists
                            keep = self.rmdup()
                            # LOGGER.debug("Pop out read to keep %s" % keep)
                            self.currentGroup = []
                            self.filteredAlignment += 1
                            flag, mlen, mis = Utils.readQuaFilter(
                                keep, matchLen, mismatch
                            )
                            self.updateCLIPinfo(keep, mlen, mis)
                            # outBAM.write(keep)
                            if keep.is_reverse:
                                outBAM_neg.write(keep)
                            else:
                                outBAM_pos.write(keep)
                        self.iniDupGroupInfo(
                            alignment, groupkey, mlen, mis
                        )  # make new group using current alignment
                else:  # there is no rmdup
                    # LOGGER.debug("Good read, update clip info %s" % read.qname)
                    self.filteredAlignment += 1
                    self.updateCLIPinfo(alignment, mlen, mis)
                    # outBAM.write(alignment)
                    if alignment.is_reverse:
                        outBAM_neg.write(alignment)
                    else:
                        outBAM_pos.write(alignment)
        # clean up the final dupGroup, if rmdup==0, there is no final dupGroup
        if len(self.currentGroup) > 0:
            keep = self.rmdup()
            self.currentGroup = []
            self.filteredAlignment += 1
            flag, mlen, mis = Utils.readQuaFilter(keep, matchLen, mismatch)
            self.updateCLIPinfo(keep, mlen, mis)
            # outBAM.write(alignment)
            if keep.is_reverse:
                outBAM_neg.write(keep)
            else:
                outBAM_pos.write(keep)
        # LOGGER CLIP sample information
        # outBAM.close()
        outBAM_pos.close()
        outBAM_neg.close()
        # pysam.index(outprefix+".filtered.bam")
        pysam.index(self.outprefix + ".pos.filtered.bam")
        pysam.index(self.outprefix + ".neg.filtered.bam")
        self.posfilteredBAM = self.outprefix + ".pos.filtered.bam"
        self.negfilteredBAM = self.outprefix + ".neg.filtered.bam"
        self.clusters = self.clusters_plus + self.clusters_minus
        self.printClusters()
        self.printMutations_chr()
        # Clean up variables
        self.originalBAM = None
        self.clusters_plus = None
        self.clusters_minus = None
        self.mutations = None
        gc.collect()
        LOGGER.debug("After filtering, %d reads left" % (self.filteredAlignment))
        LOGGER.debug("There are %d clusters in total" % (self.clusterCount))
        LOGGER.debug("There are %d mutations in total" % (self.mutationCount))

    def filter2(self, matchLen, mismatch, cliptype, duprm):
        """Filter the input BAM file according to parameters.

        Make clusters and mutations at the same time
        """
        LOGGER.info("Start to filter alignment using parameters:")
        LOGGER.info("match length:%d" % (matchLen))
        LOGGER.info("mismatch count: %d" % (mismatch))
        LOGGER.info("CLIP type:%s" % (str(cliptype)))
        LOGGER.info("Rmdup code:%s" % (str(duprm)))
        LOGGER.info(
            "There are %d reads in origianl input file" % (self.originalBAM.mapped)
        )

        self.originalMapped = self.originalBAM.mapped
        outBAM_pos = pysam.Samfile(
            self.outprefix + ".pos.filtered.bam",
            "wb",
            template=self.originalBAM,
        )
        outBAM_neg = pysam.Samfile(
            self.outprefix + ".neg.filtered.bam",
            "wb",
            template=self.originalBAM,
        )
        self.type = cliptype
        if cliptype == 3:  # make sure there is no rmdup for iCLIP data
            duprm = 0
        count = 0
        start_time = datetime.datetime.now()
        currentChr = ""
        for alignment in self.originalBAM:
            if not alignment.cigar:  # reads is unmapped
                continue
            c_chr = self.originalBAM.getrname(alignment.tid)
            if currentChr == "":
                currentChr = c_chr
            elif c_chr != currentChr:
                # LOGGER.debug("Mutation sites for Current ref %s is %d" % (currentChr,len(self.mutations.keys())))
                self.printMutationSingleChr(currentChr)
                self.mutations = None
                gc.collect()
                self.mutations = {}
                currentChr = c_chr
            count += 1
            if count % 5000000 == 0:
                stop_time = datetime.datetime.now()
                LOGGER.debug(
                    "Processed %d reads in %s" % (count, str(stop_time - start_time))
                )
                start_time = stop_time
            flag, mlen, mis = Utils.readQuaFilter(alignment, matchLen, mismatch)
            if flag:
                if duprm > 0:
                    # get group key
                    if duprm == 1:
                        groupkey = Utils.rmdupKey_Start(alignment)
                    elif duprm == 2:
                        groupkey = Utils.rmdupKey_Seq(alignment)
                    # check current group
                    if (
                        groupkey == self.currentGroupKey
                    ):  # overlap with current group, update group
                        self.updateCurrentGroup(alignment, mlen, mis)
                    else:  # get read from current group and discard it, use current read to start a new current group
                        if self.currentGroupKey != "None":  # current group exists
                            keep = self.rmdup()
                            # LOGGER.debug("Pop out read to keep %s" % keep)
                            self.currentGroup = []
                            self.filteredAlignment += 1
                            flag, mlen, mis = Utils.readQuaFilter(
                                keep, matchLen, mismatch
                            )
                            self.updateCLIPinfo(keep, mlen, mis)
                            # outBAM.write(keep)
                            if keep.is_reverse:
                                outBAM_neg.write(keep)
                            else:
                                outBAM_pos.write(keep)
                        self.iniDupGroupInfo(
                            alignment, groupkey, mlen, mis
                        )  # make new group using current alignment
                else:  # there is no rmdup
                    # LOGGER.debug("Good read, update clip info %s" % read.qname)
                    self.filteredAlignment += 1
                    self.updateCLIPinfo(alignment, mlen, mis)
                    # outBAM.write(alignment)
                    if alignment.is_reverse:
                        outBAM_neg.write(alignment)
                    else:
                        outBAM_pos.write(alignment)
        # clean up the final dupGroup, if rmdup==0, there is no final dupGroup
        if len(self.currentGroup) > 0:
            keep = self.rmdup()
            self.currentGroup = []
            self.filteredAlignment += 1
            flag, mlen, mis = Utils.readQuaFilter(keep, matchLen, mismatch)
            self.updateCLIPinfo(keep, mlen, mis)
            # outBAM.write(alignment)
            if keep.is_reverse:
                outBAM_neg.write(keep)
            else:
                outBAM_pos.write(keep)
        # clean up current mutation hash
        self.printMutationSingleChr(currentChr)
        self.mutations = None
        # LOGGER CLIP sample information
        # outBAM.close()
        outBAM_pos.close()
        outBAM_neg.close()
        pysam.index(self.outprefix + ".pos.filtered.bam")
        pysam.index(self.outprefix + ".neg.filtered.bam")
        self.posfilteredBAM = self.outprefix + ".pos.filtered.bam"
        self.negfilteredBAM = self.outprefix + ".neg.filtered.bam"
        self.clusters = self.clusters_plus + self.clusters_minus
        self.printClusters()
        # Clean up variables
        self.originalBAM = None
        self.clusters_plus = None
        self.clusters_minus = None
        self.mutations = None
        gc.collect()
        LOGGER.debug("After filtering, %d reads left" % (self.filteredAlignment))
        LOGGER.debug("There are %d clusters in total" % (self.clusterCount))
        LOGGER.debug("There are %d mutations in total" % (self.mutationCount))
