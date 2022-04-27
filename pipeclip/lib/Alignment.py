#!/usr/bin/python
# programmer: beibei.chen@utsouthwestern.edu
# Usage: definition of alignment files including BED and BAM


from . import Enrich


class BED:
    def __init__(self, chr, start, stop, name, score, strand):
        self.chr = chr
        self.start = int(start)
        self.stop = int(stop)
        self.name = name
        self.score = int(score)
        self.strand = strand

    def __str__(self):
        st = "\t".join(
            [
                str(self.chr),
                str(self.start),
                str(self.stop),
                self.name,
                str(self.score),
                self.strand,
            ]
        )
        return st

    def merge(self, read):
        self.stop = read.stop
        self.score += 1

    def overlap(self, read):
        if self.chr == read.chr and self.strand == read.strand:
            if self.start <= read.stop and self.stop >= read.start:
                return True
            else:
                return False

    def updateScore(self, s):
        self.score = int(s)

    def increaseScore(self):
        self.score += 1


class ClusterBed(BED):
    def __init__(self, chr, start, stop, name, score, strand):
        self.pvalue = 0
        self.qvalue = 0
        self.sig = False
        BED.__init__(self, chr, start, stop, name, score, strand)


class MutationBed(BED):
    def __init__(self, chr, start, stop, name, score, strand, type):
        self.type = type  # insertion,deletion,type of substitution
        self.kvalue = 0
        BED.__init__(self, chr, start, stop, name, score, strand)
        self.pvalue = 0
        self.qvalue = 0
        self.sig = False

    def updateK(self, k):
        self.kvalue = k

    def __str__(self):
        st = BED.__str__(self)
        st += "\t" + self.type + "\t" + str(self.kvalue)
        return st


class CrosslinkingBed(BED):
    def __init__(
        self,
        chr,
        start,
        stop,
        name,
        score,
        strand,
        clusterP,
        clusterQ,
        mStart,
        mName,
        mP,
    ):
        self.fisherP = 0
        self.clusterP = float(clusterP)
        self.mutationP = [mP]
        self.qvalue = float(clusterQ)
        self.mutationStarts = [str(mStart)]
        self.mutationNames = [mName]
        BED.__init__(self, chr, start, stop, name, score, strand)

    def addMutation(self, mu):  # mu is a instance
        self.mutationNames.append(mu.name)
        self.mutationStarts.append(str(mu.start))
        self.mutationP.append(mu.pvalue)

    def fishertest(self):
        fp = Enrich.fisherTest(self.clusterP, self.mutationP)
        self.fisherP = fp


class wiglist:
    def __init__(self):
        self.pos = []
        self.value = []

    def valueByPos(self, p):
        try:
            return self.value[self.pos.index(p)]
        except:
            return False

    def update(self, p, v):
        self.pos.append(p)
        self.value.append(v)
