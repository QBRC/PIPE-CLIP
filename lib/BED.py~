class BED:
		def __init__(self,chr,start,stop,name,score,strand):
			self.chr = chr
			self.start = start
			self.stop = stop
			self.name = name
			self.score = score
			self.strand = strand
	
		def __str__(self):
			st = "\t".join([self.chr,str(self.start),str(self.stop),self.name,str(self.score),self.strand])
			return st

		def merge(self,read):
				self.stop = read.stop
				self.score += 1

		def overlap(self,read):
				if self.chr == read.chr or self.strand == read.strand:
						if self.start <= read.stop and self.stop >=read.start:
								return True
				else:
						return False


