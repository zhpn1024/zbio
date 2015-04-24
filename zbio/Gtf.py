
class exon:
	def __init__(self, lst):
		self.chr = lst[0]
		self.strand = lst[6]
		self.start = long(lst[3]) - 1
		self.stop = long(lst[4])
		self.type = lst[2]
		self.attr = attr(lst[8])
	def __repr__(self):
		return self.chr + ':'+str(self.start)+'-'+str(self.stop)+':'+self.strand
		
def attr(s):
	l = s.strip(';').split('; ')
	a = {}
	#print l
	for att in l:
		l2 = att.split(' ')
		#print att, l2
		a[l2[0]] = eval(l2[1])
	return a
	
class gtf_t:
	def __init__(self, lst): 
		self.pos = exon(lst)
		self.id = self.pos.attr['transcript_id']
		self.type = lst[1]
		self.exons = []
		self.cds = []
		self.utr = []
		self.start_codon = None
		self.stop_codon = None
		self.other = []
	def add_exon(self, e): 
		#e = exon(lst)
		if e.type == 'exon': self.exons.append(e)
		elif e.type == 'CDS': self.cds.append(e)
		elif e.type == 'UTR': self.utr.append(e)
		elif e.type == 'start_codon': self.start_codon = e
		elif e.type == 'stop_codon': self.stop_codon = e
		else: self.other.append(e)
	def __repr__(self):
		s = "transcript_id " + self.id + ', ' + str(len(self.exons)) + " exons, " + self.pos.__repr__() + "\n"
		#s += self.pos.__repr__()
		return s
		
class gtf_g:
	#transcripts=[]
	def __init__(self, lst):
		self.pos = exon(lst)
		self.id = self.pos.attr['gene_id']
		self.trans = []
		self.type = lst[1]
	def add_trans(self, tr):
		#tr = gtf_t(lst)
		self.trans.append(tr)
	def __repr__(self):
		s = "gene_id " + self.id + ', ' + str(len(self.trans)) + " transcripts, " + self.pos.__repr__() + "\n"
		#s += self.pos.__repr__()
		for t in self.trans:
			s += '\t' + t.__repr__()
		return s
	
def load_gtf(file, filt = []):
	genes = {}
	trans = {}
	for l in file:
		if l[0] == '#' : continue
		lst=l.strip().split('\t')
		
		if lst[2] == 'gene':
			if filt != [] and lst[1] not in filt: 
				continue
			g = gtf_g(lst)
			genes[g.id] = g
			#print g.id
		elif lst[2] == 'transcript':
			t = gtf_t(lst)
			#g.add_trans(t)
			if t.pos.attr['gene_id'] not in genes:
				continue
			genes[t.pos.attr['gene_id']].add_trans(t)
			trans[t.id] = t
			#print t.id
		else:
			e = exon(lst)
			if e.attr['gene_id'] not in genes:
				continue
			#print e.attr['exon_number']
			trans[e.attr['transcript_id']].add_exon(e)
			
	return genes, trans
		
		