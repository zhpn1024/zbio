
class exon:
  def __init__(self, lst):
    self.chr = lst[0]
    self.strand = lst[6]
    self.start = int(lst[3]) - 1
    self.stop = int(lst[4])
    self.type = lst[2]
    self.attr = attr(lst[8])
  def __repr__(self):
    return self.chr + ':'+str(self.start)+'-'+str(self.stop)+':'+self.strand
  def __str__(self):
    return "%s:%d-%d:%s" % (self.chr, self.start, self.stop, self.strand)
  def __cmp__(self, other):
    return cmp(self.chr,other.chr) or cmp(self.start,other.start) or cmp(self.stop,other.stop)
		
def attr(s):
	l = s.strip(';').split('; ')
	a = {}
	#print l
	for att in l:
		l2 = att.split(' ')
		#print att, l2
		a[l2[0]] = eval(l2[1])
	return a
	
class gtfTrans:
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
    self.attr = self.pos.attr
  def addExon(self, e): 
    if e.type == 'exon': self.exons.append(e)
    elif e.type == 'CDS': self.cds.append(e)
    elif e.type == 'UTR': self.utr.append(e)
    elif e.type == 'start_codon': self.start_codon = e
    elif e.type == 'stop_codon': self.stop_codon = e
    else: self.other.append(e)
  def __repr__(self):
    s = "transcript_id " + self.id + ', ' + str(len(self.exons)) + " exons, " + self.pos.__repr__()
    #s += self.pos.__repr__()
    return s
  def check(self):
    for e in self.exons:
      if e.start < self.pos.start: self.pos.start = e.start
      if e.stop > self.pos.stop: self.pos.stop = e.stop
    
class gtfGene:
  def __init__(self, lst):
    self.pos = exon(lst)
    self.id = self.pos.attr['gene_id']
    self.trans = []
    self.type = lst[1]
    self.attr = self.pos.attr
  def addTrans(self, tr):
    #tr = gtf_t(lst)
    self.trans.append(tr)
  def __repr__(self):
    s = "gene_id " + self.id + ', ' + str(len(self.trans)) + " transcripts, " + self.pos.__repr__()
    #s += self.pos.__repr__()
    for t in self.trans:
      s += '\n\t' + t.__repr__()
    return s
  def check(self):
    for t in self.trans:
      t.check()
      if t.start < self.pos.start: self.pos.start = t.start
      if t.stop > self.pos.stop: self.pos.stop = t.stop
    
def loadGtf(fin, filt = []):
  genes = {}
  trans = {}
  for l in fin:
    if l[0] == '#' : continue
    lst=l.strip().split('\t')
    if lst[2] == 'gene':
      if filt != [] and lst[1] not in filt: 
        continue
      g = gtfGene(lst)
      genes[g.id] = g
      #print g.id
    elif lst[2] == 'transcript':
      t = gtfTrans(lst)
      #g.add_trans(t)
      if t.pos.attr['gene_id'] not in genes:
        g = gtfGene(lst)
        genes[g.id] = g
        #continue ##
      genes[t.pos.attr['gene_id']].addTrans(t)
      trans[t.id] = t
      #print t.id
    else:
      e = exon(lst)
      if e.attr['gene_id'] not in genes:
        g = gtfGene(lst)
        genes[g.id] = g
        #continue
      #print e.attr['exon_number']
      if e.attr['transcript_id'] not in trans:
        t = gtfGene(lst)
        genes[t.pos.attr['gene_id']].addTrans(t)
        trans[t.id] = t
      trans[e.attr['transcript_id']].addExon(e)
  return genes, trans

def gtfGeneIter(fin):
  gid = ''
  for l in fin:
    if l[0] == '#' : continue
    lst=l.strip().split('\t')
    if lst[2] == 'gene':
      if gid != '': yield g
      g = gtfGene(lst)
      gid = g.id
      trans = {}
    elif lst[2] == 'transcript':
      t = gtfTrans(lst)
      if t.pos.attr['gene_id'] != gid:
        if gid != '': yield g
        g = gtfGene(lst)
        gid = g.id
        trans = {}
      g.addTrans(t)
      trans[t.id] = t
    else:
      e = exon(lst)
      if e.attr['gene_id'] != gid:
        if gid != '': yield g
        g = gtfGene(lst)
        gid = g.id
        trans = {}
      if e.attr['transcript_id'] not in trans:
        t = gtfTrans(lst)
        g.addTrans(t)
        trans[t.id] = t
      t.addExon(e)
  if gid != '': yield g

def gtfTransIter(fin):
  tid = ''
  for l in fin:
    if l[0] == '#' : continue
    lst=l.strip().split('\t')
    if lst[2] == 'gene':
      if tid != '': yield t
    elif lst[2] == 'transcript':
      if tid != '': yield t
      t = gtfTrans(lst)
      tid = t.id
    else:
      e = exon(lst)
      if e.attr['transcript_id'] != tid:
        if tid != '': yield t
        t = gtfTrans(lst)
        tid = t.id
      t.addExon(e)
  if tid != '': yield t

