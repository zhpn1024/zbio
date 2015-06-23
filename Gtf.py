
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
  def __len__(self): 
    return self.stop - self.start
  @property
  def length(self):
    return len(self)
  def is_reverse(self): 
    return self.strand=='-'
  @property
  def end5(self): #5' end, all bed
    if self.strand != '-' : return self.start
    else : return self.stop
  @property
  def end3(self): #3' end, all bed
    if self.strand != '-' : return self.stop
    else : return self.start
		
def attr(s):
	l = s.strip(';').split('; ')
	a = {}
	#print l
	for att in l:
		l2 = att.split(' ')
		#print att, l2
		a[l2[0]] = eval(l2[1])
	return a
	
class gtfTrans(exon):
  def __init__(self, lst): 
    exon.__init__(self, lst)
    #self.pos = exon(lst)
    self.id = self.attr['transcript_id']
    #self.chr = lst[0]
    self.type = lst[1]
    self.exons = []
    self.cds = []
    self.utr = []
    self.start_codon = None
    self.stop_codon = None
    self.other = []
    self.gene = self.attr['gene_id']
    #self.attr = self.pos.attr
  def addExon(self, e): 
    if e.type == 'exon': self.exons.append(e)
    elif e.type == 'CDS': self.cds.append(e)
    elif e.type == 'UTR': self.utr.append(e)
    elif e.type == 'start_codon': self.start_codon = e
    elif e.type == 'stop_codon': self.stop_codon = e
    else: self.other.append(e)
    self.check()
  def __repr__(self):
    s = "transcript_id " + self.id + ', ' + str(len(self.exons)) + " exons, " + exon.__repr__(self)
    #s += self.pos.__repr__()
    return s
  def check(self):
    for e in self.exons:
      if e.start < self.start: self.start = e.start
      if e.stop > self.stop: self.stop = e.stop
    self.exons.sort(reverse = self.is_reverse())
    self.cds.sort(reverse = self.is_reverse())
  @property
  def cds_start(self): #Bed12
    return self.start_codon.end5
  @property
  def cds_stop(self):
    return self.stop_codon.end3
  def cdna_length(self): 
    l = 0
    for e in self.exons:
      l += len(e)
    return l
  
  def cdna_pos(self, p):
    #self.check()
    if p < self.start or p > self.stop:
      return None
    #p1 = p - self.start
    pos = 0
    for e in self.exons:
      if e.start <= p <= e.stop:
        if self.is_reverse():
          pos += e.stop - p
        else:
          pos += p - e.start
        return pos
      else:
        pos += len(e)
    return None

  def genome_pos(self, p, bias = 0):
    m = self.cdna_length()
    if p < 0 or p > m: return None
    if p == 0 : return self.end5
    if p == m: return self.end3
    p1 = p
    pos = self.start
    #l=range(len())
    for e in self.exons:
      if len(e) - p1 >= bias:
        if self.is_reverse():
          pos = e.stop - p1
        else:
          pos = p1 + e.start
        return pos
      else:
        p1 -= len(e)
    return None
    
class gtfGene(exon):
  def __init__(self, lst):
    exon.__init__(self, lst)
    #self.pos = exon(lst)
    self.id = self.attr['gene_id']
    self.trans = []
    self.type = lst[1]
    #self.attr = self.pos.attr
  def addTrans(self, tr):
    #tr = gtf_t(lst)
    self.trans.append(tr)
    self.check()
  def __repr__(self):
    s = "gene_id " + self.id + ', ' + str(len(self.trans)) + " transcripts, " + exon.__repr__(self)
    #s += self.pos.__repr__()
    for t in self.trans:
      s += '\n\t' + t.__repr__()
    return s
  def check(self):
    for t in self.trans:
      t.check()
      if t.start < self.start: self.start = t.start
      if t.stop > self.stop: self.stop = t.stop
    
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
      if t.attr['gene_id'] not in genes:
        g = gtfGene(lst)
        genes[g.id] = g
        #continue ##
      genes[t.attr['gene_id']].addTrans(t)
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
        t = gtfTrans(lst)
        genes[t.attr['gene_id']].addTrans(t)
        trans[t.id] = t
      trans[e.attr['transcript_id']].addExon(e)
  return genes, trans

def gtfGeneIter(fin, filt = []):
  genes = {}
  trans = {}
  gidlist = []
  chr = ""
  for l in fin:
    if l[0] == '#' : continue
    lst=l.strip().split('\t')
    if lst[0] != chr:
      for gid in gidlist:
        yield genes[gid]
      genes = {}
      trans = {}
      gidlist = []
      chr = lst[0]
    if lst[2] == 'gene':
      if filt != [] and lst[1] not in filt: 
        continue
      g = gtfGene(lst)
      genes[g.id] = g
      if g.id not in gidlist: gidlist.append(g.id)
    elif lst[2] == 'transcript':
      t = gtfTrans(lst)
      if t.attr['gene_id'] not in genes:
        g = gtfGene(lst)
        genes[g.id] = g
        gidlist.append(g.id)
      genes[t.attr['gene_id']].addTrans(t)
      trans[t.id] = t
    else:
      e = exon(lst)
      if e.attr['gene_id'] not in genes:
        g = gtfGene(lst)
        genes[g.id] = g
        gidlist.append(g.id)
      if e.attr['transcript_id'] not in trans:
        t = gtfTrans(lst)
        genes[t.attr['gene_id']].addTrans(t)
        trans[t.id] = t
      trans[e.attr['transcript_id']].addExon(e)
  for gid in gidlist:
    yield genes[gid]
    
def gtfTransIter(fin, filt = []):
  genes = {}
  trans = {}
  gidlist = []
  chr = ""
  for l in fin:
    if l[0] == '#' : continue
    lst=l.strip().split('\t')
    if lst[0] != chr:
      for gid in gidlist:
        for t in genes[gid].trans:
          yield t
      genes = {}
      trans = {}
      gidlist = []
      chr = lst[0]
    if lst[2] == 'gene':
      if filt != [] and lst[1] not in filt: 
        continue
      g = gtfGene(lst)
      genes[g.id] = g
      if g.id not in gidlist: gidlist.append(g.id)
    elif lst[2] == 'transcript':
      t = gtfTrans(lst)
      if t.attr['gene_id'] not in genes:
        g = gtfGene(lst)
        genes[g.id] = g
        gidlist.append(g.id)
      genes[t.attr['gene_id']].addTrans(t)
      trans[t.id] = t
    else:
      e = exon(lst)
      if e.attr['gene_id'] not in genes:
        g = gtfGene(lst)
        genes[g.id] = g
        gidlist.append(g.id)
      if e.attr['transcript_id'] not in trans:
        t = gtfTrans(lst)
        genes[t.attr['gene_id']].addTrans(t)
        trans[t.id] = t
        #print e.attr['transcript_id'], t.id
      trans[e.attr['transcript_id']].addExon(e)
  for gid in gidlist:
    for t in genes[gid].trans:
      yield t
  #return genes, trans

# def gtfGeneIter(fin):
#   gid = ''
#   for l in fin:
#     if l[0] == '#' : continue
#     lst=l.strip().split('\t')
#     if lst[2] == 'gene':
#       if gid != '': yield g
#       g = gtfGene(lst)
#       gid = g.id
#       trans = {}
#     elif lst[2] == 'transcript':
#       t = gtfTrans(lst)
#       if t.attr['gene_id'] != gid:
#         if gid != '': yield g
#         g = gtfGene(lst)
#         gid = g.id
#         trans = {}
#       g.addTrans(t)
#       trans[t.id] = t
#     else:
#       e = exon(lst)
#       if e.attr['gene_id'] != gid:
#         if gid != '': yield g
#         g = gtfGene(lst)
#         gid = g.id
#         trans = {}
#       if e.attr['transcript_id'] not in trans:
#         t = gtfTrans(lst)
#         g.addTrans(t)
#         trans[t.id] = t
#       t.addExon(e)
#   if gid != '': yield g

# def gtfTransIter(fin):
#   tid = ''
#   for l in fin:
#     if l[0] == '#' : continue
#     lst=l.strip().split('\t')
#     if lst[2] == 'gene':
#       if tid != '': yield t
#     elif lst[2] == 'transcript':
#       if tid != '': yield t
#       t = gtfTrans(lst)
#       tid = t.id
#     else:
#       e = exon(lst)
#       if e.attr['transcript_id'] != tid:
#         if tid != '': yield t
#         t = gtfTrans(lst)
#         tid = t.id
#       t.addExon(e)
#   if tid != '': yield t

