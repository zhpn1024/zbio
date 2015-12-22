
class exon:
  def __init__(self, lst, gff = False, addchr = False):
    self.chr = lst[0]
    if addchr and self.chr[0:3] != 'chr' : self.chr = 'chr' + self.chr
    self.strand = lst[6]
    self.start = int(lst[3]) - 1
    self.stop = int(lst[4])
    self.type = lst[2]
    self.score = lst[5]
    self.attrstr = lst[8]
    self.attr = attr(lst[8], gff)
    #except : 
      #print "Format error:", self.attrstr
    self.gff = gff
    self.addchr = addchr
    self.frame = lst[7]
    self.lst = lst
  def __repr__(self):
    return self.chr + ':'+str(self.start)+'-'+str(self.stop)+':'+self.strand
  def __str__(self):
    return "%s\tunknown\t%s\t%d\t%d\t.\t%s\t%s\t%s" % (self.chr, self.type, self.start+1, self.stop, self.strand,self.frame,self.attrstr)
  def short_str(self):
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
  def gid(self):
    try : 
      if self.gff : return self.attr['GeneID']
      else : return self.attr['gene_id']
    except : return ''
  @property
  def tid(self):
    try : return self.attr['transcript_id']
    except : return ""
  @property
  def id(self):
    try: return self.attr['exon_id']
    except : return ""
  @property
  def end5(self): #5' end, all bed
    if self.strand != '-' : return self.start
    else : return self.stop
  @property
  def end3(self): #3' end, all bed
    if self.strand != '-' : return self.stop
    else : return self.start
  def short_str(self):
    return "%s:%d-%d:%s" % (self.chr, self.start, self.stop, self.strand)
  def __call__(self, **args):
    lst = self.lst[:]
    new = exon(lst, self.gff, self.addchr)
    for k in args.keys():
      new.__dict__[k]=args[k]
    return new
  def __sub__(self, other):
    return sub(self, other)
  def intersect(self, other):
    return intersect(self, other)
  def union(self, other):
    return union(self, other)
    
def attr(s, gff = False):
  if type(s) == dict: return s
  if gff : 
    l = s.strip(';').split(';')
    a = {}
    for att in l:
      l2 = att.split('=')
    #print att, l2
      a[l2[0]] = l2[1]
      if l2[0] == "Dbxref" :
        l3 = l2[1].split(',')
        for att2 in l3:
          l4 = att2.split(':')
          a[l4[0]] = l4[1]
    return a
  l = s.strip(';').split('; ')
  a = {}
  #print l
  for att in l:
    att = att.strip()
    l2 = att.split(' ')
    #print att, l2
    l21 = ' '.join(l2[1:])
    a[l2[0]] = eval(l21)
    #except : print l2[0]+'|'+l2[1]
  return a
  
class gtftrans(exon):
  def __init__(self, lst, gff = False, addchr = False): 
    exon.__init__(self, lst, gff, addchr)
    #self.pos = exon(lst)
    #self.id = self.attr['transcript_id']
    self.id = self.tid
    #self.chr = lst[0]
    self.type = lst[1]
    self.exons = []
    self.cds = []
    self.utr = []
    self.start_codon = None
    self.stop_codon = None
    self.other = []
    #self.gene = self.attr['gene_id']
    #self.attr = self.pos.attr
  def add_exon(self, e): 
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
  #@property
  def cds_start(self, cdna = False): 
    if self.start_codon is None : return None
    cs = self.cdna_pos(self.start_codon.end3) - 3
    if cdna : return cs
    else : return self.genome_pos(self.cdna_pos(self.start_codon.end3) - 3, 1)
  #@property
  def cds_stop(self, cdna = False):
    if self.stop_codon is None : return None
    cs = self.cdna_pos(self.stop_codon.end5) + 3
    if cdna : return cs
    else : return self.genome_pos(self.cdna_pos(self.stop_codon.end5) + 3, 0)
  @property
  def thick_start(self):
    if self.is_reverse(): return self.cds_stop()
    else : return self.cds_start()
  @property
  def thick_stop(self):
    if self.is_reverse(): return self.cds_start()
    else : return self.cds_stop()
  def cdna_length(self): 
    l = 0
    for e in self.exons:
      l += len(e)
    return l
  def cds_length(self): 
    try : return self.cdna_pos(self.stop_codon.end5) - self.cdna_pos(self.start_codon.end3) + 6 ##
    except: return 0
  @property
  def introns(self):
    introns = []
    last = -1
    for e in self.exons:
      if last >= 0 : 
        lst = [last, e.end5]
        lst.sort()
        introns.append(exon([self.chr,'','intron',lst[0]+1,lst[1],'',self.strand,'',self.attrstr], gff=self.gff))
      last = e.end3
    return introns
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
    
class gtfgene(exon):
  def __init__(self, lst, gff = False, addchr = False):
    exon.__init__(self, lst, gff, addchr)
    #self.pos = exon(lst)
    #if gff : self.id = self.attr['GeneID']
    #else : self.id = self.attr['gene_id']
    self.id = self.gid
    self.trans = []
    self.type = lst[1]
    #self.attr = self.pos.attr
  def add_trans(self, tr):
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
  def merge_trans(self):
    es = []
    merge = gtftrans(self.lst)
    for t in self.trans:
      es += t.exons
    es.sort()
    me = es[0]
    for i in range(1, len(es)):
      if me.stop >= es[i].start:
        me = union(me, es[i])[0]
      else: 
        merge.add_exon(me)
        me = es[i]
    merge.add_exon(me)
    merge.check()
    return merge
      
    
def load_gtf(fin, filt = [], gff = False, addchr = False):
  genes = {}
  trans = {}
  for l in fin:
    if l[0] == '#' : continue
    lst=l.strip().split('\t')
    if lst[2] == 'region' : continue
    if lst[2] == 'gene':
      if filt != [] and lst[1] not in filt: 
        continue
      g = gtfgene(lst, gff, addchr)
      genes[g.id] = g
      #print g.id
    elif lst[2] == 'transcript':
      t = gtftrans(lst, gff, addchr)
      #g.add_trans(t)
      if t.gid not in genes:
        g = gtfgene(lst, gff, addchr)
        genes[g.id] = g
        #continue ##
      genes[t.gid].add_trans(t)
      trans[t.id] = t
      #print t.id
    else:
      e = exon(lst, gff, addchr)
      if e.gid not in genes:
        g = gtfgene(lst, gff, addchr)
        genes[g.id] = g
        #continue
      #print e.attr['exon_number']
      if e.tid not in trans:
        t = gtftrans(lst, gff, addchr)
        genes[t.gid].add_trans(t)
        trans[t.id] = t
      trans[e.tid].add_exon(e)
  return genes, trans

def fetch_gtf(fin, gid = '', tid = '', gff = False, addchr = False):
  genes = {}
  trans = {}
  if gid == '' and tid == '' : return genes, trans
  ch = ''
  for l in fin:
    if l[0] == '#' : continue
    lst=l.strip().split('\t')
    if ch != '' and ch != lst[0] : break
    if gid != '' and lst[8].find(gid) < 0 : continue
    if tid != '' and lst[8].find(tid) < 0 : continue
    if lst[2] == 'region' : continue
    e = exon(lst, gff, addchr)
    if e.tid != tid  and e.gid != gid : continue
    if e.gid == '' : continue
    if lst[2] == 'gene':
      g = gtfgene(lst, gff, addchr)
      if g.id == gid : genes[g.id] = g
      #print g.id
    elif lst[2] == 'transcript':
      t = gtftrans(lst, gff, addchr)
      if t.id != tid  and t.gid != gid : continue 
      #g.add_trans(t)
      if t.gid not in genes:
        g = gtfgene(lst, gff, addchr)
        genes[g.id] = g
        #continue ##
      genes[t.gid].add_trans(t)
      trans[t.id] = t
      #print t.id
    else:
      if e.gid not in genes:
        g = gtfgene(lst, gff, addchr)
        genes[g.id] = g
        #continue
      #print e.attr['exon_number']
      if e.tid not in trans:
        t = gtftrans(lst, gff, addchr)
        genes[t.gid].add_trans(t)
        trans[t.id] = t
      trans[e.tid].add_exon(e)
    ch = g.chr
  return genes, trans

def gtfgene_iter(fin, filt = [], gff = False, addchr = False):
  genes = {}
  trans = {}
  gidlist = []
  chr = ""
  for l in fin:
    if l[0] == '#' : continue
    lst=l.strip().split('\t')
    if lst[2] == 'region' : continue
    if lst[0] != chr:
      for gid in gidlist:
        genes[gid].check()
        yield genes[gid]
      genes = {}
      trans = {}
      gidlist = []
      chr = lst[0]
    e = exon(lst, gff, addchr)
    if filt != [] and e.gid not in filt : continue
    if lst[2] == 'gene':
      #if filt != [] and lst[1] not in filt: 
        #continue
      g = gtfgene(lst, gff, addchr)
      genes[g.id] = g
      if g.id not in gidlist: gidlist.append(g.id)
    elif lst[2] == 'transcript':
      t = gtftrans(lst, gff, addchr)
      if t.gid not in genes:
        g = gtfgene(lst, gff, addchr)
        genes[g.id] = g
        gidlist.append(g.id)
      genes[t.gid].add_trans(t)
      trans[t.id] = t
    else:
      #e = exon(lst, gff)
      if e.gid == '' or e.tid == '' : continue
      if e.gid not in genes:
        g = gtfgene(lst, gff, addchr)
        genes[g.id] = g
        gidlist.append(g.id)
      if e.tid not in trans:
        t = gtftrans(lst, gff, addchr)
        genes[t.gid].add_trans(t)
        trans[t.id] = t
      trans[e.tid].add_exon(e)
  for gid in gidlist:
    genes[gid].check()
    yield genes[gid]
    
def gtftrans_iter(fin, filt = [], gff = False, addchr = False):
  genes = {}
  trans = {}
  gidlist = []
  chr = ""
  for l in fin:
    if l[0] == '#' : continue
    lst=l.strip().split('\t')
    if lst[2] == 'region' : continue
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
      g = gtfgene(lst, gff, addchr)
      genes[g.id] = g
      if g.id not in gidlist: gidlist.append(g.id)
    elif lst[2] == 'transcript':
      t = gtftrans(lst, gff, addchr)
      if t.gid not in genes:
        g = gtfgene(lst, gff, addchr)
        genes[g.id] = g
        gidlist.append(g.id)
      genes[t.gid].add_trans(t)
      trans[t.id] = t
    else:
      e = exon(lst, gff, addchr)
      if e.gid not in genes:
        g = gtfgene(lst, gff, addchr)
        genes[g.id] = g
        gidlist.append(g.id)
      if e.tid not in trans:
        t = gtftrans(lst, gff, addchr)
        genes[t.gid].add_trans(t)
        trans[t.id] = t
        #print e.attr['transcript_id'], t.id
      trans[e.tid].add_exon(e)
  for gid in gidlist:
    for t in genes[gid].trans:
      yield t
  #return genes, trans

def sub(a, b): # a - b
  if a.chr != b.chr : return [a]
  if a.start >= b.stop or a.stop <= b.start : return [a]
  out = []
  if a.start < b.start : out.append(a(stop = b.start))
  if b.stop < a.stop: out.append(a(start = b.stop))
  return out

def intersect(a, b): # a & b
  out = []
  if a.chr != b.chr : return out
  if a.start >= b.stop or a.stop <= b.start : return out
  out.append(a(start = max(a.start, b.start), stop = min(a.stop, b.stop)))
  return out

def union(a, b): # a U b
  if a.chr != b.chr or a.start > b.stop or a.stop < b.start :
    return [a, b]
  return [a(start = min(a.start, b.start), stop = max(a.stop, b.stop))]

