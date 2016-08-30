def add_chr(chr):
  if chr.isdigit() or chr in ('X','Y','M',): chr = 'chr' + chr
  elif chr == 'MT' : chr = 'chrM' 
  return chr
def rm_chr(chr):
  if chr[3:].isdigit() or chr in ('chrX','chrY'): chr = chr[3:]
  elif chr == 'chrM' : chr = 'MT'
  return chr
class exon:
  def __init__(self, lst, gff = False, addchr = False):
    self.chr, self.strand, self.start, self.stop = lst[0], lst[6], int(lst[3]) - 1, int(lst[4])
    if addchr and self.chr[0:3] != 'chr' :
      if self.chr.isdigit() or self.chr in ('X','Y','M',): self.chr = 'chr' + self.chr
      elif self.chr == 'MT' : self.chr = 'chrM' 
    self.type, self.score, self.attrstr = lst[2], lst[5], lst[8]
    #self.attr = attr(lst[8], gff)
    self.gff, self.addchr, self.frame = gff, addchr, lst[7]
    self.lst = lst
  def __repr__(self):
    return self.chr + ':'+str(self.start)+'-'+str(self.stop)+':'+self.strand
  def __str__(self, sep = '\t', new = False):
    if not new : return sep.join(self.lst)
    lst = [self.chr, self.lst[1], self.type, str(self.start+1), str(self.stop), str(self.score), self.strand, str(self.frame), self.attrstr]
    return sep.join(lst)
  def string(self, sep = '\t', new = False):
    return self.__str__(sep = sep, new = new)
  def short_str(self):
    return "%s:%d-%d:%s" % (self.chr, self.start, self.stop, self.strand)
  def __cmp__(self, other):
    return cmp(self.chr,other.chr) or cmp(self.start,other.start) or cmp(self.stop,other.stop)
  def __len__(self): 
    return self.stop - self.start
  def attr(self, key):
    s = self.attrstr
    if self.gff : 
      p1 = s.find(key)
      if p1 < 0 : return ''
      else : p1 += len(key)
      p2 = s.find(';', p1)
      if s[p1] == '=' : return s[p1:p2]
      elif s[p1] == ':' : return s[p1:p2].split(',')[0]
      else : return ''
    p1 = s.find(key + ' ')
    if p1 < 0 : return ''
    else : p1 += len(key) + 1
    p2 = s.find(';', p1)
    return eval(s[p1:p2])
  @property
  def length(self):
    return len(self)
  def cdna_length(self): #Bed3 and Bed6
    return len(self)
  def is_reverse(self): 
    return self.strand=='-'
  @property
  def anti_strand(self):
    if self.strand == '+': return '-'
    elif self.strand == '-': return '+'
    else : return '.'
  @property
  def gid(self):
    try : return self.gid_c
    except : 
      if self.gff : self.gid_c = self.attr('GeneID')
      else : self.gid_c = self.attr('gene_id')
    return self.gid_c
  @property
  def tid(self):
    try : return self.tid_c
    except : 
      self.tid_c = self.attr('transcript_id')
      return self.tid_c
  @property
  def symbol(self):
    return self.attr('gene_name')
  @property
  def id(self):
    return self.attr('exon_id')
  @property
  def end5(self): #5' end, all bed
    if self.strand != '-' : return self.start
    else : return self.stop
  @property
  def end3(self): #3' end, all bed
    if self.strand != '-' : return self.stop
    else : return self.start
  @property
  def genetype(self):
    a = self.attr('gene_biotype')
    if a != '' : return a
    a = self.attr('gene_type')
    if a != '' : return a
    if self.tid[0:3] == 'NM_' : return "protein_coding"
    return a
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
  def is_contain(self,p): #if i in exon
    return self.start<=p<=self.stop
  def cdna_pos(self,p): # exon only
    if self.is_contain(p): return abs(p-self.end5)
    else: return None
  def genome_pos(self, p): # exon only
    m = self.cdna_length()
    if p < 0 or p > m: return None
    if p == 0 : return self.end5
    if p == m: return self.end3
    if not self.is_reverse() : return self.start + p
    else : return self.stop - p
  def copy(self, other):
    self.chr, self.strand, self.start, self.stop = other.chr, other.strand, other.start, other.stop
    self.type, self.score, self.attrstr = other.type, other.score, other.attrstr
    #self.attr = other.attr
    self.gff, self.addchr, self.frame = other.gff, other.addchr, other.frame
    self.lst = other.lst
  @property
  def exons(self):
    return [self]
  def cds_start(self, cdna = False): return None
  def cds_stop(self, cdna = False): return None
    
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
  for att in l:
    att = att.strip()
    p = att.find(' ')
    a[att[0:p]] = eval(att[p+1:])
  return a
def attrstr(a) :  # Not tested!
  s = ''
  for k in a :
    s += k + ' '
    if type(a[k]) == str : s += '"' + a[k] + '"'
    else : s += str(a[k])
    s += '; '
  return s
def tup2com(t): #tuple to comma string
  if type(t)==str: return t
  return ','.join(map(str,t))+','

class gtftrans(exon):
  def __init__(self, lst, gff = False, addchr = False): 
    exon.__init__(self, lst, gff, addchr)
    self.id = self.tid
    #self.type = lst[1]
    self.exons = []
    self.cds = []
    self.utr = []
    self.start_codon = None
    self.stop_codon = None
    self.other = []
  def add_exon(self, e): 
    if e.type == 'exon': self.exons.append(e)
    elif e.type == 'CDS': self.cds.append(e)
    elif e.type == 'UTR': self.utr.append(e)
    elif e.type == 'start_codon': self.start_codon = e
    elif e.type == 'stop_codon': self.stop_codon = e
    else: self.other.append(e)
    #self.check()
  def __repr__(self):
    s = self.symbol + " transcript_id " + self.id + ', ' + str(len(self.exons)) + " exons, " + exon.__repr__(self)
    return s
  def __str__(self, sep = '\t', new = False):
    arr = self.exons + self.cds + self.utr + self.other
    if self.start_codon is not None : arr.append(self.start_codon)
    if self.stop_codon is not None : arr.append(self.stop_codon)
    arr.sort(reverse = self.is_reverse())
    s = ''
    if self.type == 'transcript' : s += exon.__str__(self, new = new) + '\n'
    for e in arr : 
      s += e.__str__(sep = sep, new = new) + '\n'
    return s.strip()
  def string(self, sep = '\t', new = False):
    return self.__str__(sep = sep, new = new)
  def allRecords(self):
    arr = self.exons + self.cds + self.utr + self.other
    if self.start_codon is not None : arr.append(self.start_codon)
    if self.stop_codon is not None : arr.append(self.stop_codon)
    arr.sort(reverse = self.is_reverse())
    if self.type == 'transcript' : arr[0:0] = [self]
    return arr
  def genePredStr(self, geneName = True, extended = False, slice = None):
    lst = []
    if geneName : lst += [self.gid]
    start, stop, thickstart, thickstop = self.start, self.stop, self.thick_start, self.thick_stop
    if slice is None :
      exons = self.exons
    else :
      if start < slice[0] : start = slice[0]
      if thickstart < slice[0] : thickstart = slice[0]
      if stop > slice[1] : stop = slice[1]
      if thickstop > slice[1] : thickstop = slice[1]
      exons = []
      for e in self.exons:
        if e.stop <= slice[0] : continue
        if e.start >= slice[1] : continue
        exons.append(e(start = max(e.start, slice[0]), stop = min(e.stop, slice[1])))
    lst += [self.id, self.chr, self.strand, start, stop, thickstart, thickstop, len(exons)]
    exonStarts = [x.start for x in exons]
    exonStops = [x.stop for x in exons]
    exonStarts.sort()
    exonStops.sort()
    lst += [tup2com(exonStarts), tup2com(exonStops)]
    return '\t'.join(map(str, lst))
  def check(self):
    for e in self.exons:
      if e.start < self.start: self.start = e.start
      if e.stop > self.stop: self.stop = e.stop
    self.exons.sort(reverse = self.is_reverse())
    self.cds.sort(reverse = self.is_reverse())
  #@property
  def cds_start(self, cdna = False): 
    if self.start_codon is None and len(self.cds) == 0 : return None
    if self.start_codon is None :
      self.cds.sort(reverse = self.is_reverse())
      gs = self.cds[0].end5
      if self.stop_codon is not None : 
        cds2 = self.cds_stop(cdna = True)
        cds1 = self.cdna_pos(gs)
        cs = cds2 - (cds2 - cds1) / 3 * 3
        if cdna : return cs
        else : return self.genome_pos(cs, 1)
      try : frame = int(self.cds[0].frame)
      except : frame = 0
      if frame > 0 : 
        cs = self.cdna_pos(gs) + frame
        if cdna : return cs
        else : return self.genome_pos(cs, 1)
      if cdna : return self.cdna_pos(gs)
      else : return gs
    sc = self.start_codon
    if len(sc) != 3 : 
      for e in self.exons:
        if sc.end5 == e.end5 : 
          cs = self.cdna_pos(sc.end3) - 3
          if cdna : return cs
          else : return self.genome_pos(cs, 1)
        if sc.end3 == e.end3 :
          break
      else : 
        print 'Start codon error : %s %s %s' % (self.gid, self.id, sc.short_str())
    gs = sc.end5
    if cdna : return self.cdna_pos(gs)
    else : return gs
  def cds_stop(self, cdna = False):
    if self.stop_codon is None and len(self.cds) == 0 : return None
    if self.stop_codon is None :
      self.cds.sort(reverse = self.is_reverse())
      gs = self.cds[-1].end3
      if self.start_codon is not None : # determine by paired start codon
        cds1 = self.cds_start(cdna = True)
        cds2 = self.cdna_pos(gs)
        cs = cds1 + (cds2 - cds1) / 3 * 3
        if cdna : return cs
        else : return self.genome_pos(cs, 0)
      try : frame = int(self.cds[-1].frame) # determine by frame
      except : frame = None
      if frame is not None : 
        cds1 = self.cdna_pos(self.cds[-1].end5) + frame
        cds2 = self.cdna_pos(gs)
        cs = cds1 + (cds2 - cds1) / 3 * 3
        if cdna : return cs
        else : return self.genome_pos(cs, 1)
      if cdna : return self.cdna_pos(gs) #cs
      else : return gs # self.genome_pos(cs, 0)
    sc = self.stop_codon # Determine by stop codon
    if len(sc) != 3 : 
      for e in self.exons:
        if sc.end3 == e.end3 : 
          cs = self.cdna_pos(sc.end5) + 3
          if cdna : return cs
          else : return self.genome_pos(cs, 0)
        if sc.end5 == e.end5 :
          break
      else : 
        print 'Stop codon error : %s %s %s' % (self.gid, self.id, sc.short_str())
    gs = sc.end3
    if cdna : return self.cdna_pos(gs)
    else : return gs
  @property
  def thick_start(self):
    if self.is_reverse(): s = self.cds_stop()
    else : s = self.cds_start()
    if s is None : return self.start
    else : return s
  @property
  def thick_stop(self):
    if self.is_reverse(): s = self.cds_start()
    else : s = self.cds_stop()
    if s is None : return self.start
    else : return s
  @property
  def itemRgb(self):
    return '0,0,0'
  @property
  def blockCount(self):
    return len(self.exons)
  @property
  def blockSizes(self):
    es = self.exons[:]
    if self.is_reverse() : es.sort()
    return ','.join([str(len(e)) for e in es])
  @property
  def blockStarts(self):
    es = self.exons[:]
    if self.is_reverse() : es.sort()
    return ','.join([str(e.start - self.start) for e in es])
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
  def copytrans(self, o):
    self.exons, self.cds, self.utr, self.start_codon, self.stop_codon, self.others = o.exons, o.cds, o.utr, o.start_codon, o.stop_codon, o.others
    
class gtfgene(exon):
  def __init__(self, lst, gff = False, addchr = False):
    exon.__init__(self, lst, gff, addchr)
    self.id = self.gid
    self.trans = []
    #self.type = lst[1]
  def add_trans(self, tr):
    self.trans.append(tr)
    #self.check()
  def __repr__(self):
    s = "gene_id " + self.id + ', ' + str(len(self.trans)) + " transcripts, " + exon.__repr__(self)
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
  def max_cds(self): # NOT FINISHED!
    c1, c2 = None, None
    for t in self.trans:
      gcds1 = t.cds_start(cdna = False) 
      if gcds1 is not None : pass
        #if c1 is None or gcds1 <  : c1 = gcds1
        #else : pass
    
def load_gtf(fin, filt = [], gff = False, addchr = False):
  genes = {}
  trans = {}
  dfilt = {}
  if len(filt) > 0: 
    for gid in filt : dfilt[gid] = True
  chr = ''
  for l in fin:
    if l[0] == '#' : continue
    lst=l.strip().split('\t')
    if lst[0] != chr :
      for gid in genes : dfilt[gid] = False
      chr = lst[0]
    for gid in dfilt:
      if not dfilt[gid] : continue
      if lst[8].find(gid) >= 0 : break
    else : continue
    e = exon(lst, gff, addchr)
    if len(filt) > 0 and e.gid not in dfilt: 
      continue
    if lst[2] == 'region' : continue
    if lst[2] == 'gene':
      g = gtfgene(lst, gff, addchr)
      if g.id in genes: g.trans = genes[g.id].trans
      genes[g.id] = g
    elif lst[2] == 'transcript':
      t = gtftrans(lst, gff, addchr)
      if t.gid not in genes:
        g = gtfgene(lst, gff, addchr)
        genes[g.id] = g
      if t.id in trans: trans[t.id].copy(t)
      else : 
        genes[t.gid].add_trans(t)
        trans[t.id] = t
    else:
      e = exon(lst, gff, addchr)
      if e.gid not in genes:
        g = gtfgene(lst, gff, addchr)
        genes[g.id] = g
      if e.tid not in trans:
        t = gtftrans(lst, gff, addchr)
        genes[t.gid].add_trans(t)
        trans[t.id] = t
      trans[e.tid].add_exon(e)
  for g in genes : genes[g].check()
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
      if g.id in genes: g.trans = genes[g.id].trans
      if g.id == gid : genes[g.id] = g
    elif lst[2] == 'transcript':
      t = gtftrans(lst, gff, addchr)
      if t.id != tid  and t.gid != gid : continue 
      if t.gid not in genes:
        g = gtfgene(lst, gff, addchr)
        genes[g.id] = g
      if t.id in trans: trans[t.id].copy(t)
      else : 
        genes[t.gid].add_trans(t)
        trans[t.id] = t
    else:
      if e.gid not in genes:
        g = gtfgene(lst, gff, addchr)
        genes[g.id] = g
      if e.tid not in trans:
        t = gtftrans(lst, gff, addchr)
        genes[t.gid].add_trans(t)
        trans[t.id] = t
      trans[e.tid].add_exon(e)
    ch = g.chr
  for g in genes : genes[g].check()
  return genes, trans

def gtfgene_iter(fin, filt = [], gff = False, addchr = False, chrs = None, verbose = False):
  genes, genes2, trans, gidlist = {}, {}, {}, [] # genes2: genes with headline, for better output
  dfilt = {}
  if len(filt) > 0: 
    for k in filt : dfilt[k] = 1
  chr = ""
  for l in fin:
    if l[0] == '#' : continue
    lst=l.strip().split('\t')
    if lst[2] == 'region' : continue
    if lst[0] != chr:
      if chrs is not None and lst[0] not in chrs : continue
      if verbose : print lst[0]
      for gid in gidlist:
        if gid in genes : 
          genes[gid].check()
          yield genes[gid]
      genes, genes2, trans, gidlist = {}, {}, {}, []
      chr = lst[0]
    e = exon(lst, gff, addchr)
    if len(filt) > 0 and e.gid not in dfilt: 
        continue
    if e.gid not in genes : # check gene output
      dgs = []
      for gid in genes2 : 
        if e.start > genes2[gid].stop : 
          genes2[gid].check()
          yield genes2[gid]
          dgs.append(gid)
      for gid in dgs:
        del genes2[gid]
        del genes[gid]
    if lst[2] == 'gene':
      g = gtfgene(lst, gff, addchr)
      if g.id in genes: g.trans = genes[g.id].trans # Already exists due to disorder
      else : gidlist.append(g.id)
      genes[g.id] = g
      genes2[g.id] = g
      #if g.id not in gidlist: gidlist.append(g.id)
    elif lst[2] == 'transcript':
      t = gtftrans(lst, gff, addchr)
      if t.gid not in genes:
        g = gtfgene(lst, gff, addchr)
        genes[g.id] = g
        gidlist.append(g.id)
      if t.id in trans: trans[t.id].copy(t) # Already exists due to disorder
      else : 
        genes[t.gid].add_trans(t)
        trans[t.id] = t
    else:
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
    if gid not in genes : continue
    genes[gid].check()
    yield genes[gid]


def gtftrans_iter(fin, filt = [], gff = False, addchr = False, chrs = None, verbose = False): # To be updated
  genes, trans, trans2, tidlist = {}, {}, {}, []
  dfilt = {}
  if len(filt) > 0: 
    for k in filt : dfilt[k] = 1
  chr = ""
  for l in fin:
    if l[0] == '#' : continue
    lst=l.strip().split('\t')
    if lst[2] == 'region' : continue
    if lst[0] != chr:
      if chrs is not None and lst[0] not in chrs : continue
      if verbose : print lst[0]
      for tid in tidlist:
        if tid in trans : 
          trans[tid].check()
          yield trans[tid]
      genes, trans, trans2, tidlist = {}, {}, {}, []
      chr = lst[0]
    if lst[2] == 'gene': continue ## not gene structure
    e = exon(lst, gff, addchr)
    if len(filt) > 0 and e.tid not in dfilt: 
        continue
    if e.tid not in trans : # check gene output
      dgs = []
      for tid in trans2 : 
        if e.start > trans2[tid].stop : 
          trans2[tid].check()
          yield trans2[tid]
          dgs.append(tid)
      for tid in dgs:
        del trans2[tid]
        del trans[tid]
    if lst[2] == 'transcript':
      t = gtftrans(lst, gff, addchr)
      if len(filt) > 0 and t.tid not in dfilt: 
        continue
      if t.id in trans: trans[t.id].copy(t)
      else : 
        #genes[t.gid].add_trans(t)
        trans[t.id] = t
        trans2[t.id] = t
        tidlist.append(t.id)
    else:
      #e = exon(lst, gff, addchr)
      if e.tid not in trans:
        t = gtftrans(lst, gff, addchr)
        #genes[t.gid].add_trans(t)
        trans[t.id] = t
        tidlist.append(t.id)
        #print e.attr['transcript_id'], t.id
      trans[e.tid].add_exon(e)
  for tid in tidlist:
    if tid not in trans : continue
    trans[tid].check()
    yield trans[tid]
    #for t in genes[gid].trans:
      #yield t
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

