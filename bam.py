
import pysam
import bed

class bamfile(pysam.Samfile):
  def __repr__(self):
    return 'pysam.Samfile '+self.filename
  def fetch_reads(self, chr, start, stop, multiple_iterators=False):
    if chr not in self.references : raise StopIteration
    rds = self.fetch(reference=chr, start=start, end=stop, multiple_iterators=multiple_iterators)
    for read in rds:
      r = bam(read, bamfile)
      yield r

class bam():#AlignedRead
  
  def __init__(self, read, bamfile):
    self.read = read
    self.ref = bamfile.references
    self.lens = bamfile.lengths
  #def __getattr__(self,attr):
    #a=getattr(self,attr)
    #return a
  
  @property
  def chr(self): 
    return self.ref[self.read.tid]
  @property
  def start(self):
    return self.read.pos
  @property
  def stop(self):
    return int(self.read.aend)
  @property
  def id(self): 
    return self.read.qname
  @property
  def score(self): 
    try:
      return self.read.get_tag("AS")
    except:
      return 0.0
  @property
  def strand(self): 
    if self.read.is_reverse: return '-'
    else: return '+'
  @property
  def cigar(self): 
    return self.read.cigar
  
  def __len__(self): #All bed
    return self.read.alen
  @property
  def length(self):
    return len(self)
  
  def __repr__(self):
    return "Bam AlignedSegment " + self.id
  def __str__(self):
    return "%s:%d-%d:%s" % (self.chr, self.start, self.stop, self.strand)
  
  def cdna_length(self): 
    return self.read.qlen
  
  def center(self): #Middle point, NEED REVISE!!
    return (self.start+self.stop)/2.0
  
  def is_reverse(self): 
    return self.read.is_reverse
  
  @property
  def end5(self): #5' end, all bed
    if self.read.is_reverse: return self.stop
    else : return self.start
  @property
  def end3(self): #3' end, all bed
    if self.read.is_reverse: return self.start
    else : return self.stop
    
  def cdna_pos(self, p): #NOT READY
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
    if self.is_reverse():
      p1 = self.cdna_length() - p
      if bias == 0: bias = 1
      else: bias = 0
    pos = self.start
    #l=range(len())
    for ctype, l in self.cigar:
      if ctype in [0,7,8]:
        if l - p1 >= bias:
          pos += p1
          return pos
        else:
          pos += l
          p1 -= l
      elif ctype in [1,4]:
        if l - p1 >= bias:
          return pos
        else:
          p1 -= l
      elif ctype in [2,3]:
        pos += l
    return None
  @property
  def introns(self):
    s = []
    p = 0
    pos = self.start
    i = 1
    for ctype, l in self.cigar:
      if ctype in [0,7,8]:
        pos += l
        p += l
      elif ctype in [1,4]:
        p += l
      elif ctype in [2]:
        pos += l
      elif ctype in [3]:
        b = bed.bed6([self.chr,pos,pos + l,self.id+"_intron"+str(i),p,self.strand])
        s.append(b)
        pos += l
        i += 1
    if self.is_reverse():
      return s[::-1]
    return s
  def is_compatible(self, trans, mis = 0):
    m = 0
    last = -1
    for i in range(self.cdna_length()):
      pos = self.genome_pos(i, bias = 1)
      i2 = trans.cdna_pos(pos)
      if i2 is None : m += 1
      else :
        if last >= 0 :
          if i2 != last + 1 :
            m += min(i, self.cdna_length() - i, abs(i2 - last + 1))
        last = i2
      if m > mis : return False
    #print self, self.cigar, m
    if m > mis : return False
    else : return True
  def is_compatible0(self, trans = None, introns = [], mis = 0): # Bases of the read not in the right place
    #if trans != None : exons = trans.exons
    if trans is not None : introns = trans.introns
    m = 0
    ris = self.introns
    if len(ris) == 0 : 
      for intr in introns:
        o = self.read.get_overlap(intr.start, intr.stop)
        if o > 0 : m += o
        if m > mis : return False
    else :
      for i in range(len(introns)):
        if not ris[0].is_upstream(introns[i]) : break
      c = True
      for ri in ris:
        if ri != introns[i] : 
          c = False
          o = self.read.get_overlap(introns[i].start, introns[i].stop)
          p = int(ri.score)
          m += max(o, min(p, self.cdna_length() - p))
          if m > mis : return False
        i += 1
    return True
    
def compatible_bam_iter(bamfile, trans, mis = 0, sense = True):
  if trans.chr not in bamfile.references : raise StopIteration
  rds = bamfile.fetch(reference=trans.chr, start=trans.start, end=trans.stop)#, multiple_iterators=False)
  #introns = trans.introns
  for r in rds:
    read = bam(r, bamfile)
    if sense and read.strand != trans.strand:
      #print read.id + " not sense"
      continue
    o = read.read.get_overlap(trans.start, trans.stop)
    if o < read.cdna_length() - mis: 
      continue
    if read.is_compatible(trans = trans, mis = mis) :
      yield read
    #else:
      #print read.id, b, read.cdna_length()
      