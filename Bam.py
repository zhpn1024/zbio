
import pysam

class Bamfile(pysam.Samfile):
  def __repr__(self):
    return 'pysam.Samfile '+self.filename

class Bam():#AlignedRead
  
  def __init__(self,read,bamfile):
    self.read=read
    self.ref=bamfile.references
    
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
    
def compatibleBamIter(bamfile, trans, mis = 0, sense = True):
  rds = bamfile.fetch(reference=trans.chr, start=trans.start, end=trans.stop, multiple_iterators=False)
  exons = trans.exons
  for r in rds:
    read = Bam(r, bamfile)
    if sense and read.strand != trans.strand:
      #print read.id + " not sense"
      continue
    b = 0
    for e in exons:
      b += read.read.get_overlap(e.start, e.stop)
      #print b
    if read.cdna_length() - b <= mis:
      yield read
    #else:
      #print read.id, b, read.cdna_length()
      