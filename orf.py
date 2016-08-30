from zbio import tools

codonSize = 3
cstart = ['ATG']
cstartlike = ['TTG', 'CTG', 'GTG', 'AAG', 'AGG', 'ACG', 'ATT', 'ATA', 'ATC'] #AAG & AGG may be removed
caltstart = cstartlike
cstop = ['TGA', 'TAA', 'TAG']

senseframe = [1, 2, 3]
antiframe = [-1, -2, -3]
frame = [1, 2, 3, -1, -2, -3]

class orf:
  def __init__(self, lst = [], frame = 0, stop = -1):
    if len(lst) != 0:
      #lst = l.strip().split('\t')
      self.frame = int(lst[0])
      if '' == lst[1] : 
        self.starts = []
      else : self.starts = map(int, lst[1].split(','))
      if lst[2] == '' : self.altstarts = []
      else : self.altstarts = map(int, lst[2].split(','))
      self.stop = int(lst[3])
    else:
      self.frame = frame
      self.starts = []
      self.altstarts = []
      self.stop = stop
  def __str__(self):
    return "%d\t%s\t%s\t%d" % (self.frame, ','.join(map(str, self.starts)), ','.join(map(str, self.altstarts)), self.stop)
  def __repr__(self):
    return "orf object: " + str(self)
  def has_start(self):
    return len(self.starts) + len(self.altstarts) > 0
  def has_stop(self):
    return self.stop >= 0
  def __len__(self):
    if not self.is_complete() : return 0
    #if self.stop < self.start() : print str(self)
    return self.stop - self.start()
  def __cmp__(self, other):
    return cmp(len(self), len(other)) or cmp(self.start(),other.start)
  def length(self):
    return len(self)
  def aa_len(self):
    return len(self) / codonSize
  #@property
  def start(self, alt = True):
    if alt : return min(self.starts + self.altstarts)
    else :
      if len(self.starts) > 0 : return min(self.starts)
      else : return None
  def has_strictstart(self):
    return len(self.starts) > 0
  def is_complete(self):
    return self.has_start() and self.has_stop()
  def orf_by_ribo(self, frame, start, stop, srange = 4): #look for TIS near ribo suggested ORF
    if self.frame != frame: return None
    if self.stop < start: return None
    if min(self.starts + self.altstarts) > stop : return None
    dmin = srange * codonSize
    sm = -1
    for s in self.starts :
      d = abs(s - start)
      if d < dm :
        dm = d
        sm = s
    if sm > 0 : return (s, self.stop)
    for s in self.altstarts :
      d = abs(s - start)
      if d < dm :
        dm = d
        sm = s
    if sm > 0 : return (s, self.stop)
    return None
  def filtByLen(self, minaalen, tail = -1):
    stop = self.stop - 3
    if stop < 0 : stop = tail
    if stop < 0 : return
    th = minaalen * 3
    rm = False
    for i, s in enumerate(self.starts):
      if stop - s < th : 
        rm = True
        break
    if rm : self.starts[i:] = []
    rm = False
    for i, s in enumerate(self.altstarts):
      if stop - s < th : 
        rm = True
        break
    if rm : self.altstarts[i:] = []

def allorf(seq, strand = '+', minaalen = 0) :
  seq = seq.upper().replace('U','T')
  if strand == '+' : fr = senseframe
  elif strand == '-' : 
    fr = antiframe
    antiseq = tools.rc(seq)
  else: 
    fr = frame
    antiseq = tools.rc(seq)
  
  length = len(seq)
  for f in fr:
    if f > 0:
      s = seq
      fa = f
    else: 
      s = antiseq
      fa = -f
    o = orf(frame = f)
    for i in range(fa-1, length, codonSize):
      try: codon = s[i:i+codonSize]
      except: break
      if codon in cstart:
        o.starts.append(i)
      elif codon in cstartlike:
         o.altstarts.append(i)
      elif codon in cstop:
        o.stop = i + codonSize
        o.filtByLen(minaalen = minaalen, tail = length)
        if o.has_start():
          yield o
        o = orf(frame = f)
    o.filtByLen(minaalen = minaalen, tail = length)
    if o.has_start() : yield o
def orflist(seq, strand = '+', sort = True, minaalen = 0):
  ol = []
  for o in allorf(seq, strand = strand, minaalen = minaalen):
    ol.append(o)
  if sort : ol.sort(reverse = True)
  return ol

class fixedorf():
  def __init__(self, start, stop):
    self.start = start
    self.stop = stop
  def __len__(self):
    return self.stop - self.start
  def length(self):
    return len(self)
  def __repr__(self):
    return "Fixed ORF object: %d - %d" % (self.start, self.stop)
  def frame(self):
    return self.start % 3
  def __cmp__(self, other):
    return cmp(self.start, other.start) or cmp(self.stop, other.stop)
  def orfstr(self, seq):
    return "%s|%d-%d" % (seq[self.start:self.start+3], self.start, self.stop)
def orf_by_pos(seq, pos): ### Unkown start codon, only to find stop codon
  for i in range(pos, len(seq), codonSize):
    try: codon = seq[i:i+codonSize]
    except: break
    if codon in cstop: 
      i += codonSize
      break
  o = fixedorf(start = pos, stop = i)
  return o

def orfs_by_pos0(seq, pos): 
  orfs = []
  for f in range(codonSize):
    for i in range(pos+f, length, codonSize):
      try: codon = s[i:i+codonSize]
      except: break
      if codon in cstop: 
        orfs.append(pos+f, i+codonSize)
        break
    if len(orfs) <= f : orfs.append(pos+f, i)
  return orfs
def findorf(seq, strand = '+', altcstart = False) :
  seq = seq.upper().replace('U','T')
  if strand == '+' : fr = senseframe
  elif strand == '-' : 
    fr = antiframe
    antiseq = Tools.rc(seq)
  else: 
    fr = frame
    antiseq = Tools.rc(seq)
  if altcstart: cs = cstart + cstartlike
  else: cs = cstart
  
  length = len(seq)
  for f in fr:
    if f > 0:
      s = seq
      fa = f
    else: 
      s = antiseq
      fa = -f
    orfstart = orfstop = -1
    for i in range(fa-1, length, codonSize):
      try: codon = s[i:i+codonSize]
      except: break
      if codon in cs:
        if orfstart < 0: orfstart = i
      elif codon in cstop:
        orfstop = i + codonSize ##
        if orfstart >= 0:
          orflen = (i - orfstart) / 3
          if f > 0:
            yield orfstart, orfstop, f, orflen
          else:
            yield length - orfstop, length - orfstart, f, orflen
          orfstart = orfstop = -1
          
def orfdict(orflist, alt = True):
  od = {}
  for o in orflist:
    for s in o.starts:
      od[s] = o.stop
    if alt:
      for s in o.altstarts:
        od[s] = o.stop
  return od
def nearest_start(s, od, flank = 1):
  if s in od: return s
  for i in range(1, flank + 1):
    if s + i in od : return s + i ###
    if s - i in od : return s - i
  return None

codonTable = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L', 'TCT':'S','TCC':'S','TCA':'S','TCG':'S', 'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*', 'TGT':'C','TGC':'C','TGA':'*','TGG':'W', 
               'CTT':'L','CTC':'L','CTA':'L','CTG':'L', 'CCT':'P','CCC':'P','CCA':'P','CCG':'P', 'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q', 'CGT':'R','CGC':'R','CGA':'R','CGG':'R', 
               'ATT':'I','ATC':'I','ATA':'I','ATG':'M', 'ACT':'T','ACC':'T','ACA':'T','ACG':'T', 'AAT':'N','AAC':'N','AAA':'K','AAG':'K', 'AGT':'S','AGC':'S','AGA':'R','AGG':'R', 
               'GTT':'V','GTC':'V','GTA':'V','GTG':'V', 'GCT':'A','GCC':'A','GCA':'A','GCG':'A', 'GAT':'D','GAC':'D','GAA':'E','GAG':'E', 'GGT':'G','GGC':'G','GGA':'G','GGG':'G', 
              }

def translate(seq):
  aa = ""
  for i in range(0, len(seq), codonSize):
    try : a = codonTable[seq[i:i+codonSize].upper().replace('U','T')]
    except : a = 'X'
    aa += a
  return aa

def is_start(seq, pos, alt = False, flank = 0):
  for i in range(pos - flank, pos + flank + 1):
    try : c = seq[i:i+codonSize]
    except : continue
    if c in cstart : return True
    elif alt and c in cstartlike : return True
  return False
def is_startlike(seq, pos, flank = 0):
  return is_start(seq, pos, alt = True, flank = flank)
