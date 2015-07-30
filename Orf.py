from zbio import Tools

codonSize = 3
cstart = ['ATG']
cstartlike = ['TTG', 'CTG', 'GTG', 'AAG', 'AGG', 'ACG', 'ATT', 'ATA', 'ATC'] #AAG & AGG may be removed
cstop = ['TGA', 'TAA', 'TAG']

senseframe = [1, 2, 3]
antiframe = [-1, -2, -3]
frame = [1, 2, 3, -1, -2, -3]

class Orf:
  def __init__(self, lst = [], frame = 0, stop = -1):
    if len(lst) != 0:
      #lst = l.strip().split('\t')
      self.frame = int(lst[0])
      self.starts = map(int, lst[1].split(','))
      self.altstarts = map(int, lst[2].split(','))
      self.stop = int(lst[3])
    else:
      self.frame = frame
      self.starts = []
      self.altstarts = []
      self.stop = stop
  def __str__(self):
    return "%d\t%s\t%s\t%d" % (self.frame, ','.join(map(str, self.starts)), ','.join(map(str, self.altstarts)), self.stop)
  def __repr__(self):
    return "ORF object: " + str(self)
  def has_start(self):
    return len(self.starts) + len(self.altstarts) > 0
  def has_stop(self):
    return self.stop >= 0
  def has_strictstart(self):
    return len(self.starts) > 0
  def is_complete(self):
    return self.has_start() and self.has_stop()
  def orfByRibo(self, frame, start, stop, srange = 4):
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

def allORF(seq, strand = '+') :
  seq = seq.upper().replace('U','T')
  if strand == '+' : fr = senseframe
  elif strand == '-' : 
    fr = antiframe
    antiseq = Tools.rc(seq)
  else: 
    fr = frame
    antiseq = Tools.rc(seq)
  
  length = len(seq)
  for f in fr:
    if f > 0:
      s = seq
      fa = f
    else: 
      s = antiseq
      fa = -f
    orf = Orf(frame = f)
    for i in range(fa-1, length, codonSize):
      try: codon = s[i:i+codonSize]
      except: break
      if codon in cstart:
        orf.starts.append(i)
      elif codon in cstartlike:
         orf.altstarts.append(i)
      elif codon in cstop:
        orf.stop = i + codonSize
        if orf.has_start():
          yield orf
          orf = Orf(frame = f)

def findORF(seq, strand = '+', altcstart = False) :
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
        orfstop = i + codonSize
        if orfstart >= 0:
          orflen = (i - orfstart) / 3
          if f > 0:
            yield orfstart, orfstop, f, orflen
          else:
            yield length - orfstop, length - orfstart, f, orflen
          orfstart = orfstop = -1
          
          