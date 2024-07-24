'''
Genomic mutation processing
Copyright (c) 2017 Peng Zhang <zhpn1024@163.com>
'''
from . import fa, interval

try: cmp
except:
  def cmp(a, b):
    return (a > b) - (a < b)

def normvar(ref, alt):
  l1, l2 = len(ref), len(alt)
  if l1 == 1 or l2 == 1:
    return ref, alt, 0
  if l1 > l2:
    short, long = alt, ref
    refshort = False
  else:
    short, long = ref, alt
    refshort = True

  nt = 0 # trim length
  for i in range(1, len(short)):
    s = short[i:]
    if long.endswith(s):
      nt = len(short) - i
      break
  if nt > 0:
    ls = len(short) - nt
    short = short[0:ls]
    ls = len(long) - nt
    long = long[0:ls]

  if refshort: return short, long, nt 
  else: return long, short, nt


class Mut:
  def __init__(self, chr, pos, mutseq, reflen=1):
    self.chr = chr
    self.pos = pos
    self.reflen = reflen
    self.mutseq = mutseq
    self.indel = cmp(len(mutseq), reflen) #self.indel()
    if self.indel == 0: self.type = 'SNP'
    elif self.indel > 0: self.type = 'INS'
    else: self.type = 'DEL'

  #def indel(self):
    #return cmp(len(self.mutseq), self.reflen)

  def indel_len(self):
    return len(self.mutseq) - self.reflen

  def __str__(self):
    return '{}:{}:{}:{}:{}'.format(self.type, self.chr, self.pos, self.mutseq, self.reflen)

  def __repr__(self):
    return 'Mut object: ' + str(self)

  def __cmp__(self, other):
    return cmp(self.chr, other.chr) or cmp(self.pos, other.pos) or cmp(other.indel, self.indel)

  def __lt__(self, other):
    c = self.__cmp__(other)
    return c < 0

  def __gt__(self, other):
    c = self.__cmp__(other)
    return c > 0

  @property
  def end(self):
    return self.pos + self.reflen

def update(d, k, v, replace=True):
  if k in d:
    if replace: d[k] = v
  else:
    d[k] = v

class MutGenome(fa.Fa):
  def __init__(self, fapath, verbose = False):
    fa.Fa.__init__(self, fapath, verbose)
    self.reset_mut()

  def reset_mut(self):
    self.mutpos = {0:{}, -1:{}, 1:{}} # SNP, DEL & INS
    self.mutregion = {0:{}, -1:{}, 1:{}} # {0:interval.Interval(), -1:interval.Interval(), 1:interval.Interval()}
    self.allmut = {}

  def add_mut(self, mut, check=True, replace=False):
    chr = self.get_chrname(mut.chr)
    if chr is None :
      print("chr id {} not found in fasta file!".format(mut.chr))
      return -1
    if chr != mut.chr: mut.chr = chr
    s = str(mut)
    self.allmut[s] = mut
    indel = mut.indel
    update(self.mutpos[indel], chr, {}, replace=False)
    update(self.mutregion[indel], chr, interval.Interval(adj_merge=False), replace=False)
    if check:
      mi = interval.Interval(mut.pos, mut.end)
      for ind in self.mutregion:
        if chr in self.mutregion[ind]:
          mii = mi.intersect(self.mutregion[ind][chr])
          if len(mii) > 0:
            print('Mutation contradict with current mutation: {} {} {} {} {} {}'.format(mut, ind, chr, mii.start, mii.stop, self.mutpos[ind][chr][mii.start]))
            if not replace: return 1
    self.mutregion[indel][chr].add_itv([mut.pos, mut.end], check)
    update(self.mutpos[indel][chr], mut.pos, s)
    for p in range(mut.pos+1, mut.end):
      update(self.mutpos[indel][chr], p, s)
    return 0

  def check(self):
    for indel in self.mutregion:
      for chr in self.mutregion[indel]:
        self.mutregion[indel][chr].check()

  def get_mut(self, chr, start = 0, stop = -1, length = -1, flank = 0):
    chr, start, stop, length = self.get_pos(chr, start, stop, length)
    r = interval.Interval(start-flank, stop+flank)
    ml = []
    for indel in self.mutregion:
      if chr in self.mutregion[indel]:
        intersect = r.intersect(self.mutregion[indel][chr])
        for itv in intersect.lst:
          ml.append(self.allmut[self.mutpos[indel][chr][itv[0]]])
    ml.sort()
    return ml


  def fetch(self, chr, start = 0, stop = -1, length = -1, mutsites = None, offset = 0):
    chr, start, stop, length = self.get_pos(chr, start, stop, length)
    if chr is None :
      return None
    if length <= 0 : return ''
    sq = fa.Fa.fetch(self, chr, start, stop, length)
    r = interval.Interval(start, stop)
    md = {}
    for indel in self.mutregion:
      if chr in self.mutregion[indel]: 
        intersect = r.intersect(self.mutregion[indel][chr])
        for itv in intersect.lst:
          md[self.mutpos[indel][chr][itv[0]]] = 1
    ml = []
    for s in md:
      mut = self.allmut[s]
      #data = (mut.pos - start, -mut.indel(), mut.museq, mut.reflen)
      ml.append(mut)
    ml.sort(reverse = True)
    for mut in ml:
      p1, p2, ms = mut.pos - start, mut.end - start, mut.mutseq
      if p2 > length:
        p2, ms = length, ms[:length-p1] # not checked for p1<0 & p2>length at the same time 
      if p1 < 0:
        p1, ms = 0, ms[-p2:] # try to keep the same number of bases.
        #sq = mut.mutseq[-p2:] + sq[p2:] 
      #if mutpos is not None:
        #mutpos.append([0, mut.mutseq[-p2:], p2])
      #else:
      sq = sq[0:p1] + ms + sq[p2:]
      if mutsites is not None:
        mutsites.append([p1+offset, ms, p2-p1])
    return sq

  def intervalSeq(self, chr, regions, strand='+', mutsites = None):
    chr = self.get_chrname(chr)
    if chr is None :
      print("fa id {} not found in fasta file!".format(chr))
      return None
    s = []
    offset = 0
    for itv in regions:
      s.append(self.fetch(chr, itv[0], itv[1], mutsites=mutsites, offset=offset))
      offset += itv[1] - itv[0]
    sq = ''.join(s)
    if strand == '-': 
      sq = fa.rc(sq)
      if mutsites is not None:
        for mp in mutsites:
          mp[0] = offset - mp[0] - mp[2]
          mp[1] = fa.rc(mp[1])
    return sq

  def is_mutated(self, chr, start, stop = -1, length = -1):
    chr, start, stop, length = self.get_pos(chr, start, stop, length)
    if chr is None: return False
    r = interval.Interval(start, stop)
    for indel in self.mutregion:
      if chr in self.mutregion[indel]:
        intersect = r.intersect(self.mutregion[indel][chr])
        if len(intersect) > 0: return True
    return False

  def is_mutated_interval(self, chr, regions):
    chr = self.get_chrname(chr)
    if chr is None :
      print("chr id {} not found in fasta file!".format(chr))
      return False
    for itv in regions:
      if self.is_mutated(chr, itv[0], itv[1]): return True
    return False

