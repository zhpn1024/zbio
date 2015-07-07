import random
from zbio import Stat

def intlog2(r):
  i = 0
  e = 1
  r += 1
  while e < r:
    e *= 2
    i += 1
  return i

def first(arr):
  l = len(arr)
  m = max(arr)
  f = [0] * l
  c = 0
  for i in range(l):
    if arr[i] == m: 
      f[i] = 1.0
      c += 1
  if c <= 1: return f
  for i in range(l):
    f[i] /= c
  return f
def firstFrames(arr, bin):
  fs = []
  for i in range(0, len(arr), bin):
    fs.append(first(arr[i:i+bin]))
  return fs
def frameMean(fs):
  l = len(fs[0])
  f = [0] * l
  for i in range(len(fs)):
    for j in range(l):
      f[j] += fs[i][j]
  for j in range(l):
    f[j] /= len(fs)
  return f
def frameTestN(arr, start, stop, value, bin, n = 500): #expect no bias
  a = arr[0:len(arr)]
  c = 0
  for i in range(n):
    random.shuffle(a)
    fs = firstFrames(a[start:stop], bin)
    f = frameMean(fs)
    if value <= max(f): c += 1
  return float(c) / n
def frameTestE(arr, start, stop, value, bin, expect, n = 500):
  a = arr[0:len(arr)]
  c = 0
  af = []
  for j in range(bin):
    af.append([])
  for i in range(0, len(arr), bin):
    for j in range(bin):
      af[j].append(arr[i+j])
  for i in range(n):
    for j in range(bin):
      random.shuffle(af[j])
    fs = []
    for i in range(len(af[0])):
      ff = []
      for j in range(bin):
        ff.append(af[j][i])
      fs.append(first(ff))
    f = frameMean(fs)
    m = max(f)
    if f[expect] < m and value <= m : c += 1
  return float(c) / n
  
  
def orfRead(cnts, cds1, cds2, nhead = 12, ntail = 18):
  all1 = nhead
  all2 = len(cnts) - ntail
  if cds1 < all1 : cds1 = all1
  if cds2 > all2 : cds2 = all2
  rall = 0
  rcds = 0
  for i in range(all1, all2):
    #if cnts[i] != 0: print i, cnts[i]
    rall += cnts[i]
    if cds1 <= i < cds2 : rcds += cnts[i]
  return rall, rcds, all2-all1, cds2-cds1

class Region:
  def __init__(self, ers, start, stop, n, score = -1, p = 1):
    self.ers = ers # corrent enrichedRegions object
    self.start = start
    self.stop = stop
    self.n = n
    self.score = score
    self.p = p
  def __cmp__(self, other):
    #return cmp(self.score, other.score)
    return cmp(other.p, self.p) or cmp(self.score, other.score) or cmp(len(self), len(other))
  def overlap(self, other):
    return self.start < other.stop and self.stop > other.start
  def __len__(self): 
    return self.stop - self.start
  def copy(self):
    return Region(self.start, self.stop, self.n, self.score, self.p)
  def __str__(self):
    return "%d-%d %d %s %s" % (self.start, self.stop, self.n, str(self.score), str(self.p))
  def __repr__(self):
    return "Enriched Region object " + str(self)
  def mapback(self, nhead = 12, bin = 3):
    start = self.start * bin + nhead
    stop = self.stop * bin + nhead
    return start, stop
  def binomPval(self, l = -1):
    
    p = float(len(self))/self.ers.length
    if l > 0: p2 = float(len(self))/l
    else : p2 = p
    return Stat.binomTest(self.n, self.ers.total, p) / p2
  
  def getFrame(self, start, stop):
    #fs = firstFrames(self.ers.cnts[start, stop], self.bin)
    f = frameMean(self.ers.frames[start, stop])
    return f
  
  def frameCheck(self, n = 500, window = 20):
    f = self.getFrame(self.start, self.stop)
    fm = max(f)
    bias = None
    start, stop = self.mapback(self.ers.nhead, self.ers.bin)
    p = frameTestN(self.ers.cnts, start, stop, fm, self.crs.bin, n = n)
    print "Region frame:", f
    if p < 0.05 : 
      for i in range(len(f)):
        if f[i] == fm: bias = i
      print "Bias:", i, "p =", p
    if 
    
    #l = len(fs)
    fall = [0,0,0]
    for i in range(self.start, self.stop):
      for j in range(3):
        fall[j] += self.ers.fs[i][j]
    for j in range(3):
      fall[j] /= len(self)
    print "All:", fall
    if len(self) <= 20: return 
    fall = [0,0,0]
    for i in range(self.start, self.start+20):
      for j in range(3):
        fall[j] += self.ers.fs[i][j]
    for j in range(3):
      fall[j] /= 20
    print "Head 20:", fall
    fall = [0,0,0]
    for i in range((self.start+self.stop)/2-10,(self.start+self.stop)/2+10):
      for j in range(3):
        fall[j] += self.ers.fs[i][j]
    for j in range(3):
      fall[j] /= 20
    print "Middle 20:", fall
    fall = [0,0,0]
    for i in range(self.stop-20, self.stop):
      for j in range(3):
        fall[j] += self.ers.fs[i][j]
    for j in range(3):
      fall[j] /= 20
    print "Tail 20:", fall
      
class enrichedRegions:
  def __init__(cnts, nhead = 12, ntail = 18, bin = 3, log = True):
    self.cnts = cnts
    self.nhead = nhead
    self.ntail = ntail
    self.bin = bin
    length = len(cnts) - nhead - ntail
    if length < bin * 2: 
      print "Cdna too short!"
      return 
    self.bins = []
    self.frames = []
    self.total = 0
    for i in range(nhead, len(cnts)-ntail, bin):
      try: 
        if log: s = sum(map(intlog2, cnts[i:i+bin]))
        else: s = sum(cnts[i:i+bin])
      except: break
      self.total += s
      self.bins.append(s)
      self.frames.append(first(cnts[i:i+bin]))
    if self.total == 0 : 
      print "No reads!"
      return 
    self.length = len(self.bins)
    self.mean = float(self.total) / self.length
    self.rStarts = [0]
    self.rStops = []
    for bi in range(1, self.length):
      if self.bins[bi-1] <= self.mean and self.bins[bi] >= self.mean:
        self.rStarts.append(bi)
      if self.bins[bi-1] >= self.mean and self.bins[bi] <= self.mean:
        self.rStops.append(bi)
    self.rStops.append(self.length)

  def __len__(self):
    return self.length
  def mapback(self, pos):
    if pos > self.length: return None
    return pos * self.bin + self.nhead
  
  def findRegion(self):
    rarr = []
    rmax = Region(self, 0, 1, 0, 0)
    for bi1 in self.rStarts:
      s = 0
      rarr.append([])
      lasti2 = bi1
      for bi2 in self.rStops:
        if bi2 <= bi1: continue
        for bi in range(lasti2, bi2):
          s += self.bins[bi]
        d = bi2 - bi1
      score = float(s) / self.total - float(d) / self.length
      r = Region(self, bi1, bi2, s, score)
      if score > 0 : r.p = r.binomPval(self.length, self.total)
      #print (s, total), (bi1,bi2, len(bs)), score, r.p
      rarr[-1].append(r)
      if rmax < r : rmax = r
      lasti2 = bi2
  #smax = rmax.score
  #rout = []
    while rmax.p < 0.05:# rmax.score > smax * 0.5: # and pvalue < XXX
    #rout.append(rmax.copy())
      yield rmax.copy() ####
      for rarr2 in rarr:
        for r in rarr2:
          if rmax.overlap(r): 
            r.score = -1
            r.p = 1
      rmax = Region(self, 0, 1, 0, 0)
      for rarr2 in rarr:
        for r in rarr2:
          if rmax < r : rmax = r

    
def enrichedRegion1(cnts, nhead = 12, ntail = 18, bin = 3, log = True):
  length = len(cnts) - nhead - ntail
  if length < bin * 2: 
    print "Cdna too short!"
    return None
  bs = []
  fs = []
  total = 0
  for i in range(nhead, len(cnts)-ntail, bin):
    try: 
      if log: s = sum(map(intlog2, cnts[i:i+bin]))
      else: s = sum(cnts[i:i+bin])
    except: break
    #a = float(s) / bin
    #if log: s = intlog2(s) ###
    total += s
    bs.append(s)
    fs.append(first(cnts[i:i+bin]))
  if total == 0 : 
    print "No reads!"
    return None
  mean = float(total) / len(bs)
  r1s = [0]
  r2s = []
  for bi in range(1,len(bs)):
    if bs[bi-1] <= mean and bs[bi] >= mean:
      r1s.append(bi)
    if bs[bi-1] >= mean and bs[bi] <= mean:
      r2s.append(bi)
  r2s.append(len(bs))
  print "Number of split sites: ", len(r1s), len(r2s)
  #curve = {}
  rarr = []
  rmax = Region(0, 1, 0, 0)
  for bi1 in r1s:
    s = 0
    rarr.append([])
    lasti2 = bi1
    for bi2 in r2s:
      if bi2 <= bi1: continue
      for bi in range(lasti2, bi2):
        s += bs[bi]
      d = bi2 - bi1
      score = float(s) / total - float(d) / len(bs)
      r = Region(bi1, bi2, s, score)
      r.p = r.binomPval(len(bs), total)
      #print (s, total), (bi1,bi2, len(bs)), score, r.p
      rarr[-1].append(r)
      if rmax < r : rmax = r
      elif rmax == r and len(rmax) < len(r): rmax = r
      #if rmax.p > r.p : rmax = r
      #elif rmax.p == r.p and len(rmax) < len(r): rmax = r
      #if d not in curve or curve[d] < s : curve[d] = s
      lasti2 = bi2
  smax = rmax.score
  rout = []
  while rmax.p < 0.05:# rmax.score > smax * 0.5: # and pvalue < XXX
    rout.append(rmax.copy())
    for rarr2 in rarr:
      for r in rarr2:
        if rmax.overlap(r): 
          r.score = -1
          r.p = 1
    rmax = Region(0, 1, 0, 0)
    for rarr2 in rarr:
      for r in rarr2:
        if rmax < r : rmax = r
        elif rmax == r and len(rmax) < len(r): rmax = r
  for r in rout:
    print r
    mb = r.mapback()
    print "DNA pos:", mb, "Genome pos:", (t.genome_pos(mb[0]), t.genome_pos(mb[1]))
    r.frame(fs)
  