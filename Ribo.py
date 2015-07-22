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
  if m <= 0 : return f
  c = 0
  for i in range(l):
    if arr[i] == m: 
      f[i] = 1.0
      c += 1
  if c <= 1: return f
  for i in range(l):
    f[i] /= c
  return f
def firstFrames(arr, bin = 3):
  fs = []
  for i in range(0, len(arr), bin):
    fs.append(first(arr[i:i+bin])) # May out of index
  return fs
def frameMean(fs):
  l = len(fs[0])
  f = [0] * l
  c = 0
  for i in range(len(fs)):
    for j in range(l):
      f[j] += fs[i][j]
      #if sum(fs[i]) > 0 : c += 1
    c += sum(fs[i])
    #print fs[i],c
  if c > 0 : 
    for j in range(l): f[j] /= c
  return f
def frameTestN(arr, length, value, bin = 3, n = 1000, show = False): #expect no bias
  a = arr[:]
  c = 0
  for i in range(n):
    random.shuffle(a)
    fs = firstFrames(a[0:length], bin)
    f = frameMean(fs)
    if max(f) >= value: c += 1
    if show : print a, fs, f, c
    if i == 9 and c > 5 : return float(c) / 10
    if i == 99 and c > 15 : return float(c) / 100
  return float(c) / n
def frameTestE(arr, length, value, expect, bin = 3, n = 1000, show = False):
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
    for i in range(length/bin):
      ff = []
      for j in range(bin):
        ff.append(af[j][i])
      fs.append(first(ff))
    f = frameMean(fs)
    m = max(f)
    if f[expect] < m and m >= value : c += 1
    if show : print fs, f, c
    if i == 9 and c > 5 : return float(c) / 10
    if i == 99 and c >15 : return float(c) / 100 ###
  return float(c) / n
def bias(f):
  m = max(f)
  for i in range(len(f)):
    if f[i] == m : return i
  
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
    self.start = start #binned 
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
    return Region(self.ers, self.start, self.stop, self.n, self.score, self.p)
  def __str__(self):
    return "%d-%d n=%d score=%s p=%s" % (self.start, self.stop, self.n, str(self.score), str(self.p))
  def __repr__(self):
    return "Enriched Region object " + str(self)
  def mapback(self, start = -1, stop = -1, nhead = -1, bin = -1):
    if nhead < 0 : nhead = self.ers.nhead
    if bin < 0 : bin = self.ers.bin
    if start < 0 : start = self.start 
    if stop < 0 : stop = self.stop
    start = start * bin + nhead
    stop = stop * bin + nhead
    return start, stop
  def binomPval(self, l = -1):
    p = float(len(self))/self.ers.length
    if l > 0: p2 = float(len(self))/l
    else : p2 = p
    return Stat.binomTest(self.n, self.ers.total, p) / p2
  
  def getFrame(self, start = -1, stop = -1):
    if start < 0 : start = self.start 
    if stop < 0 : stop = self.stop
    f = frameMean(self.ers.frames[start:stop])
    return f
  def localPos(self, r, window = 20):
    l = len(self)
    rs = int(r * l - window / 2)
    re = rs + window
    if rs < 0 :
      re -= rs
      rs = 0
    if re >= l:
      rs -= re - l
      re = l
    rs += self.start
    re += self.start
    return (rs, re)
  
  def frameCheck(self, n = 1000, window = 20, local = [0.0,0.5,1.0]):
    f = self.getFrame(self.start, self.stop)
    fm = max(f)
    tstart, tstop = self.mapback()
    p = frameTestN(self.ers.cnts[tstart:tstop], tstop-tstart, fm, self.ers.bin, n = n)
    #print "Region frame:",
    if p < 0.05 : 
      b = bias(f)
      #print "Bias =", b, "p =", p, "f =", f
    #else: print "No significant Bias, p =", p, "f =", f
    regbias = (f, p)
    locbias = {}
    l = len(self)
    if l > window :
      for r in local:
        (rs, re) = self.localPos(r, window)
        (ts, te) = self.mapback(rs, re)
        rf = self.getFrame(rs, re)
        rfm = max(rf)
        rb = bias(rf)
        if p < 0.05 : 
          if round(rfm, 3) == round(rf[b], 3) : continue
          rpn = frameTestN(self.ers.cnts[ts:te], te-ts, rfm, self.ers.bin, n = n)
          if rpn < 0.05 :
            rpe = frameTestE(self.ers.cnts[tstart:tstop], te-ts, rfm, b, self.ers.bin, n = n)
          else : rpe = 1
          rp = max(rpn, rpe)
          if rp < 0.05 : tp = '1-1'
          else : tp = '1-0'
        else: 
          rp = frameTestN(self.ers.cnts[ts:te], te-ts, rfm, self.ers.bin, n = n)
          if rp < 0.05 : tp = '0-1'
          else : tp = '0-0'
        if rp < 0.05 : 
          #print "Local bias: r =", r, ", rBias =", rb, ", rp =", rp, ", rf =", rf
          locbias[r] = (rf, rp)
    return (regbias, locbias)

class enrichedRegions:
  def __init__(self, cnts, nhead = 12, ntail = 18, bin = 3, log = True):
    self.cnts = cnts
    self.nhead = nhead
    self.ntail = ntail
    self.bin = bin
    self.bins = []
    self.frames = []
    self.total = 0
    length = len(cnts) - nhead - ntail
    if length < bin * 2: 
      #print "Transcript too short!"
      return 
    for i in range(nhead, len(cnts)-ntail, bin):
      try: 
        if log: s = sum(map(intlog2, cnts[i:i+bin]))
        else: s = sum(cnts[i:i+bin])
      except: break
      self.total += s
      self.bins.append(s)
      self.frames.append(first(cnts[i:i+bin]))
    if self.total == 0 : 
      #print "No reads!"
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
    if self.total == 0 : raise StopIteration
    rarr = []
    rmax = Region(self, 0, 1, 0, 0)
    #print self.rStarts, self.rStops
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
        if score > 0 : r.p = r.binomPval()
        rarr[-1].append(r)
        if rmax < r : rmax = r
        lasti2 = bi2
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
