import random, numpy, math
from zbio import stat, bam, gtf, bed, interval, exp, orf, tools
from scipy.stats import ranksums, mannwhitneyu

codonSize = 3
maxNH = 5
secondary = False
minMapQ = 1 
minRlen = 6
minTransLen = 50
defOffset = 12
nhead, ntail = defOffset, 30 - defOffset
bin = 3
#offset = {}
#for i in range(26,35):
  #offset[i] = 12
def offset(read = 30, offdict = None):
  if offdict is None : return defOffset
  if type(read) == int : length = read # compatible with old versions!
  else : length = read.cdna_length()
  if 'm0' not in offdict or not read.is_m0() : od = offdict
  else : od = offdict['m0']
  if length in od : return od[length]
  else : return None

def bin_counts(arr, bin = bin):
  l = len(arr)
  lb = l / bin
  barr = [0] * lb
  for i in range(lb):
    p = i * bin
    barr[i] = sum(arr[p:(p+bin)])
  return barr
def rstest(x, y): # rank sum test p value of x > y
  st1, p1 = ranksums(x, y)
  try : st, p = mannwhitneyu(x, y)
  except : return 0.5
  if st1 > 0 : return p
  else: return 1 - p ###
class ribo: #ribo seq profile in transcript
  def __init__(self, trans, ribobam, offset = offset, offdict = None, maxNH = maxNH, minMapQ = minMapQ, secondary = secondary, compatible = True, mis = 2):
    self.length = trans.cdna_length()
    self.cnts = [0] * self.length
    self.nhead = nhead
    self.ntail = ntail
    self.bin = bin
    self.total = 0
    self.trans = trans
    for r in ribobam.fetch_reads(trans.chr, trans.start, trans.stop):
      if r.strand != trans.strand : continue ## Must be the same strand?
      try: 
        if r.read.get_tag('NH') > maxNH : continue
      except: pass
      if r.read.mapping_quality < minMapQ : continue
      if not secondary and r.is_secondary : continue
      if compatible : c = r.is_compatible(trans, mis = mis)
      else: c = r.is_inside(trans, mis = mis)
      if not c : continue ## 
      l = r.cdna_length()
      #if l not in offset: continue
      off = offset(r, offdict)
      if off is None: continue
      i = trans.cdna_pos(r.genome_pos(off)) ## Default bias
      if i is None : continue
      #if i not in self.cnts: self.cnts[i] = 0
      self.cnts[i] += 1
      self.total += 1
    #if total == 0 : continue
  def abdscore(self, norm = 1000):
    return math.log(float(self.total) / (self.length - 30) * norm)
  def enrich_test(self, start, stop): # Do not use
    inarr = self.cnts[start:stop]
    outarr = self.cnts[self.nhead:start] + self.cnts[stop:self.length-self.ntail]
    if len(outarr) <= 0: return None
    inbc = bin_counts(inarr)
    outbc = bin_counts(outarr)
    p = rstest(inbc, outbc)
    return p
  def frame_test_region(self, region, frame): # Frame test in a region
    inarr, outarr = [], []
    for i in region.num_iter():
      i = int(i)
      #print i
      if i % codonSize == frame : inarr.append(self.cnts[i])
      else : outarr.append(self.cnts[i])
    if len(inarr) <= 0 or len(outarr) <= 0 : return None
    p = rstest(inarr, outarr)
    return p
  def enrich_test_region(self, r1, r2, frame): # Enrich test between regions
    inarr, outarr = [], []
    for i in r1.num_iter():
      i = int(i)
      if i % codonSize == frame : inarr.append(self.cnts[i])
    for i in r2.num_iter():
      i = int(i)
      if i % codonSize == frame : outarr.append(self.cnts[i])
    if len(inarr) <= 0 or len(outarr) <= 0 : return None
    p = rstest(inarr, outarr)
    return p
  def frame_test(self, start, stop): # old frame test, to be replaced
    l = (stop - start) / bin
    inarr = [0] * l
    outarr = [0] * 2 * l
    for i in range(l):
      p = start + i * bin
      inarr[i] = self.cnts[p]
      outarr[i*2] = self.cnts[p + 1]
      outarr[i*2 + 1] = self.cnts[p + 2]
    if len(inarr) <= 0 or len(outarr) <= 0 : return 1
    p = rstest(inarr, outarr)
    return p
  def multi_orf_test(self, orflist): # Main function of multiple ORF detection for riboseq data
    tid = self.trans.id
    blank = []
    for i in range(codonSize): ## Blank regions in different frames
      b = interval.interval(start=self.nhead, stop=self.length-self.ntail, id=tid)
      blank.append(b)
    eps, fps = [1] * len(orflist), [1] * len(orflist)
    #indeps = [None] * len(orflist)
    gid = [None] * len(orflist) # Group ID
    grp = []
    grplst = [] # Group list
    # Group ORFs
    for i, o in enumerate(orflist): # o is orf.fixedORF object ##(start, stop)
      found = False
      for j, og in enumerate(grp):
        if o.stop == og.stop : # If belong to the group of og
          if og.start > o.start : grp[j] = o # Represented by the longest ORF
          gid[i] = j
          grplst[j].append(o)
          found = True
          break
      if not found : 
        grp.append(o)
        grplst.append([o])
        gid[i] = len(grp) - 1
    # Calculate ORF representative regions 
    for j, og in enumerate(grp):
      blank[og.frame()].sub_itv([og.start, og.stop]) # not blank in given frame
      grplst[j].sort()
      for k, o in enumerate(grplst[j]):
        if k == 0 : o.prev = None # To compare with upstream ORF in the same group
        else : o.prev = grplst[j][k-1]
        try : o.next = grplst[j][k+1] # Downstream ORF
        except : o.next = None
        #if o.next is None : stop = o.stop
        #else : stop = o.next.start  ## Inside group regions
        o.region = interval.interval(start = o.start, stop = o.stop)
        o.indr = indep_region(o, blank)
    blankall = blank[0]
    for i in range(1, codonSize):
      blankall = blankall.intersect(blank[i])
    blank.append(blankall)
    for i, o in enumerate(orflist):
      #print self.trans.id, o, o.region, blank
      eps[i], fps[i] = self.efpvalues(o, blank) #calculate enrichment and frame p-values
    return eps, fps
    
  def efpvalues(self, orf, blank): # For multi_orf_test
    if orf.indr.rlen() < minRlen : r1 = orf.region
    else : r1 = orf.indr
    fp = self.frame_test_region(r1, orf.frame())
    r2 = blank[codonSize]
    if r2.rlen() < minRlen or orf.indr.rlen() < minRlen : r2 = blank[orf.frame()]
    ep = self.enrich_test_region(r1, r2, orf.frame())
    ''' More complicated tests
    if orf.prev is None : 
      r2 = blank[codonSize]
      if r2.rlen() < minRlen or orf.indr.rlen() < minRlen : r2 = blank[orf.frame()]
      ep = self.enrich_test_region(r1, r2, orf.frame())
    else : 
      ep = None
      op = orf.prev
      while op is not None : # Should be larger than all upstream ORF regions
        if op.indr.rlen() < minRlen : r2 = op.region
        else : r2 = op.indr
        ep1 = self.enrich_test_region(r1, r2, orf.frame())
        if ep1 > ep : ep = ep1
        op = op.prev
    '''
    return ep, fp

  def orf_test(self, orflist): # to be replaced by multi_orf_test
    tid = self.trans.id
    blanklist = [bed.bed3(tid, self.nhead, self.ntail)]
    indeplist = []
    result = []
    for o in orflist:
      b = bed.bed3(tid, o.start, o.stop)
      b.frame = start % 3
      bl1 = [] #new blank list
      ib = [] # overlap of b and blank
      for blk in blanklist:
        bl1 += blk - b #update blank list
        ib += b.intersect(blk)
      blanklist = bl1
      pes = []
      pes.append(self.enrich_testl(ib, blanklist)) ## to be added
      pf = self.frame_test1(ib)
      pes = [pe0]
      for i in range(len(indeplist)):
        idp1 = [] #new independent region for each orf
        ovl = [] # overlap of b and independent region of orf[i]
        for bidp in indeplist[i]:
          idp1 += bidp - b
          ovl += b.overlap(bidp)
        if len(ovl) > 0: pes.append(self.enrich_testl(ovl, idp1))
        indeplist[i] = idp1  
      pe = stat.fisher_method(pes)
      indeplist.append(ib)
      result.append((pe, pf))
    return result
  
  def tis_test(self, start, r, p):
    zt = stat.ztnb(r, p)
    p = zt.pvalue(self.cnts[start])
    return p
  def is_summit(self, p, flank = 3):
    if p < nhead or p > self.length - ntail: return False
    if self.cnts[p] == 0: return False
    for i in range(1, flank + 1):
      if self.cnts[p] < self.cnts[p + i]: return False
      if self.cnts[p] < self.cnts[p - i]: return False
    return True
  def cnts_dict(self):
    cd = exp.readdict() #{}
    for i in range(nhead, len(self.cnts) - ntail):
      if self.cnts[i] > 0 : cd[i] = self.cnts[i]
    return cd
  def top_summits(self, n = 10, minratio = 0, flank = 3):
    slist = []
    for i in range(len(self.cnts)):
      if self.cnts[i] == 0 : continue
      if not self.is_summit(i, flank = flank): continue
      slist.append((i, self.cnts[i]))
    l = len(slist)
    if l <= 1 : return slist
    slist.sort(key = lambda x: x[1], reverse = True)
    minc = slist[0][1] * minratio ## min reads cutoff, not used
    if minc <= 0 : return slist[0:n]
    for i in range(l):
      if slist[i][1] < minc : return slist[0:min(i,n)]
    return slist[0:n]
  def top_summits_iter(self, minratio = 0, flank = 3): #generate all possible sites, do not specify total number
    slist = []
    for i in range(len(self.cnts)):
      if self.cnts[i] == 0 : continue
      if not self.is_summit(i, flank = flank): continue
      slist.append((i, self.cnts[i]))
    l = len(slist)
    if l < 1 : return 
    slist.sort(key = lambda x: x[1], reverse = True)
    minc = slist[0][1] * minratio ## min reads cutoff, not used
    #if minc <= 0 : return slist[0:n]
    for i in range(l):
      if slist[i][1] >= minc : yield slist[i]
    #return slist[0:n]

#### For multiple ORF detection
def indep_region(orf, blank):
  of = orf.frame
  r = orf.region
  for i in range(codonSize):
    if i != of : r = r.intersect(blank[i])
  return r

#### For tis
def pidx(value, lst, parts):
  l = len(lst) - 1
  for i in range(len(parts)):
    if lst[int(l * parts[i])] >= value : break
  return i
def pidx_uplim(i, lst, parts):
  l = len(lst) - 1
  return lst[int(l * parts[i])]

def estimate_tis_bg(gtfpath, bampath, genome, parts = [0.25, 0.5, 0.75], offset = offset, offdict = None, whole = True, maxqt = 0.95, skip_tis = True, alt_tis = True, tis_flank = 1, addchr = False):
  parts.sort()
  if parts[-1] < 1 : parts.append(1)
  bamfile = bam.bamfile(bampath, "rb")
  gtffile = open(gtfpath, 'r')
  fulldata, genes, sl = {}, [], []
  data = [exp.readdict() for i in range(len(parts))]
  #for i in range(len(parts)):
    #data.append({})
  lastchr = ''
  for g in gtf.gtfgene_iter(gtffile, addchr = addchr):
    if g.chr not in genome: continue
    if g.chr != lastchr: 
      print g.chr
      lastchr = g.chr
    merge = g.merge_trans()
    ml = merge.cdna_length()
    if ml < minTransLen : continue ##
    mtis = ribo(merge, bamfile, offset = offset, offdict = offdict, compatible = False)
    if mtis.total == 0 : continue
    mcds1, mcds2 = [], []
    for t in g.trans:
      cds1, cds2 = merge.cdna_pos(t.cds_start(cdna = False)), merge.cdna_pos(t.cds_stop(cdna = False))
      if cds1 is not None : mcds1.append(cds1)
      if cds2 is not None : mcds2.append(cds2)
    
    score = mtis.abdscore()
    sl.append(score)
    fulldata[g.id] = mtis.cnts_dict(), score, mcds1, mcds2
    genes.append(merge)
    #print g.id, fulldata[g.id]
  sl.sort()
  print 'Group data...'
  for m in genes:
    cnts, score, mcds1, mcds2 = fulldata[m.gid]
    #print m.gid, cnts
    ip = pidx(score, sl, parts)
    if whole : start = nhead ##
    elif len(mcds2) > 0: start = max(mcds2) ##The last stop codon
    else : continue ##
    ml = m.cdna_length()
    msq = tools.gtf2seq(genome, m)
    for i in range(start, ml - ntail):
      if i not in cnts: continue # ignore 0s
      if skip_tis and orf.is_start(msq, i, alt = alt_tis, flank = tis_flank) : continue # i in mcds1: continue
      #if cnts[i] not in data[ip]:
        #data[ip][cnts[i]] = 0
      #data[ip][cnts[i]] += 1
      data[ip].record(cnts[i])
      #print m.gid, ip, i, cnts[i]
      
  print 'Estimate ZTNB parameters...'
  paras = []
  l = len(sl) - 1
  ip = 0
  for ip, da in enumerate(data):
    maxcnt = da.quantile(maxqt) + 5 ##
    d = {}
    for i in range(maxcnt): ## max range 
      if i in da: d[i] = da[i]
    zt = stat.ztnb()
    zt.estimate(d, max_iter=10000)#, nlike=100)
    paras.append((zt.r, zt.p))
    #print >>estfile, zt.r, zt.p, sl[int(l * parts[ip])], d
    #ip += 1
  gtffile.close()
  bamfile.close()
  return paras, sl, data

def frame_bias(arr): # ['0', '1', '2', '01', '02', '12', '012']
  m = max(arr[0:codonSize])
  s = ''
  if m == 0 : return s
  for i in range(codonSize):
    if arr[i] == m : s += str(i)
  return s

################### Quality Control ###########################

# Calculate several distributions for quality control
def lendis(gtfpath, bampath, lens = [25,35], dis = [-40,20], ccds = False, maxNH = maxNH, minMapQ = minMapQ, secondary = secondary, minR = 1, m0 = True):
  #minR = 100
  bamfile = bam.bamfile(bampath, "rb")
  gtffile = open(gtfpath,'r')
  d = dis[1] - dis[0]
  dis1, dis2 = {}, {} #sum of reads near start or stop condons
  disf = {} # distribution of frame
  #fbias = {} # frame bias
  fbkeys = ['0', '1', '2', '01', '02', '12', '012']
  lendis = {} # reads length distribution
  if m0 : 
    dis1m0, dis2m0 = {}, {} #sum of reads near start or stop condons
    disfm0 = {} # distribution of frame
    lendism0 = {}
  for l in range(lens[0], lens[1]):
    lendis[l] = 0
    dis1[l] = [exp.readdict() for i in range(d)]
    dis2[l] = [exp.readdict() for i in range(d)]
    disf[l] = [exp.readdict() for i in range(codonSize)]
    #fbias[l] = {}
    if m0 :
      lendism0[l] = 0
      dis1m0[l] = [exp.readdict() for i in range(d)]
      dis2m0[l] = [exp.readdict() for i in range(d)]
      disfm0[l] = [exp.readdict() for i in range(codonSize)]
    #for k in fbkeys:
      #fbias[l][k] = 0
  for g in gtf.gtfgene_iter(gtffile):
    maxlen = 0
    mt = None
    for t in g.trans:
      if ccds and ('ccds_id' not in t.attr) : continue
      cds1, cds2 = t.cds_start(cdna = True), t.cds_stop(cdna = True) 
      try : cdslen = cds2 - cds1 
      except : continue
      if cdslen % 3 != 0 : continue
      if cdslen > maxlen : maxlen, mt = cdslen, t
    if mt is None : continue
    t = mt
    tl = t.cdna_length()
    cds1, cds2 = t.cds_start(cdna = True), t.cds_stop(cdna = True) - codonSize
    tdis1, tdis2 = {}, {}
    cnts = {}
    if m0 : 
      tdis1m0, tdis2m0 = {}, {}
      cntsm0 = {}
    for l in range(lens[0], lens[1]): 
      tdis1[l] = [0] * d
      tdis2[l] = [0] * d
      cnts[l] = [0] * tl
      if m0 :
        tdis1m0[l] = [0] * d
        tdis2m0[l] = [0] * d
        cntsm0[l] = [0] * tl
    tr = 0 # Total reads
    
    for r in bam.compatible_bam_iter(bamfile, t, mis = 2, maxNH = maxNH, minMapQ = minMapQ, secondary = secondary):
      ism0 = r.is_m0()
      #if m0 :
        #if r.strand == '+' and r.get_tag('MD')[0] == '0' : ism0 = True # mismatch at 0
        #if r.strand != '+' and r.get_tag('MD')[-1] == '0' : 
          #if not r.get_tag('MD')[-2].isdigit() : ism0 = True
      l = r.cdna_length()
      if l not in tdis1: continue # not in given length range
      i = t.cdna_pos(r.genome_pos(0)) # 5' end 
      if i is None : continue
      tr += 1
      if not ism0 : 
        lendis[l] += 1
        cnts[l][i] += 1
        ir = i - cds1
        if dis[0] <= ir < dis[1] : tdis1[l][ir - dis[0]] += 1
        ir = i - cds2
        if dis[0] <= ir < dis[1] : tdis2[l][ir - dis[0]] += 1
      else : 
        lendism0[l] += 1
        cntsm0[l][i] += 1
        ir = i - cds1
        if dis[0] <= ir < dis[1] : tdis1m0[l][ir - dis[0]] += 1
        ir = i - cds2
        if dis[0] <= ir < dis[1] : tdis2m0[l][ir - dis[0]] += 1
    if tr < minR : continue
    for l in tdis1: 
      for di in range(d): 
        dis1[l][di].record(tdis1[l][di])
        dis2[l][di].record(tdis2[l][di])
        if m0 :
          dis1m0[l][di].record(tdis1m0[l][di])
          dis2m0[l][di].record(tdis2m0[l][di])
      for i in range(cds1, cds2, codonSize):
        io = i-defOffset
        if io < 0 : continue
        for i2 in range(codonSize):
          disf[l][i2].record(cnts[l][io+i2]) #
          if m0 : disfm0[l][i2].record(cntsm0[l][io+i2])
        #s = frame_bias(cnts[l][io:io+codonSize])
        #if s != '' : fbias[l][s] += 1
  results = [lendis, dis1, dis2, disf]
  if m0 : results += [lendism0, dis1m0, dis2m0, disfm0]
  return results
def meandis(dis, geomean = False, add = 100):
  mdis = {}
  for l in dis: 
    mdis[l] = [0] * (d)
    for di in range(d): 
      if geomean: mdis[l][di] = dis[l][di].geomean(add = add)
      else : mdis[l][di] = dis[l][di].mean()
  return mdis
#print lendis
#print mdis1
#print mdis2
def lendisM0(gtfpath, bampath, lens = [26,35], dis = [-40,20], maxNH = maxNH, minMapQ = minMapQ, minR = 1):
  bamfile = bam.bamfile(bampath, "rb")
  gtffile = open(gtfpath,'r')
  d = dis[1] - dis[0]
  dis1, dis2 = {}, {} #sum of reads near start or stop condons
  disf = {} # distribution of frame
  fbias = {} # frame bias
  fbkeys = ['0', '1', '2', '01', '02', '12', '012']
  lendis = {} # reads length distribution
  for l in range(lens[0], lens[1]):
    lendis[l] = 0
    dis1[l] = [None] * d
    dis2[l] = [None] * d
    disf[l] = [None] * codonSize
    fbias[l] = {}
    for di in range(d): 
      dis1[l][di] = exp.readdict()
      dis2[l][di] = exp.readdict()
    for i in range(codonSize):
      disf[l][i] = exp.readdict()
    for k in fbkeys:
      fbias[l][k] = 0
  for g in gtf.gtfgene_iter(gtffile):
    maxlen = 0
    mt = None
    for t in g.trans:
      cds1, cds2 = t.cds_start(cdna = True), t.cds_stop(cdna = True) 
      try : cdslen = cds2 - cds1 
      except : continue
      if cdslen % 3 != 0 : continue
      if cdslen > maxlen : maxlen, mt = cdslen, t
    if mt is None : continue
    t = mt
    tl = t.cdna_length()
    cds1, cds2 = t.cds_start(cdna = True), t.cds_stop(cdna = True) - codonSize
    tdis1, tdis2 = {}, {}
    cnts = {}
    for l in range(lens[0], lens[1]): 
      tdis1[l] = [0] * d
      tdis2[l] = [0] * d
      cnts[l] = [0] * tl
    tr = 0
    
    for r in bam.compatible_bam_iter(bamfile, t, mis = 2, maxNH = maxNH, minMapQ = minMapQ):
      if r.get_tag('MD')[0] != '0' : continue ## Only reads with mismatch at 1st nt
      l = r.cdna_length()
      if l not in tdis1: continue
      i = t.cdna_pos(r.genome_pos(0)) # 5' end 
      if i is None : continue
      tr += 1
      lendis[l] += 1
      cnts[l][i] += 1
      ir = i - cds1
      if dis[0] <= ir < dis[1] : tdis1[l][ir - dis[0]] += 1
      ir = i - cds2
      if dis[0] <= ir < dis[1] : tdis2[l][ir - dis[0]] += 1
    if tr < minR : continue
    for l in tdis1: 
      for di in range(d): 
        dis1[l][di].record(tdis1[l][di])
        dis2[l][di].record(tdis2[l][di])
      for i in range(cds1, cds2, codonSize):
        io = i-defOffset
        if io < 0 : continue
        for i2 in range(codonSize):
          disf[l][i2].record(cnts[l][io+i2]) #
        s = frame_bias(cnts[l][io:io+codonSize])
        if s != '' : fbias[l][s] += 1
  return lendis, dis1, dis2, disf, fbias

def quality(arr, thresholds = [0.6, 0.7]): # Quality estimation by RPF frame distribution
  m = max(arr)
  for i, a in enumerate(arr):
    if a == m : 
      frame = i
      break
  txt = "%.2f " % (m)
  use = True
  if m < thresholds[0] : 
    txt += 'fail'
    use = False
  elif m < thresholds[1] : txt += 'pass'
  else : txt += 'good'
  return use, frame, txt

# Estimate RPF P site offset distance
def get_offset(arr, dis = [-40,20], frame = 0, defOffset = defOffset, flank = 6, tis = False): 
  a0 = [x for i, x in enumerate(arr) if (i + dis[0]) % codonSize == frame and i + dis[0] <= - defOffset - flank]
  a1 = [x for i, x in enumerate(arr) if (i + dis[0]) % codonSize == frame and i + dis[0] > - defOffset - flank]
  #a0.sort()
  #a1.sort()
  #a0m = 1.0 * sum(a0) / len(a0)
  #a1m = 1.0 * sum(a1) / len(a1)
  a0m = max(a0) # / 2.0
  a1m = max(a1) # / 2.0
  #th = int(((a0m + 1) * (a1m + 1)) ** 0.5)
  th = int((a0m + a1m) / 6.0) ###
  #print th
  for p in range(-defOffset - flank + 1, -defOffset + flank):
    if p % codonSize != frame : continue
    if tis :
      if arr[p - dis[0]] == a1m : return -p, th
    elif arr[p - dis[0]] > th : return -p, th
  return -p, th


offsetFuncStr = '''
def offset(length = 30):
  if length not in offdict : return None
  else : return offdict[length]
'''
def write_off_para(parafile, offdict): # Generate python code of offset function for parameter file
  print >>parafile, 'offdict =', offdict
  #print >>parafile, offsetFuncStr

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
def first_frames(arr, bin = 3):
  fs = []
  for i in range(0, len(arr), bin):
    fs.append(first(arr[i:i+bin])) # May out of index
  return fs
def frame_mean(fs):
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
def frame_test_n(arr, length, value, bin = 3, n = 1000, show = False): #expect no bias
  a = arr[:]
  c = 0
  for i in range(n):
    random.shuffle(a)
    fs = first_frames(a[0:length], bin)
    f = frame_mean(fs)
    if max(f) >= value: c += 1
    if show : print a, fs, f, c
    if i == 19 and c > 10 : return float(c) / 20
    if i == 99 and c > 20 : return float(c) / 100
  return float(c) / n
def frame_test_e(arr, length, value, expect, bin = 3, n = 1000, show = False):
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
    f = frame_mean(fs)
    m = max(f)
    if f[expect] < m and m >= value : c += 1
    if show : print fs, f, c
    if i == 19 and c > 10 : return float(c) / 20
    if i == 99 and c > 20 : return float(c) / 100 ###
  return float(c) / n
def bias(f):
  m = max(f)
  for i in range(len(f)):
    if f[i] == m : return i
  
def orf_read(cnts, cds1, cds2, nhead = 12, ntail = 18):
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

class region:
  def __init__(self, ers, start, stop, n, score = -1, p = 1):
    self.ers = ers # corrent enrichedregions object
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
    return region(self.ers, self.start, self.stop, self.n, self.score, self.p)
  def __str__(self):
    return "%d-%d n=%d score=%s p=%s" % (self.start, self.stop, self.n, str(self.score), str(self.p))
  def __repr__(self):
    return "Enriched region object " + str(self)
  def mapback(self, start = -1, stop = -1, nhead = -1, bin = -1):
    if nhead < 0 : nhead = self.ers.nhead
    if bin < 0 : bin = self.ers.bin
    if start < 0 : start = self.start 
    if stop < 0 : stop = self.stop
    start = start * bin + nhead
    stop = stop * bin + nhead
    return start, stop
  def binom_pvalue(self, l = -1):
    p = float(len(self))/self.ers.length
    if l > 0: p2 = float(len(self))/l
    else : p2 = p
    #print self.n, p
    return stat.binom_test(self.n, self.ers.total, p) / p2
  
  def get_frame(self, start = -1, stop = -1):
    if start < 0 : start = self.start 
    if stop < 0 : stop = self.stop
    f = frame_mean(self.ers.frames[start:stop])
    return f
  def local_pos(self, r, window = 20):
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
  
  def frame_check(self, n = 1000, window = 20, local = [0.0,0.5,1.0]):
    f = self.get_frame(self.start, self.stop)
    fm = max(f)
    tstart, tstop = self.mapback()
    p = frame_test_n(self.ers.cnts[tstart:tstop], tstop-tstart, fm, self.ers.bin, n = n)
    #print "region frame:",
    if p < 0.05 : 
      b = bias(f)
      #print "Bias =", b, "p =", p, "f =", f
    #else: print "No significant Bias, p =", p, "f =", f
    regbias = (f, p)
    locbias = {}
    l = len(self)
    if l > window :
      for r in local:
        (rs, re) = self.local_pos(r, window)
        (ts, te) = self.mapback(rs, re)
        rf = self.get_frame(rs, re)
        rfm = max(rf)
        rb = bias(rf)
        if p < 0.05 : 
          if round(rfm, 3) == round(rf[b], 3) : continue
          rpn = frame_test_n(self.ers.cnts[ts:te], te-ts, rfm, self.ers.bin, n = n)
          if rpn < 0.05 :
            rpe = frame_test_e(self.ers.cnts[tstart:tstop], te-ts, rfm, b, self.ers.bin, n = n)
          else : rpe = 1
          rp = max(rpn, rpe)
          if rp < 0.05 : tp = '1-1'
          else : tp = '1-0'
        else: 
          rp = frame_test_n(self.ers.cnts[ts:te], te-ts, rfm, self.ers.bin, n = n)
          if rp < 0.05 : tp = '0-1'
          else : tp = '0-0'
        if rp < 0.05 : 
          #print "Local bias: r =", r, ", rBias =", rb, ", rp =", rp, ", rf =", rf
          locbias[r] = (rf, rp)
    return (regbias, locbias)

class enrichedregions:
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
  
  def find_region(self):
    if self.total == 0 : raise StopIteration
    rarr = []
    rmax = region(self, 0, 1, 0, 0)
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
        r = region(self, bi1, bi2, s, score)
        if score > 0 : r.p = r.binom_pvalue()
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
      rmax = region(self, 0, 1, 0, 0)
      for rarr2 in rarr:
        for r in rarr2:
          if rmax < r : rmax = r
