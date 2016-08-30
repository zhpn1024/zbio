import math
from scipy.stats import nbinom, chisquare, chisqprob
logarr = [None] # log(N)
logsumarr = [0] # log(N!)
def logarr_ext(n, logarr = logarr, logsumarr = logsumarr):
  l = len(logarr)
  if l < n + 1 : 
    logarr += [None] * (n + 1 - l)
    logsumarr += [None] * (n + 1 - l)
    for i in range(l, n + 1):
      logarr[i] = math.log(i)
      #print i, logarr[i]
      logsumarr[i] = logsumarr[i-1] + logarr[i]
  return logarr, logsumarr
def data_count(data):
  total, cnt = 0, 0
  for k in data:
    total += k * data[k]
    cnt += data[k]
  return total, cnt
def mean_var(data):
  total, tsq, cnt = 0, 0, 0
  for k in data:
    total += k * data[k]
    cnt += data[k]
    tsq += k * k * data[k]
  mean = float(total) / cnt
  var = float(tsq) / cnt - mean ** 2
  return mean, var
def load_data(arr):
  data = {}
  for k in arr:
    if k not in data : data[k] = 0
    data[k] += 1
  return data
def fisher_method(ps):
  n = 0
  fs = 0
  for p in ps:
    if p is None : continue
    if p == 0 : return 0.0, -1 ###
    fs += - 2 * math.log(p)
    n += 1
  fp = chisqprob(fs, 2 * n)
  return fp, fs

def combination_log(n, k, logarr = logarr, show = False): # N choose K, log combination number, NATURAL LOG!
  if n < 0 : return None
  if k > n or k < 0: return None #None is log(0)
  #if k * 2 > n: k = n - k
  #lpr = 0.0
  nk = n - k
  logarr_ext(n, logarr = logarr)
  lpr = logsumarr[n] - logsumarr[k] - logsumarr[nk]
  '''for i in range(n, nk, -1):
    lpr += logarr[i]
  for i in range(k, 0, -1):
    lpr -= logarr[i]'''
  return lpr

def ACprob(x, y, r = 1): # p(y|x) = r^y C(x+y,y) / (1+r)^(x+y+1) , r = N2/N1 by Audic and Claverie 
  lp = combination_log(x+y, x)
  lp += math.log(r) * y
  lp -= math.log(1 + r) * (1 + x + y)
  return math.exp(lp)
def ACtest(x, y, r = 1, alt = 'auto', double = True): # Diff expression test by Audic and Claverie 
  if alt in ('g', 'greater') : n1, n2, nr = x, y, r  # if x > y ?
  elif alt in ('l', 'less') : n1, n2, nr = y, x, 1.0/r
  elif x * r < y : n1, n2, nr = y, x, 1.0/r # n2 is smaller than n1
  else : n1, n2, nr = x, y, r
  pv = 0
  for i in range(n2 + 1):
    pv += ACprob(n1, i, nr)
    if double and pv >= 0.5 : return 1
    #print pv
  if double : pv *= 2 # Two tailed
  if pv > 1 : pv = 1
  return pv

def FCtest(x, y, r = 1, fc = 1.5, alt = 'auto', double = True): # Test whether exp diff > fold change (binom test). For data with no replicate|dispersion 
  n = x + y
  if y > 0 : fcr = 1.0 * x * r / y
  else : fcr = fc + 1 ###
  if alt in ('g', 'greater') and fcr <= fc : return 1
  elif alt in ('l', 'less') and fcr >= 1.0 / fc : return 1
  if 1.0 / fc <= fcr <= fc : return 1
  if fcr > fc : 
    p = fc / (fc + r)
    pv = binom_test(x, n, p = p, alt = "g")
  elif fcr < 1.0 / fc : 
    p = 1.0 / (1 + r * fc)
    pv = binom_test(x, n, p = p, alt = "l")
  if double : pv *= 2  # Doubling the smaller tail
  if pv > 1 : pv = 1
  return pv

def hypergeo0(N, K, n, k):
  p = 1.0
  nk = n - k
  Nk = N - k
  NK = N - K
  NKnk1 = N - K - n + k + 1
  Nn1 = N - n + 1
  #pmax = 0
  #pmin = 1
  for i in range(k):
    p *= float(K - i) / (N - i)
    p *= float(n - i) / (k - i)
    if i < nk: # j = nk - i - 1
      p *= float(NK - i) / (Nk - i)
      #p *= float(NKnk1 + i) / (Nn1 + i)
    #if p > pmax: pmax = p
    #if p < pmin: pmin = p
  #print p, pmax, pmin
  for i in range(k, nk):
    p *= float(NK - i) / (Nk - i)
    #p *= float(NKnk1 + i) / (Nn1 + i)
  return p

def hypergeo1(N, K, n, k):
  p = 1.0
  NK = N - K
  nk = n - k
  Knk = K + n - k 
  Nnk = N - n + k
  for i in range(nk):
    p *= float(NK - i) / (N - i)
    p *= float(n - i) / (nk - i)
    if i < k:
      p *= float(K - i) / (Nnk - i)
  #print p
  for i in range(nk, k):
    p *= float(K - i) / (Nnk - i)
  #for i in range(nk, n):
    #p *= float(Knk - i) / (N - i)
    #p *= float(n - i) / (n - i)
  return p
    
def hypergeo(N, K, n, k):
  if k >= n / 2 : return hypergeo0(N, K, n, k)
  else: return hypergeo1(N, K, n, k)
  
def binomial(k, n, p = 0.5, show=False):
  if k > n or k < 0: return 0
  if p < 0 : p = 0
  if p > 1 : p = 1
  if k * 2 > n:
    k = n - k
    p = 1 - p
  q = 1 - p
  pr = 1.0
  nk = n - k
  if k == 0 : return q ** n
  t = nk / k
  if t * k < nk: t += 1
  qi = 0
  for i in range(k):
    pi = 1.0
    pi *= n - i
    pi /= k - i
    pi *= p
    for ti in range(t):
      if qi >= nk: break
      pi *= q
      qi += 1
    pr *= pi
    if show: print pr
  return pr
def binom_log(k, n, p = 0.5, logarr = logarr, show = False): #log probability value
  if n < 0 : return None
  if k > n or k < 0: return None #None is log(0)
  if p < 0 : p = 0
  if p > 1 : p = 1
  if k * 2 > n:
    k = n - k
    p = 1 - p
  q = 1 - p
  if p == 0 :
    if k == 0 : return 0
    else : return None
  elif q == 0 : return None
  lpr = 0.0
  lp = math.log(p)
  lq = math.log(q)
  nk = n - k
  #if k == 0 : return  lq * n
  lpr += lp * k + lq * nk
  lpr += combination_log(n, k)
  '''logarr_ext(n, logarr = logarr)
  for i in range(k):
    lpr += logarr[n - i] - logarr[k - i]'''
  return lpr
def binom_test(k, n, p = 0.5, alt = "g", log = True, logarr = logarr, show=False): # No two sided yet!
  if not log : return binomTest0(k, n, p, alt, show) # if log, p are calculated with log10 values
  #logarr = [None] * (n + 1)
  #logarr_ext(n, logarr = logarr)
  lpk = binom_log(k, n, p, logarr)
  if show : print lpk
  if lpk is None : 
    if alt[0] == 'g' and p >= 1: return 1
    if alt[0] != 'g' and p <= 0: return 1
    return 0
  elif lpk == 0 : return 1
  #if k == 0 and alt[0] == 'g' : return 1
  #if k == n and alt[0] != 'g' : return 1
  pv = math.exp(lpk) ### log reverse
  q = 1 - p
  lp = math.log(p)
  lq = math.log(q)
  if alt[0] == 'g':
    for i in range(k, n):
      lpk += lp + logarr[n - i] - lq - logarr[i + 1]
      #r = p * (n - i) / q / (i + 1)
      #pk *= r
      pv += math.exp(lpk)
  else:
    for i in range(k, 0, -1):
      lpk += lq + logarr[i] - lp - logarr[n - i + 1]
      #r = q * i / p / (n - i + 1)
      #pk *= r
      pv += math.exp(lpk)
  return pv
def binomTest0(k, n, p = 0.5, alt = "g", show=False): # No two sided yet!
  pk = binomial(k, n, p)
  if show : print pk
  if pk == 1 : return 1
  if pk == 0 : 
    if alt[0] == 'g' and p >= 1: return 1
    if alt[0] != 'g' and p <= 0: return 1
    return 0
  #if k == 0 and alt[0] == 'g' : return 1
  #if k == n and alt[0] != 'g' : return 1
  pv = pk
  q = 1 - p
  if alt[0] == 'g': 
    for i in range(k, n):
      #r = 1.0
      r = p * (n - i) / q / (i + 1)
      pk *= r
      pv += pk
  else:
    for i in range(k, 0, -1):
      r = q * i / p / (n - i + 1)
      pk *= r
      pv += pk
  return pv
class negbinom: #Number of 'failures' before 'r' 'successes' with success probability 'p'
  rMax = 1e8
  rMin = 1e-8
  Delta = 1e-8
  def __init__(self, r = 1.0, p = 0.5):
    self.p = p
    self.r = r
  @property
  def q(self):
    return 1 - self.p
  def expect(self):
    return self.q * self.r / self.p
  def variance(self):
    return self.expect() / self.p
  def logpmf(self, k = 0):
    return nbinom.logpmf(k, self.r, self.p)
  def pmf(self, k = 0):
    return nbinom.pmf(k, self.r, self.p)
  def cdf(self, k = 0):
    return nbinom.cdf(k, self.r, self.p)
  def pvalue(self, k = 0):
    p = 1 - nbinom.cdf(k-1, self.r, self.p)
    if p > self.Delta : return p
    p = nbinom.pmf(k, self.r, self.p)
    if p == 0 : return p
    ka = k + 1
    pa = nbinom.pmf(ka, self.r, self.p)
    p += pa
    while pa / p > self.Delta:
      #p += pa
      ka += 1
      pa = nbinom.pmf(ka, self.r, self.p)
      p += pa
    #p += pa
    return p
      
  def estimate(self, data): #data dict value:counts
    total, cnt = data_count(data)
    rmax, rmin = self.rMax, self.rMin
    rmid = math.sqrt(rmax * rmin)
    while (rmax - rmin) / rmid >= self.Delta:
      #print rmax, rmin
      score = self.r_log_like_score(data, rmid)
      #print score, cnt, rmid
      if score > 0 : rmin = rmid
      elif score < 0 : rmax = rmid
      else : break
      #rmid = (rmax + rmin) / 2
      rmid = math.sqrt(rmax * rmin)
    self.r = rmid
    #self.p = total / (self.r * cnt + total)
    self.p = self.r / (self.r + 1.0*total/cnt)
    return self.r, self.p

  def log_likelihood(self, data):
    score = 0.0
    for k in data:
      score += data[k] * self.logpmf(k)
    return score
  def r_log_like_score(self, data, r = -1):
    if r < 0 : r = self.r
    total, cnt = data_count(data)
    #dr = scipy.special.digamma(r)
    #score = math.log(self.q) - dr
    s1, d = 0, 0
    for i in range(max(data.keys()) + 1):
      #d += 1.0 / (r + i)
      if i in data : s1 += d * data[i]
      d += 1.0 / (r + i)
    score = s1 / cnt + math.log(r / (r + 1.0*total/cnt))
    return score
  def expected(self, k, size):
    p = self.pmf(k)
    return size * p
  def estimate_truncated(self, data, size, max_iter = 1e4, nlike = 10, report = False): # Not finished!
    total, cnt = data_count(data)
    lastllh = 0
    i = 0
    km = max(data)
    #d = {}
    #for k in data : d[k] = data[k]
    for j in range(int(max_iter)) :
      ps = 0
      for k in data : ps += self.pmf(k)
      #size = round(total / ps)
      k = 0
      d = {}
      #ek = self.expected(k, size)
      while k <= km or d[k-1] >= 1 : 
        if k in data : d[k] = data[k]
        else : d[k] = round(self.expected(k, size))
        k += 1
        #if k > km : print k, d[k]
      #ez = self.expected_zeros(cnt)
      #data[0] = ez
      self.estimate(d)
      #data[0] = 0
      i += 1
      if i < nlike : continue
      i = 0
      llh = self.log_likelihood(d)
      diff = abs(2 * (llh - lastllh) / (llh + lastllh) / nlike)
      if diff < self.Delta : break
      if report : 
        print diff, llh, self.r, self.p, k, size,
        for i in range(15) : print '%d:%d' % (i, d[i]),
        print '...'
      lastllh = llh
    return self.r, self.p
  def estimate_by_012(self, n0, n1, n2, start = 0) : 
    r1 = (start + 1) * float(n1) / n0 
    r2 =  (start + 2) * float(n2) / n1
    p1 = r2 - r1
    self.p = 1 - p1 # 
    self.r = (r1 - start * p1) / p1
    return self.r, self.p
  def fit_lower(self, data, nmax = 40, pth = 0.01, start = 0) : 
    r, p = self.estimate_by_012(data[start], data[start+1], data[start+2], start = start)#, nlike=100)
    if r > 0 and 0 < p < 0.9 : 
      #self.r, self.p = r, p
      return r, p
    else : 
      import numpy as np
      from scipy.optimize import leastsq
      def res(p, x, y) : 
        return p[0] * x + p[1] - y
      m = max(data)
      if m > nmax : m = nmax # max
      s = data[start] + data[start+1] + data[start+2]
      y = [(start + 1) * float(data[start+1]) / data[start], (start + 2) * float(data[start+2]) / data[start+1]]
      w = [(data[start+1] * data[start])**0.5, (data[start+1] * data[start+2])**0.5]
      for i in range(start+3, m): 
        y.append(i * float(data.value(i)) / data.value(i-1))
        w.append((data[i] * data[i-1])**0.5)
        wm = min([wi for wi in w if wi > 0])
        w2 = [int(round(wi / wm)) for wi in w]
        s += data.value(i)
        xl = sum([[j] * w2[j] for j in range(i)], [])
        yl = sum([[y[j]] * w2[j] for j in range(i)], [])
        xa = np.array(xl, dtype = float)
        ya = np.array(yl, dtype = float)
        rst = leastsq(res, [1,1], args=(xa, ya))
        p = 1 - rst[0][0]
        r = rst[0][1] / rst[0][0]
        print i, (r, p), len(xl), w#, y
        if r <= 0 or p >= 1 or p <= 0 : continue
        if 0.2 < p < 1 : break
        self.r, self.p = r, p
        #f1 = self.pmf(i) / self.cdf(i)
        f1 = self.pmf(i+1) / (self.cdf(i+1) - self.cdf(start-1))
        f2 = self.pmf(i+2) / (self.cdf(i+2) - self.cdf(start-1))
        #pv1 = binom_test(data.value(i), s, p = f1, alt = 'g')
        pv1 = binom_test(data.value(i+1), s+data.value(i+1), p = f1, alt = 'g')
        pv2 = binom_test(data.value(i+2), s+data.value(i+1)+data.value(i+2), p = f2, alt = 'g')
        #for j in range(i) : print j, self.pmf(j),
        print data.value(i+1), s+data.value(i+1), f1, 'pv1 =', pv1
        print data.value(i+2), s+data.value(i+1)+data.value(i+2), f2, 'pv2 =', pv2
        #fp, fs = fisher_method([pv1, pv2])
        fp = max(pv1, pv2)
        #print fp
        if fp < pth : break #and pv2 < pth : break
      self.r, self.p = r, p
      return r, p
  def fit_linear(self, data, poisson = False, maxi = 40, minc = 3, start = 0, total = None, show = False) : 
    if poisson : N2 = 0
    else : 
      if total is None : total = sum(data.values())
      N2 = 2.0 / total
    rate, rvar = {}, {}
    for i in data : 
      if i - 1 < start : continue
      if i-1 not in data or data[i-1] <= 0 : continue
      if i not in data or data[i] <= 0 : continue
      if data[i] < minc : break
      #if data.value(i-1) == 0 or data.value(i) == 0 : continue
      if i >= maxi : break
      rate[i] = float(i) * data[i] / data[i-1]
      rvar[i] = rate[i] ** 2 * (1.0 / data[i] + 1.0 / data[i] - N2)
    #print 'rate =', rate
    #print 'rvar =', rvar
    c1, c2, c3, c4, c5 = 0,0,0,0,0
    a, b = {}, {}
    r, p = {}, {}
    ia = rate.keys()
    ia.sort()
    for n, i in enumerate(ia) : 
      c1 += 1 / rvar[i]
      c2 += i / rvar[i]
      c3 += i * i / rvar[i]
      c4 += rate[i] / rvar[i]
      c5 += i * rate[i] / rvar[i]
      if n > 0 : 
        a[i] = (c3 * c4 - c2 * c5) / (c1 * c3 - c2 * c2)
        b[i] = (c2 * c4 - c1 * c5) / (c2 * c2 - c1 * c3)
        p[i] = 1 - b[i]
        r[i] = a[i] / b[i] + 1
        if show : print i, a[i], b[i], r[i], p[i]
    self.r, self.p = r[i], p[i]
    return a, b, r, p, rate, rvar
    r, p = self.estimate_by_012(data[start], data[start+1], data[start+2], start = start)#, nlike=100)
    if r > 0 and 0 < p < 0.9 : 
      #self.r, self.p = r, p
      return r, p
    else : 
      import numpy as np
      from scipy.optimize import leastsq
      def res(p, x, y) : 
        return p[0] * x + p[1] - y
      m = max(data)
      if m > nmax : m = nmax # max
      s = data[start] + data[start+1] + data[start+2]
      y = [(start + 1) * float(data[start+1]) / data[start], (start + 2) * float(data[start+2]) / data[start+1]]
      w = [(data[start+1] * data[start])**0.5, (data[start+1] * data[start+2])**0.5]
      for i in range(start+3, m): 
        y.append(i * float(data.value(i)) / data.value(i-1))
        w.append((data[i] * data[i-1])**0.5)
        wm = min([wi for wi in w if wi > 0])
        w2 = [int(round(wi / wm)) for wi in w]
        s += data.value(i)
        xl = sum([[j] * w2[j] for j in range(i)], [])
        yl = sum([[y[j]] * w2[j] for j in range(i)], [])
        xa = np.array(xl, dtype = float)
        ya = np.array(yl, dtype = float)
        rst = leastsq(res, [1,1], args=(xa, ya))
        p = 1 - rst[0][0]
        r = rst[0][1] / rst[0][0]
        print i, (r, p), len(xl), w#, y
        if r <= 0 or p >= 1 or p <= 0 : continue
        if 0.2 < p < 1 : break
        self.r, self.p = r, p
        #f1 = self.pmf(i) / self.cdf(i)
        f1 = self.pmf(i+1) / (self.cdf(i+1) - self.cdf(start-1))
        f2 = self.pmf(i+2) / (self.cdf(i+2) - self.cdf(start-1))
        #pv1 = binom_test(data.value(i), s, p = f1, alt = 'g')
        pv1 = binom_test(data.value(i+1), s+data.value(i+1), p = f1, alt = 'g')
        pv2 = binom_test(data.value(i+2), s+data.value(i+1)+data.value(i+2), p = f2, alt = 'g')
        #for j in range(i) : print j, self.pmf(j),
        print data.value(i+1), s+data.value(i+1), f1, 'pv1 =', pv1
        print data.value(i+2), s+data.value(i+1)+data.value(i+2), f2, 'pv2 =', pv2
        #fp, fs = fisher_method([pv1, pv2])
        fp = max(pv1, pv2)
        #print fp
        if fp < pth : break #and pv2 < pth : break
      self.r, self.p = r, p
      return r, p
  def chisquare_test(self, data):
    total, cnt = data_count(data)
    obs, exs = [], []
    ob, ex = 0, 0
    i0 = 0
    for i in range(max(data.keys()) + 1) :
      if i in data : ob += data[i]
      ex += cnt * self.pmf(i)
      if ex >= 5 : 
        obs.append(ob)
        exs.append(ex)
        ob, ex = 0, 0
        i0 = i + 1
    ex = cnt * self.pvalue(i0)
    #obs.append(ob)
    #exs.append(ex)
    obs[-1] += ob
    exs[-1] += ex
    print obs, exs, len(obs) - 1, len(exs), sum(obs), sum(exs)
    return chisquare(obs, exs)
  
class ztnb(negbinom):
  def logpmf(self, k = 1):
    if k < 1 : return negbinom.pmf(self, -1)
    p0 = negbinom.pmf(self, 0)
    lp = negbinom.logpmf(self, k)
    return lp - math.log(1 - p0)
  def pmf(self, k = 1):
    if k < 1 : return 0
    p0 = negbinom.pmf(self, 0)
    p = negbinom.pmf(self, k)
    return p / (1- p0)
  def expected_zeros(self, size):
    p0 = negbinom.pmf(self, 0)
    return size * p0 / (1 - p0)
  def estimate(self, data, max_iter = 1e4, nlike = 10, report = False):
    d1 = {}
    for i in data : 
      if i > 0 : d1[i] = data[i]
    total, cnt = data_count(d1)
    lastllh = 0
    i = 0
    for j in range(int(max_iter)) :
      ez = self.expected_zeros(cnt)
      d1[0] = ez
      negbinom.estimate(self, d1)
      d1[0] = 0
      i += 1
      if i < nlike : continue
      i = 0
      llh = self.log_likelihood(d1)
      d = abs(2 * (llh - lastllh) / (llh + lastllh) / nlike)
      if d < self.Delta : break
      if report : print d, llh, self.r, self.p, ez
      lastllh = llh
    return self.r, self.p
  def pvalue(self, k = 1):
    if k <= 1 : return 1
    p0 = negbinom.pmf(self, 0)
    p = negbinom.pvalue(self, k)
    return p / (1- p0)
  def expect(self):
    p0 = negbinom.pmf(self, 0)
    return negbinom.expect(self) / (1 - p0)
  def variance(self):
    nb = negbinom(self.r, self.p)
    p0 = nb.pmf(0)
    return (nb.variance() + nb.expect() ** 2) / (1 - p0) - self.expect() ** 2
    #return self.expect / self.p

class poisson: # Poisson distribution
  lMax = 1e8
  lMin = 1e-8
  Delta = 1e-8
  def __init__(self, l = 1.0):
    self.l = l
  def expect(self):
    return self.l
  def variance(self):
    return self.l
  def logpmf(self, k = 0, logarr = logarr):
    logarr_ext(k, logarr = logarr)
    lpr = k * math.log(self.l) - self.l
    lpr -= logsumarr[k]
    #for i in range(k):
      #lpr -= logarr[k - i]
    return lpr
  def pmf(self, k = 0):
    return math.exp(self.logpmf(k))
  def cdf(self, k = 0, logarr = logarr):
    if k < 0 : return 0
    logarr_ext(k, logarr = logarr)
    lpr = self.logpmf(0)
    logl = math.log(self.l)
    cdf = math.exp(lpr)
    for i in range(1, k + 1):
      lpr += logl - logarr[i]
      cdf += math.exp(lpr)
    return cdf
  def pvalue(self, k = 0):
    return 1 - self.cdf(k-1)
  def estimate(self, data): #data dict value:counts
    total, cnt = data_count(data)
    self.l = float(total) / cnt
    return self.l
  def log_likelihood(self, data):
    score = 0.0
    for k in data:
      score += data[k] * self.logpmf(k)
    return score
  def r_log_like_score(self, data, l = -1):
    if l < 0 : l = self.l
    total, cnt = data_count(data)
    return float(total) / l / cnt -1
  def chisquare_test(self, data):
    total, cnt = data_count(data)
    obs, exs = [], []
    ob, ex = 0, 0
    i0 = 0
    for i in range(max(data.keys()) + 1) :
      if i in data : ob += data[i]
      ex += cnt * self.pmf(i)
      if ex >= 5 : 
        obs.append(ob)
        exs.append(ex)
        ob, ex = 0, 0
        i0 = i + 1
    ex = cnt * self.pvalue(i0)
    obs[-1] += ob
    exs[-1] += ex
    print obs, exs, len(obs) - 1, len(exs), sum(obs), sum(exs)
    return chisquare(obs, exs)

class ztpoisson: #Zero truncated poisson distribution
  def expect(self):
    p0 = poisson.pmf(self, 0)
    return self.l / (1 - p0)
  def variance(self):
    p0 = poisson.pmf(self, 0)
    return (self.l + self.l ** 2) / (1 - p0) - self.expect() ** 2
  def logpmf(self, k = 0):
    if k < 1 : return negbinom.pmf(self, -1)
    p0 = math.exp(poisson.logpmf(self, 0))
    lp = poisson.logpmf(self, k)
    return lp - math.log(1 - p0)
  #def pmf(self, k = 0):
    #return math.exp(self.logpmf(k))
  def cdf(self, k = 1, logarr = logarr):
    if k <= 0 : return 0
    logarr_ext(k, logarr = logarr)
    lpr = self.logpmf(1)
    logl = math.log(self.l)
    cdf = math.exp(lpr)
    for i in range(2, k + 1):
      lpr += logl - logarr[i]
      cdf += math.exp(lpr)
    return cdf
  #def pvalue(self, k = 0):
    #return 1 - self.cdf(k - 1)
  def estimate(self, data, max_iter = 1e4, nlike = 10): ######
    total, cnt = data_count(data)
    lastllh = 0
    i = 0
    for j in range(max_iter) :
      ez = self.expected_zeros(cnt)
      data[0] = ez
      negbinom.estimate(self, data)
      data[0] = 0
      i += 1
      if i < nlike : continue
      i = 0
      llh = self.log_likelihood(data)
      d = abs(2 * (llh - lastllh) / (llh + lastllh) / nlike)
      if d < self.Delta : break
      print d, llh, self.r, self.p, ez
      lastllh = llh
    return self.r, self.p
  def log_likelihood(self, data):
    score = 0.0
    for k in data:
      score += data[k] * self.logpmf(k)
    return score
  def r_log_like_score(self, data, r = -1):
    if r < 0 : r = self.r
    total, cnt = data_count(data)
    #dr = scipy.special.digamma(r)
    #score = math.log(self.q) - dr
    s1, d = 0, 0
    for i in range(max(data.keys()) + 1):
      #d += 1.0 / (r + i)
      if i in data : s1 += d * data[i]
      d += 1.0 / (r + i)
    score = s1 / cnt + math.log(r / (r + 1.0*total/cnt))
    return score
  def chisquare_test(self, data):
    total, cnt = data_count(data)
    obs, exs = [], []
    ob, ex = 0, 0
    i0 = 0
    for i in range(max(data.keys()) + 1) :
      if i in data : ob += data[i]
      ex += cnt * self.pmf(i)
      if ex >= 5 : 
        obs.append(ob)
        exs.append(ex)
        ob, ex = 0, 0
        i0 = i + 1
    ex = cnt * self.pvalue(i0)
    #obs.append(ob)
    #exs.append(ex)
    obs[-1] += ob
    exs[-1] += ex
    print obs, exs, len(obs) - 1, len(exs), sum(obs), sum(exs)
    return chisquare(obs, exs)

class rankSumTiesExact: # exact rank sum test of x < y (one-sided) for ribo
  def __init__(self, x, y, show = False):
    self.a = list(x) + list(y)
    self.x = list(x)
    self.N, self.n = len(self.a), len(x)
    self.cd = countDict(self.a) # count dict
    self.xcd = countDict(x)
    self.rd = rankDict(self.cd)
    self.xrs =  self.rankSum(self.xcd)
    self.ks = self.cd.keys()
    self.ks.sort()
    self.l = len(self.ks)
    self.vs = [self.count(i) for i in range(self.l)]
    if show : 
      print self.cd, self.xcd
      print self.tieRatio()
      print self.complexity()
      #print self.numStats()
  def copy(self, other):
    self.a = other.a
    self.x = other.x
    self.N, self.n = other.N, other.n
    self.cd = other.cd # count dict
    self.xcd = other.xcd
    self.rd = other.rd
    self.xrs =  other.xrs
    self.ks = other.ks
    self.l = other.l
    self.vs = other.vs
  def complexity(self):
    logarr_ext(max(self.cd.values()) + 1)
    complog = sum([logarr[self.count(i)+1] for i in range(self.l)])
    return complog
  def shuffleTest(self, n = 1000, show = False):
    import random
    c = 0
    a = self.a[:]
    for i in range(n):
      random.shuffle(a)
      rcd = countDict(a[0:self.n])
      rs =  self.rankSum(rcd)
      if rs >= self.xrs: c += 1
      if show : print a, rs, xrs, c
    #if i == 19 and c > 10 : return float(c) / 20
    #if i == 99 and c > 20 : return float(c) / 100
    return float(c) / n
  def count(self, i):
    return self.cd[self.ks[i]]
  def rank(self, i):
    return self.rd[self.ks[i]]
  def numStats(self):
    logarr_ext(max(self.cd.values()) + 1)
    return self._numStats(self.N, self.n, 0)
  def _numStats(self, N, n, i):
    n1 = N - n
    if n < 0 or n1 < 0 : return 0
    if n == 0 or n1 == 0 : return 1
    if i == self.l - 2 : return min(n, n1) + 1
    if max(self.vs[i:]) == 1 : return math.exp(combination_log(N,n))
    t1, t2 = self.count(i), N - self.count(i)
    if t1 >= n >= t2 or t1 >= n1 >= t2 : 
      d = math.exp(sum([logarr[self.count(k)+1] for k in range(i+1, self.l)]))
      return d
      #print 'multi all', d, N, n, i
    s = 0
    for ni in range(max(0, n-t1),min(self.count(i), n) + 1) :
      d = self._numStats(N - self.count(i), n - ni, i+1)
      s += d
      #print d,s,N - ni, n - ni, i+1
    return s
  def numStatsRaw(self):
    p = float(self.n) / self.N
    r = p * (1-p)
    s = 1
    for i in range(self.l) : 
      #var = self.count(i) * r
      #sd = var ** 0.5
      #s *= sd * 2
      s *= self.count(i) + 1
    return s
  def tieRatio(self):
    s = sum([v**3 - v for v in self.cd.values() if v > 1])
    return s / float(self.N ** 3 - self.N)
  def test(self, th = 20, delta = 1e-4):
    if self.n <= th or self.N - self.n <= th : return self.fastTest(delta = delta)
    if self.complexity() <= 2 * th * logarr[2] : return self.fastTest(delta = delta)
    return self.mwtest()
  def mwtest(self, use_continuity = True, show = False):
    from scipy.stats import norm
    n1, n2, n = self.n, self.N - self.n, self.N
    mu = n1 * (n + 1) / 2.0
    s = sum([v**3 - v for v in self.cd.values() if v > 1])
    if show : print 'effect size', n - s / float(n*(n-1))
    var = n1 * n2 * (n ** 3 - n - s) / 12.0 / n / (n - 1)
    #print n, s, var, n1, self.cd
    if use_continuity : z = (self.xrs - mu - 0.5) / var ** 0.5
    else : z = (self.xrs - mu) / var ** 0.5
    p = norm.sf(abs(z))
    if show : print self.xrs, mu, self.xrs-mu+n1*(n1+1)/2.0, var, z, p
    if z >= 0 : return p
    else : return 1 - p # one sided
  def isExtreme(self, rs, twotailed = False, alt = 'g', delta = 1e-5): # if the given rank sum is farther than xrs
    if twotailed : 
      if not hasattr(self, 'mu'):
        n1, n2, n = self.n, self.N - self.n, self.N
        self.mu = n1 * (n + 1) / 2.0
        self.th = abs(self.xrs - self.mu)
      return abs(prs - self.mu) >= self.th - delta
    elif alt == 'g' : return rs >= self.xrs - delta
    else : return rs <= self.xrs + delta
  def exactTest(self, twotailed = False, show = False): 
    pval = 0
    if twotailed : 
      n1, n2, n = self.n, self.N - self.n, self.N
      mu = n1 * (n + 1) / 2.0
      th = abs(self.xrs - mu)
    for pcd in self.multiHypergeoIter(self.n, (), 0, self.N):
      prs =  self.rankSum(pcd)
      #if show : print prs, pcd
      if twotailed : 
        if abs(prs - mu) < th - 0.0001 : continue
      elif prs < self.xrs - 0.0001 : continue
      p = multiHypergeoProb(pcd, self.cd)
      pval += p
      if show : print prs, pcd, p
    if pval > 1.0 : pval = 1.0
    return pval
  def fastTest(self, show = False, delta = 1e-4):
    self.pval = 0
    self.lp0 = - combination_log(self.N, self.n)
    for pcd in self.multiHypergeoFastIter(self.n, (), 0, self.N, 0, delta=delta):
      p = multiHypergeoMergeProb(pcd, self.cd, self.N, self.n, lp0 = self.lp0) # merged probability
      self.pval += p
      if show : print p, pcd
    if self.pval > 1.0 : self.pval = 1.0
    return self.pval
  def multiHypergeoIter(self, n, vs, i, N):
    if n == 0 : yield zipDict(self.ks, vs) # no more, return
    elif n > N : return # impossible
    elif i + 1 == len(self.ks) : # the last key
      if n <= self.cd[self.ks[i]] : yield zipDict(self.ks, vs + (n,))
      return
    else :
      N1 = N - self.count(i)
      for j in range(max(0, n-N1), min(n, self.count(i)) + 1):
        for pcd in self.multiHypergeoIter(n-j, vs+(j,), i+1, N1):
          yield pcd
  def multiHypergeoFastIter(self, n, vs, i, N, rs, delta=1e-4): # fast iter, only select rank sum higher than x rank sum conditions
    #if i <= 2 : print n, vs, i, N, rs, self.pval
    if n > N : return # impossible
    #rsu, rsd = self.rankSumUpDownLimit(vs, rs)
    #print rsu, rsd, rs, self.xrs
    #if rsu < self.xrs : return # lower conditions
    #elif rsd >= self.xrs : yield zipDict(self.ks, vs) # good enough
    elif n == 0 : yield zipDict(self.ks, vs) # no more, return
    elif i + 1 == len(self.ks) : # the last key
      if n <= self.count(i) : yield zipDict(self.ks, vs + (n,))
      return
    else :
      N1 = N - self.count(i)
      jd, ju = max(0, n-N1), min(n, self.count(i)) + 1
      rsu, rsd = {}, {}
      j1, j2 = jd, ju
      while j1 < j2 : # looking for threshold of all lower conditions
        j = (j1 + j2) / 2 # + 1 for upper int
        rsu[j] = self.rankSumUpLimit(vs+(j,), rs+j*self.rank(i))
        if rsu[j] < self.xrs - 0.0001 : j2 = j # only keep possible conditions
        else : j1 = j + 1
      ju = j1
      j1, j2 = jd, ju
      while j1 < j2 : # looking for threshold of all higher conditions
        j = (j1 + j2) / 2 # + 1 for upper int
        rsd[j] = self.rankSumDownLimit(vs+(j,), rs+j*self.rank(i))
        if rsd[j] >= self.xrs - 0.0001 : j1 = j + 1 # only keep partial possible conditions
        else : j2 = j
      for j in range(jd, j1) : yield zipDict(self.ks, vs+(j,))
      #allgood = False
      #print j1, jd, ju
      ps = [multiHypergeoMergeProb(zipDict(self.ks, vs+(j,)), self.cd, self.N, self.n, lp0=self.lp0) for j in range(j1, ju)]
      for jo in orderIter(ps, reverse = True): #range(j1, ju):
        j = j1 + jo
        '''if allgood : yield zipDict(self.ks, vs+(j,))
        else : 
          vsj, rsj = vs+(j,), rs+j*self.rank(i)
          rsd = self.rankSumDownLimit(vsj, rsj)
          if rsd >= self.xrs : 
            allgood = True
            yield zipDict(self.ks, vsj)
          else :'''
        vsj, rsj = vs+(j,), rs+j*self.rank(i)
        #p = multiHypergeoMergeProb(zipDict(self.ks, vsj), self.cd, self.N, self.n, lp0 = self.lp0) 
        if self.pval > 0 and ps[jo] / self.pval < delta : # break # whether prob. high enough to calculate in detail
          #if j not in rsu : rsu[j] = self.rankSumUpLimit(vs+(j,), rs+j*self.rank(i))
          #if j not in rsd : rsd[j] = self.rankSumDownLimit(vs+(j,), rs+j*self.rank(i))
          #ratio = (rsu[j] + 1 - self.xrs) / (rsu[j] + 1 - rsd[j])
          #self.pval += ps[jo] * ratio # / 2 # * ratio # 1st order proximation
          self.pval += ps[jo] / 2
        else :
          for pcd in self.multiHypergeoFastIter(n-j, vsj, i+1, N1, rsj, delta=delta):
            yield pcd
  def rankSum(self, cd):
    return sum([self.rd[k] * cd[k] for k in cd])
  def rankSumUpLimit(self, vs, rs): # rank sum up & down limit
    n, i = self.n - sum(vs), len(self.ks) - 1
    while n >= 1 :
      if n > self.count(i) : d = self.count(i)
      else : d = n
      rs += d * self.rank(i)
      n -= d
      i -= 1
    return rs
  def rankSumDownLimit(self, vs, rs): # rank sum up & down limit
    n, i = self.n - sum(vs), len(vs)
    while n >= 1 : 
      #print n, i
      if n > self.count(i) : d = self.count(i)
      else : d = n
      rs += d * self.rank(i)
      n -= d
      i += 1
    return rs
def countDict(arr):
  data = {}
  for d in arr:
    if d not in data : data[d] = 0
    data[d] += 1
  return data
def rankDict(data):
  ks = data.keys()
  ks.sort()
  total = 0
  rd = {}
  for k in ks:
    total += data[k]
    rd[k] = total - (data[k] - 1) / 2.0
  return rd

def zipDict(ks, vs):
  data = {}
  for i, v in enumerate(vs) : data[ks[i]] = v
  return data

def orderIter(arr, reverse = False):
  ad = {}
  for i, v in enumerate(arr): 
    if v not in ad : ad[v] = []
    ad[v].append(i)
  #print ad
  a = ad.keys()
  a.sort(reverse = reverse)
  for v in a : 
    for i in ad[v] : yield i

def order(arr, reverse = False):
  return [i for i in orderIter(arr, reverse = reverse)]

def multiHypergeoProb(pcd, cd, show = False):
  lp = - combination_log(sum(cd.values()), sum(pcd.values()))
  if show : print sum(cd.values()), sum(pcd.values()), lp
  for k in pcd :
    lp += combination_log(cd[k], pcd[k])
    if show : print cd[k], pcd[k], lp
  return math.exp(lp)
def multiHypergeoMergeProb(pcd, cd, N, n, lp0 = None, show = False):
  if lp0 is None : lp = - combination_log(N, n)
  else : lp = lp0
  if show : print lp
  #m = max(pcd)
  #for k in cd :
    #if k > m : break

  for k in pcd :
    lp += combination_log(cd[k], pcd[k])
    n -= pcd[k]
    N -= cd[k]
    if show : print cd[k], pcd[k], N, n, lp
  lp += combination_log(N, n)
  return math.exp(lp)

def glmNBTest(x, y):
  import statsmodels.formula.api as smf
  import statsmodels.api as sm
  if max(x) == 0 : return 1.
  data = {}
  data['data'] = list(x) + list(y)
  data['grp'] = [1]*len(x) + [0]*len(y)
  model=smf.glm("data ~ grp", data=data, family=sm.families.NegativeBinomial()).fit()
  p = model.pvalues[1] / 2
  if model.tvalues[1] >= 0 : return p
  else : return 1-p