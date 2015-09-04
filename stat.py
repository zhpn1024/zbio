import math
from scipy.stats import nbinom, chisquare
logarr = [None]
def logarrExt(n, logarr = logarr):
  l = len(logarr)
  if l < n + 1 : 
    logarr += [None] * (n + 1 - l)
    for i in range(l, n + 1):
      logarr[i] = math.log(i)
  return logarr
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
  logarrExt(n, logarr = logarr)
  #if len(logarr) >= n + 1 : arr = False
  #else : arr = True
    #logarr = [None] * (n + 1)
    #for i in range(1, n + 1):
      #logarr[i] = math.log10(i)
  for i in range(k):
    lpr += logarr[n - i] - logarr[k - i]
    #else : lpr += math.log10(n - i) - math.log10(k - i)
  return lpr
def binom_test(k, n, p = 0.5, alt = "g", log = True, logarr = logarr, show=False): # No two sided yet!
  if not log : return binomTest0(k, n, p, alt, show) # if log, p are calculated with log10 values
  #logarr = [None] * (n + 1)
  logarrExt(n, logarr = logarr)
  lpk = binom_log(k, n, p, logarr)
  if show : print lpk
  if lpk is None : 
    if alt[0] == 'g' and p >= 1: return 1
    if alt[0] != 'g' and p <= 0: return 1
    return 0
  elif lpk == 0 : return 1
  #if k == 0 and alt[0] == 'g' : return 1
  #if k == n and alt[0] != 'g' : return 1
  pv = math.exp(lpk) ### log10 reverse
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
  def pvalue(self, k = 0):
    return 1 - nbinom.cdf(k-1, self.r, self.p)
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
  def estimate(self, data, max_iter = 1e4, nlike = 10):
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