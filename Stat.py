import math
from scipy.stats import nbinom
logarr = [None]
def logarrExt(n, logarr = logarr):
  l = len(logarr)
  if l < n + 1 : 
    logarr += [None] * (n + 1 - l)
    for i in range(len(logarr), n + 1):
      logarr[i] = math.log(i)
  return logarr
def dataCount(data):
  total, cnt = 0, 0
  for k in data:
    total += k * data[k]
    cnt += data[k]
  return total, cnt
def loadData(arr):
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
def binomLog(k, n, p = 0.5, logarr = logarr, show = False): #log probability value
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
def binomTest(k, n, p = 0.5, alt = "g", log = True, logarr = logarr, show=False): # No two sided yet!
  if not log : return binomTest0(k, n, p, alt, show) # if log, p are calculated with log10 values
  #logarr = [None] * (n + 1)
  logarrExt(n, logarr = logarr)
  lpk = binomLog(k, n, p, logarr)
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