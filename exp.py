class exp():
  def __init__(self, sample, exp): #sample and expression list
    self.sample = sample
    self.exp = exp
    self.value = [] #sort index
  def __str__(self):
    return '\t'.join(map(str, self.exp))
  def __repr__(self):
    return self.headerline + "\n" + str(self)
  def __len__(self):
    return len(self.sample)
  def __cmp__(self, other):
    c = 0
    for i in range(len(self.value)):
      c = c or cmp(self.value[i], other.value[i])
      if c != 0: break
    return c
  def headerline(self,sep='\t'):# Header string, fit all bed
    sep=str(sep)
    return sep.join(map(str, self.sample))
    
class trans(exp):
  def __init__(self, tid, sample, exp):
    self.id = tid
    exp.__init__(self, sample, exp)
  def __str__(self):
    return self.id + '\t' + exp.__str__(self)
  def __repr__(self):
    return 'tid\t' + '\t'.join(map(str, self.sample)) + "\n" + str(self)
  @property
  def tid(self):
    return self.id
  def headerline(self,sep='\t'):# Header string, fit all bed
    sep=str(sep)
    return 'tid' + sep + sep.join(map(str, self.sample))
    
class gene(exp):
  def __init__(self, gid, sample, exp):
    self.id = gid
    exp.__init__(self, sample, exp)
    self.trans = []
  def __str__(self):
    if len(self.trans) == 0:
      return self.id + '\t' + exp.__str__(self)
    else:
      s = ''
      for t in self.trans:
        s += self.id + '\t' + str(t)
      return s
  def __repr__(self):
    if len(self.trans) == 0:
      return 'gid\t' + '\t'.join(map(str, self.sample)) + "\n" + str(self)
    else:
      return 'gid\ttid\t' + '\t'.join(map(str, self.sample)) + "\n" + str(self)
  @property
  def gid(self):
    return self.id
  def add_trans(self, t):
    self.trans.append(t)
  def headerline(self,sep='\t'):# Header string, fit all bed
    sep=str(sep)
    if len(self.trans) == 0:
      return 'gid' + sep + sep.join(map(str, self.sample))
    else:
      return 'gid' + sep + 'tid' + sep + sep.join(map(str, self.sample))
  
def subarr(lst, ids = [], head = 0):
  if len(ids) == 0:
    l = len(lst)
    return lst[head:l]
  a = []
  for i in ids:
        a.append(lst[i])
  return a
def gtexp_iter(expfile, gi = 0, ti = 1, sep = '\t', ei = [], sample = []):
  if len(sample) == 0 :
    l = expfile.next()
    lst = l.strip().split(sep)
    sample = subarr(lst, ei, ti+1)
  for l in expfile:
    lst = l.strip().split(sep)
    n = len(lst)
    gid = lst[gi]
    tid = lst[ti]
    exp = map(float, subarr(lst, ei, ti+1))
    t = Trans(tid, sample, exp)
    t.gid = gid
    yield t
def gtexp_load(expfile, gi = 0, ti = 1, sep = '\t', ei = [], sample = []):
  gs = {}
  for t in gtexp_iter(expfile, gi, ti, sep, ei, sample):
    if t.gid not in gs:
      gs[t.gid] = gene(t.gid, t.sample, t.exp)
    gs[t.gid].add_trans(t)
  return gs

def texp_iter(expfile, ti = 0, sep = '\t', ei = [], sample = []):
  if len(sample) == 0 :
    l = expfile.next()
    lst = l.strip().split(sep)
    sample = subarr(lst, ei, ti+1)
  for l in expfile:
    lst = l.strip().split(sep)
    n = len(lst)
    #gid = lst[gi]
    tid = lst[ti]
    exp = map(float, subarr(lst, ei, ti+1))
    t = Trans(tid, sample, exp)
    #t.gid = gid
    yield t

def gexp_iter(expfile, gi = 0, sep = '\t', ei = [], sample = []):
  if len(sample) == 0 :
    l = expfile.next()
    lst = l.strip().split(sep)
    sample = subarr(lst, ei, gi+1)
  for l in expfile:
    lst = l.strip().split(sep)
    n = len(lst)
    #gid = lst[gi]
    gid = lst[gi]
    exp = map(float, subarr(lst, ei, gi+1))
    g = gene(gid, sample, exp)
    #t.gid = gid
    yield g

