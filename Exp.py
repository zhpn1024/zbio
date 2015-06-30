class Exp():
  def __init__(self, sample, exp): #sample and expression list
    self.sample = sample
    self.exp = exp
  def __str__(self):
    return '\t'.join(map(str, self.exp))
  def __repr__(self):
    return '\t'.join(map(str, self.sample)) + "\n" + str(self)
  def __len__(self):
    return len(self.sample)
    
class Trans(Exp):
  def __init__(self, tid, sample, exp):
    self.id = tid
    Exp.__init__(self, sample, exp)
  def __str__(self):
    return self.id + '\t' + Exp.__str__(self)
  def __repr__(self):
    return 'tid\t' + '\t'.join(map(str, self.sample)) + "\n" + str(self)
  @property
  def tid(self):
    return self.id
    
class Gene(Exp):
  def __init__(self, gid, sample, exp):
    self.id = gid
    Exp.__init__(self, sample, exp)
    self.trans = []
  def __str__(self):
    if len(self.trans) == 0:
      return self.id + '\t' + Exp.__str__(self)
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
  def addTrans(self, t):
    self.trans.append(t)
  
def subarr(lst, ids = [], head = 0):
  if len(ids) == 0:
    l = len(lst)
    return lst[head:l]
  a = []
  for i in ids:
        a.append(lst[i])
  return a
def gtExpIter(expfile, gi = 0, ti = 1, sep = '\t', ei = [], sample = []):
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
def gtExpLoad(expfile, gi = 0, ti = 1, sep = '\t', ei = [], sample = []):
  gs = {}
  for t in gtExpIter(expfile, gi, ti, sep, ei, sample):
    if t.gid not in gs:
      gs[t.gid] = Gene(t.gid, t.sample, t.exp)
    gs[t.gid].addTrans(t)
  return gs

def tExpIter(expfile, ti = 0, sep = '\t', ei = [], sample = []):
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

def gExpIter(expfile, gi = 0, sep = '\t', ei = [], sample = []):
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
    g = Gene(gid, sample, exp)
    #t.gid = gid
    yield g

class libGid():
  def __init__(self, species = 'human', sep = '\t'):
    self.t2g = {}
    self.alias = {}
    if species in ('human', 'hs', 'hg19', 'hg38'):
      t2gfname = 'trans2gene_hg19.txt'
      alifname = 'ali_human.txt'
    t2gfile = open(t2gfname, 'r')
    for l in t2gfile:
      lst = l.strip().split(sep)
      self.t2g[lst[0]] = lst[1]
    alifile = open(alifname, 'r')
    for l in alifile:
      lst = l.strip().split(sep)
      self.alias[lst[0]] = lst[1]
  def checkGid(self, t):
    return t.gid == self.t2g[t.id]
  def isAli(self, gid):
    return gid in self.alias