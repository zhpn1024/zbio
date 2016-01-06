class exp(): #values for one gene/trans/probe
  def __init__(self, id, sample, data, anno = ''): #sample and expression list
    self.id = id
    self.sample = sample
    self.data = data
    self.anno = anno
    self.value = [] #sort index
  def __str__(self, showanno = False, sep = '\t'):
    s = self.id + sep
    if showanno : s += self.anno + sep
    return s + sep.join(map(str, self.data))
  def string(self, showanno = False, sep = '\t'):
    return self.__str__(showanno, sep)
  def __repr__(self):
    return self.headerline() + "\n" + str(self)
  def __len__(self):
    return len(self.sample)
  def __cmp__(self, other):
    c = 0
    for i in range(len(self.value)):
      c = c or cmp(self.value[i], other.value[i])
      if c != 0: break
    return c
  def headerline(self, sep='\t'):# Header string, fit all bed
    sep=str(sep)
    return 'id' + sep + sep.join(map(str, self.sample))
    
class trans(exp):
  def __init__(self, tid, sample, data):
    self.id = tid
    exp.__init__(self, tid, sample, data)
  #def __str__(self):
    #return self.id + '\t' + exp.__str__(self)
  #def __repr__(self):
    #return 'tid\t' + '\t'.join(map(str, self.sample)) + "\n" + str(self)
  @property
  def tid(self):
    return self.id
  def headerline(self,sep='\t'):# Header string, fit all bed
    sep=str(sep)
    return 'tid' + sep + sep.join(map(str, self.sample))
    
class gene(exp):
  def __init__(self, gid, sample, data):
    self.id = gid
    exp.__init__(self, gid, sample, data)
    self.trans = []
  #def __str__(self): #gene only or gene & trans
  #  #if len(self.trans) == 0:
  #    #return self.id + '\t' + exp.__str__(self)
  #  #else:
  #    #s = ''
  #    #for t in self.trans:
  #      #s += self.id + '\t' + str(t)
  #    #return s
  #def __repr__(self):
  #  if len(self.trans) == 0:
  #    return 'gid\t' + '\t'.join(map(str, self.sample)) + "\n" + str(self)
  #  else:
  #    return 'gid\ttid\t' + '\t'.join(map(str, self.sample)) + "\n" + str(self)
  @property
  def gid(self):
    return self.id
  def add_trans(self, t):
    self.trans.append(t)
  def headerline(self,sep='\t'):# not finished!
    sep=str(sep)
    if len(self.trans) == 0:
      return 'gid' + sep + sep.join(map(str, self.sample))
    else:
      return 'gid' + sep + 'tid' + sep + sep.join(map(str, self.sample))
  
def subarr(lst, ids = [], head = 0): # selected items in the list
  if len(ids) == 0:
    #l = len(lst)
    return lst[head:]
  a = []
  for i in ids:
        a.append(lst[i])
  return a
def gtexp_iter(expfile, gi = 0, ti = 1, sep = '\t', ei = [], sample = []):
  if len(sample) == 0 :
    l = expfile.next()
    lst = l.strip('\n').split(sep)
    sample = subarr(lst, ei, ti+1)
  for l in expfile:
    lst = l.strip('\n').split(sep)
    n = len(lst)
    gid = lst[gi]
    tid = lst[ti]
    data = map(float, subarr(lst, ei, ti+1))
    t = trans(tid, sample, data)
    t.gid = gid
    yield t
def gtexp_load(expfile, gi = 0, ti = 1, sep = '\t', ei = [], sample = []):
  gs = {}
  for t in gtexp_iter(expfile, gi, ti, sep, ei, sample):
    if t.gid not in gs:
      gs[t.gid] = gene(t.gid, t.sample, t.data)
    gs[t.gid].add_trans(t)
  return gs

def exp_iter(expfile, ii = 0, dsi = -1, annoi = -1, sep = '\t', header = True, ei = [], sample = [], skip = 0, innerskip = []):
  if dsi < 0 : dsi = max(ii, annoi) + 1 # supposed data start id
  i = 0
  while i < skip: 
    l = expfile.next()
    i += 1
  if header :
    l = expfile.next()
    i += 1
    lst = l.strip('\n').split(sep)
    sample = subarr(lst, ei, dsi)
  for l in expfile:
    i += 1
    if i in innerskip : continue
    lst = l.strip('\n').split(sep)
    n = len(lst)
    id = lst[ii]
    anno = ''
    if annoi >= 0 : anno = lst[annoi]
    data = map(float, subarr(lst, ei, dsi))
    e = exp(id, sample, data, anno)
    #t.gid = gid
    yield e

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
    data = map(float, subarr(lst, ei, gi+1))
    g = gene(gid, sample, data)
    #t.gid = gid
    yield g

class profile():
  def __init__(self):
    self.exps = {}
    self.ids = []
  def add_exp(self, e):
    self.exps[e.id] = exp
    self.ids.append(e.id)
  def __len__(self):
    return len(self.exps)
  def BHcorrection(self, pid = -1, total = -1):
    lst = self.exps.values()
    n = len(lst)
    if total < 0: total = n
    for e in lst:
      if len(e.value) < 1: e.value.append(1)
      e.value[0] = e.data[pid]
    lst.sort()
    qc = 1
    for i in range(n-1, -1, -1):
      q = float(lst[i].value[0]) * total / (i+1)
      if q > qc : q = qc
      lst[i].q = q
      qc = q
    return lst
  def write(outfile, header = True, showanno = False, sep = '\t'):
    for id in self.ids:
     outfile.write(self.exps[id].string(showanno, sep) + '\n')
    