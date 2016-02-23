def is_overlap(i1, i2):
  return i1[0] < i2[1] and i1[1] > i2[0]

class interval: # all intervals are supposed to be [start, end) and start should be <= end
  def __init__(self, itvs = [], start = 0, stop = 0, id = ''):
    self.id = id
    self.lst = []
    if itvs == [] and start != stop : self.lst.append([start, stop])
    else :
      for itv in itvs:
        self.lst.append(itv[:])
    self.check()
  def check(self, adj_merge = False) :
    for itv in self.lst[:]:
      if itv[0] > itv[1] : self.lst.remove(itv)
        #temp = itv[0]
        #itv[0] = itv[1]
        #itv[1] = temp
    self.lst.sort()
    #print self.lst
    dl = []
    j = 0
    for i in range(1, len(self.lst)):
      if self.lst[i][0] < self.lst[j][1] or (adj_merge and self.lst[i][0] == self.lst[j][1]) :
        self.lst[j][1] = self.lst[i][1]
        dl.append(i)
      else : j = i
    for i in dl[::-1]: del self.lst[i]
  def __len__(self):
    return len(self.lst)
  def rlen(self):
    l = 0
    for itv in self.lst: l += itv[1] - itv[0]
    return l
  def __repr__(self):
    s = 'interval '+ self.id + ':'
    if len(self.lst) == 0 : s += ' empty!'
    for itv in self.lst:
      s += ' '+str(itv[0])+'-'+str(itv[1])
    return s
  def __getitem__(self, i): #All bed
    return self.lst[i]
  def __add__(self, other): # do not change self
    new = interval(self)
    #print new
    for itv in other:
      new.lst.append(itv[:])
    new.check()
    return new
  def add(self, other):
    return self + other
  def add_itv(self, itv): # input single interval, change the object
    self.lst.append(itv[:])
    self.check()
    return self
  def sub_itv(self, itv):
    lst = []
    for i in self.lst:
      if is_overlap(i, itv) : 
        if i[0] < itv[0] : lst.append([i[0],itv[0]])
        if i[1] > itv[1] : lst.append([itv[1],i[1]])
      else : lst.append(i)
    self.lst = lst
    self.check()
    return self
  def __sub__(self, other): # not optimized! m * n
    new = interval(self)
    for itv in other:
      new.sub_itv(itv[:])
    #new.check()
    return new
  def sub(self, other):
    return self - other
  def ints_itv(self, itv): # intersect 
    lst = []
    for i in self.lst:
      if is_overlap(i, itv): lst.append(max(i[0], itv[0]), min(i[1], itv[1]))
    self.lst = lst
    self.check()
    return self
  def intersect(self, other): # optimized
    lst = []
    j = 0
    for itv in self.lst:
      while j < len(other) and other[j][1] < itv[0] : j += 1
      j1 = j
      while j1 < len(other) and is_overlap(itv, other[j1]) : 
        lst.append([max(itv[0], other[j1][0]), min(itv[1], other[j1][1])])
        j1 += 1
    new = interval(lst)
    return new
    #return (self + other) - (self - other) - (other - self) # hehe
  def is_inside(self, p, left = True, right = False) : # is p inside interval
    for itv in self.lst:
      if itv[0] < p < itv[1] : return True
      if left and itv[0] == p : return True
      if right and itv[1] == p : return True
    return False
  def num_iter(self, start = None, step = 1):
    if self.rlen() <= 0 : return
    if start is None : start = self.start
    i = start
    while i < self.stop : 
      if self.is_inside(i) : yield i
      i += step
    #for i in range(start, self.stop, step):
      #if self.is_inside(i) : yield i
  @property
  def start(self):
    return self.lst[0][0]
  @property
  def stop(self):
    return self.lst[-1][1]
  def is_empty(self):
    return len(self.lst) == 0
