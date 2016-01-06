class table():
  def __init__(self, lst = [], sep = '\t', byrow = False):
    self.matrix = []
    self.colnames = []
    self.rownames = []
    self.ncol = 0
    self.nrow = 0
    self.is_neat = True
    if len(lst) > 0 :
      for s in lst:
        self.matrix.append(s.split(sep))
      self.check()
    if byrow : self.transpose()
  def load(self, infile, header = True, rowname = False, sep = '\t', n = -1, skip = 0, innerskip = []):
    i = 0
    while i < skip: 
      l = infile.next()
      i += 1
    if header :
      l = infile.next()
      i += 1
      lst = l.strip('\n').split(sep)
      self.colnames = lst
    for l in infile:
      i += 1
      if n >= 0 and i > n : break
      if i in innerskip : continue
      lst = l.strip('\n').split(sep)
      if rowname : 
        self.rownames.append(lst[0])
        self.matrix.append(lst[1:])
      else : self.matrix.append(lst)
    #print self.matrix
    self.check()
  def check(self, silent = True):
    self.nrow = len(self.matrix)
    if not silent : print 'nrows =', self.nrow,
    cols = [len(lst) for lst in self.matrix]
    self.ncol = max(cols)
    mincol = min(cols)
    if mincol != self.ncol :
      self.is_neat = False
      if not silent : print ', mincols =', mincol, ', maxcols =', self.ncol, ', NOT neat!'
    else :
      self.is_neat = True
      if not silent : print ', ncols =', self.ncol, ', neat.'
  def headerline(self, sep = '\t'):
    return sep.join([str(x) for x in self.colnames])
  def write(self, outfile, header = True, rowname = False, sep = '\t'):
    if header : outfile.write(self.headerline(sep) + '\n')
    for i in range(self.nrow):
      s = ''
      if rowname :
        try : s += str(self.rownames[i]) + sep
        except : s += str(i) + sep
      s += sep.join([str(x) for x in self.matrix[i]])
      outfile.write(s + '\n')
  def rowidx(self, i):
    if type(i) == int : return i
    else : return self.rownames.index(i)
  def colidx(self, j):
    if type(j) == int : return j
    else : return self.colnames.index(j)
  def row(self, i): #the ith row
    if type(i) == int : return self.matrix[i]
    else : return self.matrix[self.rownames.index(i)]
  def col(self, j): # the ith column
    if type(j) == int : idx = j
    else : idx = self.colnames.index(j)
    lst = [None] * self.nrow
    for i in range(self.nrow):
      try : lst[i] = self.matrix[i][idx]
      except : pass
    return lst
  def transpose(self):
    pass
  def subtable(self, rows = [], cols = []):
    subt = table()
    if len(rows) == 0 : rows = range(self.nrow)
    if len(cols) == 0 : cols = range(self.ncol)
    for i in rows:
      #print self.matrix[self.rowidx(i)]
      lst = [self.matrix[self.rowidx(i)][self.colidx(j)] for j in cols]
      subt.matrix.append(lst)
    if len(self.colnames) > 0 : subt.colnames = [self.colnames[self.colidx(j)] for j in cols]
    if len(self.rownames) > 0 : subt.rownames = [self.rownames[self.rowidx(j)] for j in rows]
    subt.check()
    return subt