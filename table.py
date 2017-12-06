'''
Table processing
Copyright (c) 2016 Peng Zhang <zhpn1024@163.com>
'''

indexTypes = (int, str, float)
def intersect(x, y):
  sy = set(y)
  return [i for i in x if i in sy]
class Table():
  def __init__(self, lst = [], sep = '\t', byrow = False, default = None):
    self.matrix = []
    self.colnames = []
    self.rownames = []
    self.ncol = 0
    self.nrow = 0
    self.is_neat = True
    self.col_indexed = False
    self.row_indexed = False
    self.ridx = {}
    self.cidx = {}
    self.default = default
    if len(lst) > 0 :
      for s in lst:
        self.matrix.append(s.split(sep))
      self.check()
    #if byrow : self = self.transpose()
  def __call__(self, i = None, j = None):
    if type(i) in indexTypes and type(j) in indexTypes : return self.matrix[self.rowidx(i)][self.colidx(j)]
    return self.subtable(i, j)
  def load(self, infile, header = True, rowname = False, sep = '\t', n = -1, skip = 0, innerskip = []):
    i = 0
    while i < skip: 
      l = next(infile)
      i += 1
    if header :
      l = next(infile)
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
    if not silent : print ('nrows = {}'.format(self.nrow))
    cols = [len(lst) for lst in self.matrix]
    self.ncol = max(cols)
    mincol = min(cols)
    if mincol != self.ncol :
      self.is_neat = False
      if not silent : print (', mincols =', mincol, ', maxcols =', self.ncol, ', NOT neat!')
    else :
      self.is_neat = True
      if not silent : print (', ncols =', self.ncol, ', neat.')
  def __repr__(self):
    s = "nrows = %d, ncols = %d" % (self.nrow, self.ncol)
    if self.is_neat: s += ', neat!\n'
    else: s += ', NOT neat.\n'
    s += self.string(nrow = 5, ncol = 5, rowname = True)
    return s
  def __str__(self): pass
  def string(self, nrow = None, ncol = None, header = True, rowname = False, sep = '\t'):
    if nrow is None : nrow = self.nrow
    else : nrow = min(nrow, self.nrow)
    if ncol is None : ncol = self.ncol
    else : ncol = min(ncol, self.ncol)
    s = []
    if header : 
      if len(self.colnames) == 0 :colnames = list(range(ncol))
      else : colnames = self.colnames
      if rowname : s.append(sep)
      s.append(sep.join([str(x) for x in colnames[0:ncol]]) + '\n')
    for i in range(nrow):
      if rowname : 
        try : s.append(str(self.rownames[i]) + sep)
        except : s.append(str(i) + sep)
      s.append(sep.join([str(x) for x in self.matrix[i][0:ncol]]) + '\n')
    return ''.join(s)
  def headerline(self, sep = '\t', rowname = False):
    s = sep.join([str(x) for x in self.colnames])
    if rowname: return 'Rownames\t' + s
    else: return s # sep.join([str(x) for x in self.colnames])
  def write(self, outfile, header = True, rowname = False, sep = '\t'):
    if header : outfile.write(self.headerline(sep, rowname) + '\n')
    for i in range(self.nrow):
      s = ''
      if rowname :
        try : s += str(self.rownames[i]) + sep
        except : s += str(i) + sep
      s += sep.join([str(x) for x in self.matrix[i]])
      outfile.write(s + '\n')
  def indexrow(self):
    #if len(self.rownames) == 0 : return
    #self.ridx = {}
    for i, key in enumerate(self.rownames):
      self.ridx[key] = i
    i = len(self.rownames)
    if i < self.nrow:
      for i2 in range(i, self.nrow):
        self.ridx[str(i2)] = i2
        self.rownames.append(str(i2))
    self.row_indexed = True
  def indexcol(self):
    #if len(self.colnames) == 0 : return
    #self.cidx = {}
    for i, key in enumerate(self.colnames):
      self.cidx[key] = i
    i = len(self.colnames)
    if i < self.ncol:
      for i2 in range(i, ncol):
        self.cidx[str(i2)] = i2
        self.colnames.append(str(i2))
    self.col_indexed = True
  def rowidx(self, i):
    if type(i) == str : 
      if self.row_indexed: return self.ridx[i]
      else : return self.rownames.index(i)
    else : return int(i)
  def colidx(self, j):
    if type(j) == str : 
      if self.col_indexed: return self.cidx[j]
      else : return self.colnames.index(j)
    else : return int(j)
  def row(self, i): #the ith row
    return self.matrix[self.rowidx(i)]
    #if type(i) == int : return self.matrix[i]
    #else : return self.matrix[self.rownames.index(i)]
  def col(self, j): # the ith column
    idx = self.colidx(j)
    #if type(j) == int : idx = j
    #else : idx = self.colnames.index(j)
    lst = [self.default] * self.nrow
    for i in range(self.nrow):
      try : lst[i] = self.matrix[i][idx]
      except : pass
    return lst
  def transpose(self):
    self.check()
    t = table()
    for j in range(self.ncol) :
      t.matrix.append(self.col(j))
    if len(self.colnames) > 0 : t.rownames = self.colnames[:]
    if len(self.rownames) > 0 : t.colnames = self.rownames[:]
    t.check()
    return t
  def subtable(self, rows = None, cols = None):
    subt = table()
    if rows is None : rows = range(self.nrow)
    if cols is None : cols = range(self.ncol)
    if type(rows) in indexTypes : rows = [rows]
    if type(cols) in indexTypes : cols = [cols]
    for i in rows:
      #print self.matrix[self.rowidx(i)]
      lst = [self.matrix[self.rowidx(i)][self.colidx(j)] for j in cols]
      subt.matrix.append(lst)
    if len(self.colnames) > 0 : subt.colnames = [self.colnames[self.colidx(j)] for j in cols]
    if len(self.rownames) > 0 : subt.rownames = [self.rownames[self.rowidx(j)] for j in rows]
    subt.check()
    return subt
  def put(self, row, col, value):
    if type(row) is int:
      if row >= 0:
        i = row
        if self.nrow < i + 1 :
          for i2 in range(self.nrow, i+1):
            self.matrix.append([])
          self.nrow = i + 1
          #radd = True
      else:
        i = self.nrow + i
        if i < 0: 
          print('Invalid negative row: {}'.format(row))
          return
    else:
      if not self.row_indexed:
        self.indexrow()
      if row in self.ridx:
        i = self.ridx[row]
      else:
        self.rownames.append(row)
        i = self.nrow
        self.nrow += 1
        self.ridx[row] = i
        self.matrix.append([])
        #radd = True
    
    if type(col) is int:
      if col >= 0:
        j = col
        if self.ncol < j + 1 :
          self.ncol = j + 1
          #radd = True
      else:
        j = self.ncol + j
        if i < 0:
          print('Invalid negative col: {}'.format(col))
          return
    else:
      if not self.col_indexed:
        self.indexcol()
      if col in self.cidx:
        j = self.cidx[col]
      else:
        self.colnames.append(col)
        j = self.ncol
        self.ncol += 1
        self.cidx[col] = j

    l = len(self.matrix[i])
    if l < j+1:
      self.matrix[i] += [self.default] * (j+1-l)
    self.matrix[i][j] = value

    return i, j, value


