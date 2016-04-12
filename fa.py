from os.path import isfile

def fa_iter(file):
  id = sq = ""
  for l in file:
    l = l.strip()
    if l == '' : continue
    if l[0] == '>':
      if id != '':
        yield (id,sq)
      sq = ''
      id = l[1:]
    else:
      sq += l.replace('U','T').replace('u','t')
  yield (id,sq)

class faidx:
  def __init__(self, fid, length = 0, pos = 0, ls = 0, ll = 0):
    self.id = fid
    self.length = length
    self.pos = pos
    self.ls = ls
    self.ll = ll
    self.is_endls = False
    self.is_endll = False
  def __str__(self):
    return "%s\t%d\t%d\t%d\t%d" % (self.id, self.length, self.pos, self.ls, self.ll)
  
  def updatels(self, ls):
    if ls <= 0 : raise Exception("Empty sequence!")
    if self.ls == 0 : self.ls = ls
    elif self.is_endls : raise Exception("Sequence has ended: %d" % (ls))
    elif ls > self.ls :
      raise Exception("Length error: %d -> %d" % (self.ls, ls))
    elif ls < self.ls : self.is_endls = True
  def updatell(self, ll):
    #print ll
    if ll <= 0 : raise Exception("Empty line!")
    if self.ll == 0 : self.ll = ll
    elif self.is_endll : raise Exception("Sequence has ended: %d" % (ll))
    elif ll > self.ll :
      raise Exception("Length error: %d -> %d" % (self.ll, ll))
    elif ll < self.ll : self.is_endll = True
      #print 'End: ', self.ll, ll

class fa:
  def __init__(self, fapath):
    self.file = open(fapath, 'r')
    idxpath = fapath + '.fai'
    if isfile(idxpath) : self.load_idx(idxpath)
    else : self.make_idx(idxpath)
  def make_idx(self, idxpath):
    self.file.seek(0)
    self.idxarr = []
    pos, flen = 0, 0
    for l in self.file:
      pos += len(l)
      if l[0] == '>' : 
        if flen > 0 : self.idxarr[-1].length = flen
        fid = l.strip().split()[0][1:]
        #print fid
        self.idxarr.append(faidx(fid, pos = pos))
        flen = 0
      else : 
        ll, ls = len(l), len(l.strip())
        flen += ls
        self.idxarr[-1].updatels(ls)
        self.idxarr[-1].updatell(ll)
        #if ll > self.idxarr[-1][4] : self.idxarr[-1][4] = ll
    self.idxarr[-1].length = flen
    outfile = open(idxpath, 'w')
    self.idx = {}
    for i, idx in enumerate(self.idxarr):
      print >>outfile, idx
      self.idx[idx.id] = idx
    outfile.close()
  def load_idx(self, idxpath):
    self.idx = {}
    for l in open(idxpath):
      lst = l.strip().split('\t')
      if len(lst) < 5 : continue
      idx = faidx(lst[0])
      idx.length, idx.pos, idx.ls, idx.ll = map(int, lst[1:5])
      self.idx[idx.id] = idx

  def seekpos(self, fid, p):
    idx = self.idx[fid]
    pos = idx.pos # chr start
    pos += int(p / idx.ls) * idx.ll # line start
    pos += p % idx.ls # base start
    return pos
  def fetch(self, fid, start = 0, stop = -1, length = -1):
    if start < 0 : start = start % self.idx[fid].length
    if length < 0 : 
      if stop < 0 : stop = stop % self.idx[fid].length
      length = stop - start
    if length <= 0 : return ''
    self.file.seek(self.seekpos(fid, start))
    seq = ''
    idx = self.idx[fid]
    for l in self.file:
      l = l.strip()
      if length > len(l) :
        seq += l
        length -= len(l)
      else : 
        seq += l[0:length]
        break
    return seq
  def __iter__(self):
    for fid in self.idx:
      yield fid
    #return self.idx.itervalues()
    
