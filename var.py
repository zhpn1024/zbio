'''
Genomic variant processing
Copyright (c) 2022 Peng Zhang <zhpn1024@163.com>
'''

chr_order = {'1':1, '2':2, '3':3, '4':4, '5':5, '6':6, '7':7, '8':8, '9':9, '10':10, '11':11, '12':12, '13':13, '14':14, '15':15, '16':16,
             '17':17, '18':18, '19':19, '20':20, '21':21, '22':22, 'X':23, 'Y':24, 'M':25}

def cmpchr(c1, c2):
  if c1 not in chr_order: return 1
  if c2 not in chr_order: return -1
  return cmp(chr_order[c1], chr_order[c2])

def cmp3(a, b):
  ''' cmp for python3'''
  return (a > b) - (a < b)

class Variant:
  '''Single variant record'''
  def __init__(self, chr, pos, ref, alt, id = '.'):
    self.chr, self.pos, self.id = chr, int(pos), id
    self.ref = ref
    if type(alt) == str:
      self.alt = alt.split(',')
    else:
      self.alt = list(alt)
    self.alleles = [ref] + self.alt

  def get_allele(self, i):
    return self.alleles[int(i)]

  def __cmp__(self, other):
    return cmpchr(self.chr, other.chr) or cmp(self.pos, other.pos)

  def cmp(self, other):
    return cmpchr(self.chr, other.chr) or cmp3(self.pos, other.pos)
  def __eq__(self, other):
    return self.cmp(other) == 0
  def __lt__(self, other):
    return self.cmp(other) < 0
  def __gt__(self, other):
    return self.cmp(other) > 0

  def __str__(self):
    return '{}\t{}\t{}\t{}\t{}'.format(self.chr, self.pos, self.id, self.ref, ','.join(self.alt))
  def __repr__(self):
    return '{}:{}:{}:{}>{}'.format(self.chr, self.pos, self.id, self.ref, ','.join(self.alt))

class GenoType:
  '''Genotype for one variant'''
  def __init__(self, var, gtstr, phased = False):
    self.var = var
    self.phased = phased
    if '|' in gtstr:
      self.phased = True
      self.gt = gtstr.split('|')
    else:
      self.gt = gtstr.split('/')

  def __str__(self):
    if self.phased: sep = '|'
    else: sep = '/'
    return sep.join(self.gt)

  def __repr__(self):
    return '{} {}'.format(sef.var.__repr__(), self.__str__())

  def __len__(self):
    return len(self.gt)

  def has_unknown(self):
    return '.' in self.gt

  def is_contain(self, a):
    return str(a) in self.gt
  def is_homo(self):
    if self.has_unknown(): return False
    return self.gt[0] == self.gt[1]
  def other_allele(self, a):
    if str(a) == self.gt[0]: return self.gt[1]
    else: return self.gt[0]


class VarList:
  '''Variant list for all variants'''
  def __init__(self):
    self.vars = {}
    self.lst = []

  def add(self, var):
    if var.chr not in self.vars: self.vars[var.chr] = {}
    if var.pos not in self.vars[var.chr]:
      self.vars[var.chr][var.pos] = var
      self.lst.append(var)

  def is_contain(self, chr, pos):
    if chr not in self.vars: return False
    if pos not in self.vars[chr]: return False
    return True

  def get(self, chr, pos):
    if not self.is_contain(chr, pos): return None
    return self.vars[chr][pos]

  def sort(self, reverse = False):
    self.lst.sort(reverse = reverse)


def getGT(s, na, igq = None, iad = None, gqth = 10, dpth = 5, adth = 0.2):
  d = [0] * na
  lst = s.split(':')
  if lst[0].startswith('.'): return d
  lst2 = lst[0].split('/')
  if len(lst2) < 2:
    lst2 = lst[0].split('|')
    if len(lst2) < 2:
      return d
  if igq is not None:
    if len(lst) <= igq: print(s, igq)
    GQ = lst[igq]
    if GQ == '.' or int(GQ) < gqth:
      return d
  if iad is not None: ###
    AD = eval('[' + lst[iad] + ']')
    adsum = sum(AD)
    if adsum < dpth:
      return d
  for i in lst2:
    if i == '.': continue
    i = int(i)
    if iad is not None and AD[i] < adsum * adth:
      continue
    d[i] += 1
  return d

def countGT(lst):
  n = len(lst)
  var = Variant(lst[0], lst[1], lst[3], lst[4])
  vn = len(var.alleles)
  data = [getGT(lst[i], vn) for i in range(9, n)]
  cd = {}
  for d in data:
    d = tuple(d)
    if d not in cd: cd[d] = 0
    cd[d] += 1
  return cd
