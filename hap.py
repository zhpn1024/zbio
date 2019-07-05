'''
Genomic haplotype processing
Copyright (c) 2017 Peng Zhang <zhpn1024@163.com>
'''

class SNP:
  def __init__(self, chr, pos, ref, alt, a1 = '.', a2 = '.', score = 1):
    self.chr, self.pos = chr, pos
    self.ref, self.alt = ref, alt
    self.a1, self.a2 = a1, a2
    self.score = score
    self.id = '{}_{}'.format(chr, pos)

  def __gt__(self, other):
    return self.pos > other.pos
  def __lt__(self, other):
    return self.pos < other.pos
  def __cmp__(self, other):
    return cmp(self, other)

  def copy(self, other, reverse = False):
    if not reverse: self.a1, self.a2 = other.a1, other.a2
    else: self.a1, self.a2 = other.a2, other.a1

  def is_homo(self):
    return self.a1 == self.a2

class HapBlock:
  def __init__(self, chr = None, start = None, stop = None):
    self.chr = chr
    self.start = start
    self.stop = stop
    self.snps = []
    self.allsnps = {}

  def add(self, snp, check = False):
    if snp.a1 not in ('0', '1') or snp.a2 not in ('0', '1'):
      #print('Unknown phasing: {} {} {} {}'.format(snp.chr, snp.pos, snp.a1, snp.a2))
      return
    if self.chr is None: self.chr = snp.chr
    elif snp.chr != self.chr:
      print('chr not match: {} {}'.format(snp.chr, self.chr))
      return
    self.snps.append(snp)
    self.allsnps[snp.id] = snp
    if self.start is None: self.start = snp.pos
    self.stop = snp.pos
    if check: self.check()

  def check(self):
    if len(self.snps) == 0: return
    self.snps.sort()
    self.start = self.snps[0].pos #: self.start = snp.pos
    self.stop = self.snps[-1].pos #: self.stop = snp.pos

  def clear(self):
    self.chr = None
    self.start = None
    self.stop = None
    self.snps = []
    self.allsnps = {}

  def is_contain(self, snp, silent = False):
    if snp.id not in self.allsnps: return False
    ss = self.allsnps[snp.id]
    if snp.ref != ss.ref or snp.alt != ss.alt:
      if not silent: print('SNP not match: {} {}->{} {}->{}'.format(snp.id, snp.ref, snp.alt, ss.ref, ss.alt))
      return False
    return True

  def match(self, other, rate = True, skip_homo = True, count_um = False, silent = False):
    s, n = 0.0, 0
    um = 0
    for snp in other.snps:
      #if skip_homo and snp.is_homo(): continue
      if not self.is_contain(snp, silent = silent):
        if count_um:
          if snp.id in self.allsnps: um += 1
        continue
      if skip_homo and snp.is_homo(): continue
      ss = self.allsnps[snp.id]
      if skip_homo and ss.is_homo(): continue
      ms = snp.score * ss.score
      n += ms
      if snp.a1 == ss.a1: s += ms / 2
      else: s -= ms / 2
      if snp.a2 == ss.a2: s += ms / 2 ##
      else: s -= ms / 2
    if n <= 1: return None
    if rate: return s / n
    else: return s, n, um

  def copy(self, other, match = 1, use_score = True, skip_homo = True):
    if match < 0:
      reverse = True
      m2 = - match # (0.5 - match) / 0.5
    else:
      reverse = False
      m2 = match # (match - 0.5) / 0.5
    for snp in other.snps:
      if skip_homo and snp.is_homo(): continue
      if not self.is_contain(snp): continue
      ss = self.allsnps[snp.id]
      if skip_homo and ss.is_homo(): continue
      if use_score:
        ms = snp.score - m2 * ss.score ##
        if ms < 0: continue
      ss.copy(snp, reverse = reverse)

def rmchr(s):
  if s.startswith('chr'): return s[3:]
  else: return s

def hapcut2_block_iter(hcfile, sep = '\t'):
  #hb = HapBlock()
  for l in hcfile:
    lst = l.rstrip().split(sep)
    if lst[0].startswith('BLOCK'):
      #if len(hb.snps) > 0: yield hb
      hb = HapBlock()
      continue
    elif lst[0].startswith('*'):
      if len(hb.snps) > 0:
        hb.check()
        #print(hb.chr, hb.start, hb.stop)
        yield hb
      continue
    snp = SNP(rmchr(lst[3]), int(lst[4]), lst[5], lst[6], lst[1], lst[2], float(lst[10])/100)
    hb.add(snp)
  if len(hb.snps) > 0:
    hb.check()
    yield hb

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

class Hap:
  '''Haplotype fragment'''
  def __init__(self, chr, varlist):
    self.chr = chr
    self.varlist = varlist
    self.data = {}
    self.qual = {}

  def __repr__(self):
    s = ['Hap object {}'.format(self.chr)]
    ps = sorted(self.data)
    for p in ps:
      s.append('{}\t{}\t{}'.format(p, self.data[p], self.qual[p]))
    return '\n'.join(s)

  def add(self, var, a, qual = 1, replace = False):
    if var.chr != self.chr:
      print('Hap.add: chr not match {} {}'.format(self.chr, var.chr))
      return
    if var.pos in self.data:
      if not replace:
        print('Hap.add: duplicate var {} {} {}'.format(var.chr, var.pos, a))
        return
    else:
      self.varlist.add(var)
    self.data[var.pos] = a
    self.qual[var.pos] = qual

  def match(self, other, rate = False):
    m, n = 0, 0
    for p in other.data:
      if p not in self.data: continue
      if self.data[p] == '.' or other.data[p] == '.': continue
      n += 1
      if self.data[p] == other.data[p]:
        m += 1
    if rate:
      if n == 0: return 1
      else: return float(m) / n
    else: return m, n

  def merge(self, other, replace = False):
    for p in other.data:
      if p in self.data and not replace: continue
      self.data[p] = other.data[p]
      self.qual[p] = other.qual[p]
      #self.varlist.append(other.vars[p]) ## share varlist
    #self.varlist.sort()

class SampleGTsChr:
  '''Genotypes in one sample chr'''
  def __init__(self, chr, varlist, sid = ''):
    self.sid = sid
    self.chr = chr
    self.varlist = varlist
    self.gts = {}

  def add(self, gt):
    if gt.var.chr != self.chr:
      print('SampleGTsChr.add: chr not match {} {}'.format(gt.var.chr, self.chr))
      return
    self.varlist.add(gt.var)
    self.gts[gt.var.pos] = gt

  def is_compatible(self, hap, mis = 0):
    if self.chr != hap.chr:
      print('SampleGTsChr.is_compatible: chr not match {} {}'.format(hap.chr, self.chr))
      return True
    m = 0
    for p, a in hap.data.items():
      if self.gts[p].is_contain(a): continue
      m += hap.qual[p]
      if m > mis: return False
    return True

  def other_hap(self, hap):
    oh = Hap(hap.chr, hap.varlist)
    for p, a in hap.data.items():
      oa = self.gts[p].other_allele(a)
      oh.add(hap.varlist.get(hap.chr, p), oa, hap.qual[p])
    return oh
