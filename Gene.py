Gid = ['ID','GENEID','GI','GID','ENTREZ']
Sym = ['SYM','SYMBOL','NAME']
Ali = ['ALI','ALIAS','ALIASES','SYNONYM','SYNONYMS']
Ens = ['ENS','ENSEMBL','ENSG']
Full = ['FULL','FULLNAME']
class Gene():
  def __init__(self, l, sep = '\t', idx = [1,2,4,5,8,0]):
    lst = l.strip().split(sep)
    self.gid = lst[idx[0]]
    self.sym = lst[idx[1]]
    self.ali = lst[idx[2]].split('|')
    self.attr = {}
    alst = lst[idx[3]].split('|')
    for a in alst:
      al = a.split(':')
      if len(al) > 1:
        self.attr[al[0]] = a[len(al[0])+1:]
    #print lst[idx[3]],al,self.attr
    self.full = lst[idx[4]]
    self.taxid = lst[idx[5]]
  def __str__(self):
    return self.sym
  def __repr__(self):
    return "Gene object: %s %s %s %s %s" % (self.gid, self.sym, self.ens, self.ali, self.full)
  def __getattr__(self, name):
    nu = name.upper()
    if nu in Gid:
      return self.gid
    elif nu in Sym:
      return self.sym
    elif nu in Ali:
      return self.ali
    elif nu in Ens:
      try: return self.attr['Ensembl']
      except: return ''
    elif nu in Full:
      return self.full
    elif name in self.attr:
      return self.attr[name]
    else:
      raise AttributeError, name

class GeneDict():
  identifier = ('GeneID', 'Symbol', 'Ensembl')
  def __init__(self, ginfofile, sep = '\t'): # gene_info file for one species from NCBI
    self.gid = {}
    self.sym = {}
    self.ali = {}
    self.ens = {}
    self.symUp = {} #Upper case
    self.aliUp = {}
    lastsp = ''
    for l in ginfofile:
      if l[0] == '#' : continue
      g = Gene(l, sep)
      if lastsp != g.taxid and lastsp != '':
        print "Warning: There may be more than one species in gene_info file!"
      lastsp = g.taxid
      self.gid[g.gid] = g
      self.sym[g.sym] = g
      self.symUp[g.sym.upper()] = g
      for a in g.ali:
        if a not in ['','-']: 
          self.ali[a] = g
          self.aliUp[a.upper()] = g
      if g.ens not in ['','-']: self.ens[g.ens] = g
  
  def __call__(self, *name, **attr):
    namelist = []
    gids = {} # Most supported names
    for n in name:
      if type(n) in [list, tuple]:
        for nn in n : namelist.append(nn)
      elif type(n) == dict:
        for nn in n : attr[nn] = n[nn]
      else: namelist.append(n)
    for a in attr:
      au = a.upper()
      if au in Gid:
        g = self.byGid(attr[a])
        if g != None: return g
      elif au in Ens:
        g = self.byEns(attr[a])
        if g != None: return g
      else: 
        namelist.append(attr[a])
    for n in namelist:
        g = self.bySym(n)
        if g == None: continue
        if g.gid not in gids: gids[g.gid] = 0
        gids[g.gid] += 1
    if len(gids) == 0 : return None
    gid = sorted(gids.items(), key=lambda d:d[1])[-1][0]
    return self.gid[gid]

  def bySym(self, sym):
    sym = str(sym)
    if sym in self.sym: return self.sym[sym]
    if sym in self.ali: return self.ali[sym]
    su = sym.upper()
    if su in self.symUp: return self.symUp[su]
    if su in self.aliUp: return self.aliUp[su]
    return None
  def byGid(self, gid):
    gid = str(gid)
    if gid in self.gid: return self.gid[gid]
    return None
  def byEns(self, ens):
    ens = str(ens)
    if ens in self.ens: return self.ens[ens]
    return None
  def isAli(self, gid):
    return gid in self.alias