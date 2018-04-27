'''
Gene symbol annotation processing
Copyright (c) 2016 Peng Zhang <zhpn1024@163.com>
'''

import sys
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
      raise AttributeError(name)

class geneDict():
  identifier = ('GeneID', 'Symbol', 'Ensembl')
  def __init__(self, ginfofile, sep = '\t'): # gene_info file for one species from NCBI
    self.gid = {}
    self.sym = {}
    self.ali = {}
    self.ens = {}
    self.sym_up = {} #Upper case
    self.ali_up = {}
    lastsp = ''
    neg = {}
    for l in ginfofile:
      if l[0] == '#' : continue
      g = Gene(l, sep)
      if lastsp != g.taxid and lastsp != '':
        if g.taxid not in neg:
          sys.stderr.write("Warning: There may be more than one species in gene_info file! Neglected: "+g.taxid+"\n")
          neg[g.taxid] = 1
        continue
      lastsp = g.taxid
      self.gid[g.gid] = g
      self.sym[g.sym] = g
      self.sym_up[g.sym.upper()] = g
      for a in g.ali:
        if a not in ['','-']: 
          self.ali[a] = g
          self.ali_up[a.upper()] = g
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
        g = self.by_gid(attr[a])
        if g != None : return g
      elif au in Ens:
        g = self.by_ens(attr[a])
        if g != None : return g
      else: 
        namelist.append(attr[a])
    for n in namelist:
        g = self.by_sym(n)
        if g == None : 
          if self.sym_rectify and n in self.sym_err:
            g = self.by_sym(self.sym_err[n])
          if g == None : continue
        if g.gid not in gids: gids[g.gid] = 0
        gids[g.gid] += 1
    if len(gids) == 0 : return None
    gid = sorted(gids.items(), key=lambda d:d[1])[-1][0]
    return self.gid[gid]

  def by_sym(self, sym):
    sym = str(sym)
    if sym in self.sym: return self.sym[sym]
    if sym in self.ali: return self.ali[sym]
    su = sym.upper()
    if su in self.sym_up: return self.sym_up[su]
    if su in self.ali_up: return self.ali_up[su]
    if su[0:3] == 'LOC' and su[3:] in self.gid:
      return self.gid[su[3:]]
    return None
  def by_gid(self, gid):
    gid = str(gid)
    if gid in self.gid: return self.gid[gid]
    return None
  def by_ens(self, ens):
    ens = str(ens)
    if ens in self.ens: return self.ens[ens]
    return None
  sym_rectify = True
  sym_err = {'MT-ATP6':'MTATP6','MT-CYB':'MTCYB','MT-ND4':'MTND4','MT-CO2':'MTCO2','MT-CO3':'MTCO3','MT-ND3':'MTND3','MT-ND5':'MTND5',
            'PERPL':'PERP','MT-ATP8':'MTATP8','MT-ND1':'MTND1','MT-CO1':'MTCO1','MT-ND2':'MTND2',
            '1-Sep':'SEPT1','2-Sep':'SEPT2','3-Sep':'SEPT3','4-Sep':'SEPT4','5-Sep':'SEPT5','6-Sep':'SEPT6','7-Sep':'SEPT7',
            '8-Sep':'SEPT8','9-Sep':'SEPT9','10-Sep':'SEPT10','11-Sep':'SEPT11','12-Sep':'SEPT12','14-Sep':'SEPT14','15-Sep':'SEP15',
            '1-Mar':'MARCH1','2-Mar':'MARCH2','3-Mar':'MARCH3','4-Mar':'MARCH4','5-Mar':'MARCH5','6-Mar':'MARCH6','7-Mar':'MARCH7',
            '8-Mar':'MARCH8','9-Mar':'MARCH9','10-Mar':'MARCH10','11-Mar':'MARCH11','1-Feb':'FEB1','2-Feb':'FEB2','3-Feb':'FEB3',
            '4-Feb':'FEB4','6-Feb':'FEB6','5-Feb':'FEB5','7-Feb':'FEB7','9-Feb':'FEB9','10-Feb':'FEB10','11-Feb':'FEB11',
            '1-Apr':'APR-1','2-Apr':'APR-2','3-Apr':'APR-3','1-May':'MAY1','1-Oct':'OCT1','2-Oct':'OCT2','3-Oct':'OCT3',
            '4-Oct':'OCT4','6-Oct':'OCT6','1-Nov':'NOV1','2-Nov':'NOV2','1-Dec':'DEC1','2-Dec':'DEC2',
            '41153':'SEPT1','41154':'SEPT2','41155':'SEPT3','41156':'SEPT4','41157':'SEPT5','41158':'SEPT6','41159':'SEPT7',
            '41160':'SEPT8','41161':'SEPT9','41162':'SEPT10','41163':'SEPT11','41164':'SEPT12','41166':'SEPT14','41167':'SEP15',
            '40969':'MARCH1','40970':'MARCH2','40971':'MARCH3','40972':'MARCH4','40973':'MARCH5','40974':'MARCH6','40975':'MARCH7',
            '40976':'MARCH8','40977':'MARCH9','40978':'MARCH10','40979':'MARCH11','42248':'SEPT1','42249':'SEPT2','42250':'SEPT3',
            '42251':'SEPT4','42252':'SEPT5','42253':'SEPT6','42254':'SEPT7','42255':'SEPT8','42256':'SEPT9','42257':'SEPT10',
            '42258':'SEPT11','42259':'SEPT12','42261':'SEPT14','42262':'SEP15','42064':'MARCH1','42065':'MARCH2','42066':'MARCH3',
            '42067':'MARCH4','42068':'MARCH5','42069':'MARCH6','42070':'MARCH7','42071':'MARCH8','42072':'MARCH9','42073':'MARCH10',
            '42074':'MARCH11','42036':'FEB1','42037':'FEB2','42038':'FEB3','42039':'FEB4','42041':'FEB6','42040':'FEB5',
            '42042':'FEB7','42044':'FEB9','42045':'FEB10','42046':'FEB11','42095':'APR-1','42096':'APR-2','42097':'APR-3',
            '42125':'MAY1','42278':'OCT1','42279':'OCT2','42280':'OCT3','42281':'OCT4','42283':'OCT6','42309':'NOV1','42310':'NOV2',
            '42339':'DEC1','42340':'DEC2','43352':'SEPT9','43345':'SEPT2','43350':'SEPT7'
           }
