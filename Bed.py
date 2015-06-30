

class Bed3:
  Header=('chr','start','stop')
  Format=(str,int,int)
  n=3
  def __init__(self, x, bin = False): #Fit all bed
    if type(x)==str:
      l = x.strip().split('\t')
      if bin:
        l[0:1] = []
    elif type(x)==dict:
      l=[x[i] for i in self.Header]
    else:
      l=x
    lst=[]
    for i in range(self.n):
      try:
        lst.append(self.Format[i](l[i]))
      except:
        lst.append('.')
    self.items=tuple(lst)
  def __str__(self): #Bed3 and Bed6
    return '\t'.join(map(str,self.items))
    #return self.chr+'\t'+str(self.start)+'\t'+str(self.stop)
  def __repr__(self): #All bed
    return 'Bed'+str(self.n)+' object:\n'+str(self)+'\n'
  def shortStr(self):
    return "%s:%d-%d:%s" % (self.chr, self.start, self.stop, self.strand)
  @property
  def chr(self): #All bed
    return self.items[0]
  @property
  def start(self): #All bed
    return self.items[1]
  @property
  def stop(self): #All bed
    return self.items[2]
  @property
  def id(self): #Bed3 only
    return "noname"
  @property
  def score(self): #Bed3 only
    return 0.0
  @property
  def strand(self): #Bed3 only
    return "."
  def __len__(self): #All bed
    return self.stop-self.start
  @property
  def length(self):
    return len(self)
  @property
  def end(self): #All bed
    return self.stop
  @property
  def dict(self): #All bed
    return dict([(self.Header[i],self.items[i]) for i in range(self.n)])
  @property
  def end5(self): #5' end, all bed
    if self.strand != '-' : return self.start
    else : return self.stop
  @property
  def end3(self): #3' end, all bed
    if self.strand != '-' : return self.stop
    else : return self.start
  def is_reverse(self): #All bed
    return self.strand=='-'
  def bed(self): #Bed copy, Bed3 only
    return Bed3(self)
  def __getitem__(self,i): #All bed
    return self.items[i]
  def __cmp__(self,other): #All bed
    return cmp(self.chr,other.chr) or cmp(self.start,other.start) or cmp(self.stop,other.stop)
  def headerline(self,sep='\t'):# Header string, fit all bed
    sep=str(sep)
    return sep.join(self.Header)
  def head(self,n=10): #Bed3 only
    if n>len(self):
      n=len(self)
    return Bed3([self.chr,self.start,self.start+n])
  def tail(self,n=10): #Bed3 only
    if n>len(self):
      n=len(self)
    return Bed3([self.chr,self.stop-n,self.stop])
  def __call__(self,**args): #New bed that modified from self. All bed
    cp=self.bed()
    if args=={}:
      return cp
    d=self.dict
    for k in args.keys():
      if k in d.keys():
        d[k]=args[k]
      #else: assert 0, k+' not defined!\n'
    cp.__init__(d)
    return cp
  def center(self): #Middle point, Bed3 and Bed6
    return (self.start+self.stop)/2.0
  def cdna_length(self): #Bed3 and Bed6
    return len(self)
  @property
  def exons(self): #Bed3 and Bed6
    l=[]
    l.append(self(id=self.id+"_Exon_1"))
    return l
  def is_contain(self,p): #if i in bed, all bed
    return self.start<=p<=self.stop
  def is_exon(self,p): #Bed3 and Bed6
    return self.is_contain(p)
  def is_intron(self,p): #False, Bed3 and Bed6
    return False
  def is_sense(self,other): #In same strand? All bed
    return self.strand==other.strand
  def cdna_pos(self,p): #Bed3 and Bed6
    if self.is_contain(p): return abs(p-self.end5)
    else: return None

class Bed6(Bed3):
  Header=('chr','start','stop','id','score','strand')
  Format=(str,int,int,str,str,str)
  n=6
  #def __init__(self,x):
  #def __str__(self):

  @property
  def id(self): #Bed6 and Bed12
    return self.items[3]
  @property
  def score(self): #Bed6 and Bed12
    return self.items[4]
  @property
  def strand(self): #Bed6 and Bed12
    return self.items[5]
  
  def bed(self): #Only Bed6
    return Bed6(self)
  
  def head(self, n=10): #Bed6 and Bed12
    if n>len(self):
      n=len(self)
    if self.strand != '-':
      start=self.start
      stop=start+n
    else:
      stop=self.stop
      start=stop-n
      if start<0: start=0
    return Bed6([self.chr,start,stop,self.id+'_Head_'+str(n),self.score,self.strand])
  
  def tail(self,n=10): #Bed6 and Bed12
    if n>len(self):
      n=len(self)
    if self.strand == '-':
      start=self.start
      stop=start+n
    else:
      stop=self.stop
      start=stop-n
      if start<0: start=0
    return Bed6([self.chr,start,stop,self.id+'_Tail_'+str(n),self.score,self.strand])
  

def comTotup(s): #comma string to tuple
  if type(s)==tuple: return s
  return tuple(map(int,s.strip().strip(',').split(',')))
def tupTocom(t): #tuple to comma string
  if type(t)==str: return t
  return ','.join(map(str,t))+','

class Bed12(Bed6):
  
  Header=('chr','start','stop','id','score','strand',"cds_start","cds_stop","itemRgb","blockCount","blockSizes","blockStarts")
  Format=(str,int,int,str,str,str,int,int,comTotup,int,comTotup,comTotup)
  n=12
  
  @property
  def cds_start(self): #Bed12
    return self.items[6]
  @property
  def cds_stop(self):
    return self.items[7]
  @property
  def itemRgb(self):
    return self.items[8]
  @property
  def blockCount(self):
    return self.items[9]
  @property
  def blockSizes(self):
    return self.items[10]
  @property
  def blockStarts(self):
    return self.items[11]

  @property
  def blockStops(self):
    return tuple(map(lambda x,y:x+y,self.blockStarts,self.blockSizes))
  def blockStop(self,i): 
    return self.blockStarts[i]+self.blockSizes[i]
  
  def __str__(self): #Bed string, Bed12 only
    l=list(self.items)
    l[8]=tupTocom(l[8])
    l[10]=tupTocom(l[10])
    l[11]=tupTocom(l[11])
    return '\t'.join(map(str,l))
  
  def bed(self): #Bed copy, Bed12 only
    return Bed12(self)
  
  def cdna_length(self): #Bed12
    l = 0
    for i in self.blockSizes:
      l += i
    return l
  
  def center(self): #Middle point, Bed12 only
    len=self.cdna_length()/2.0
    for i in range(self.blockCount):
      len-=self.blockSizes[i]
      if len<0:
        break
    return self.start+self.blockStarts[i]+self.blockSizes[i]+len
  @property
  def exons(self): #Bed12
    a=[]
    if self.strand=="-":
      step=-1
      j=self.blockCount
    else:
      step=1
      j=1
    for i in range(self.blockCount):
      start=self.start+self.blockStarts[i]
      end=self.start+self.blockStarts[i]+self.blockSizes[i]
      id=self.id+"_Exon_"+str(j)
      j+=step
      a.append(Bed6((self.chr,start,end,id,self.score,self.strand)))
    if self.strand=="-":
      return a[::-1]
    else:
      return a
  @property
  def introns(self):
    a=[]
    if self.strand=="-":
      step=-1
      j=self.blockCount-1
    else:
      step=1
      j=1
    for i in range(self.blockCount-1):
      start=self.start+self.blockStarts[i]+self.blockSizes[i]
      end=self.start+self.blockStarts[i+1]
      id=self.id+"_Intron_"+str(j)
      j+=step
      a.append(Bed6((self.chr,start,end,id,self.score,self.strand)))
    if self.strand=="-":
      return a[::-1]
    else:
      return a
  
  def cdna_pos(self, p):
    if p<self.start or p>self.stop:
      return None
    p1=p-self.start
    pos=0
    l=range(self.blockCount)
    if self.is_reverse():
      l=l[::-1]
    for i in l:
      if self.blockStarts[i]<=p1<=self.blockStop(i):
        if self.is_reverse():
          pos+=self.blockStop(i)-p1
        else:
          pos+=p1-self.blockStarts[i]
        return pos
      else:
        pos+=self.blockSizes[i]
    else:
      return None
  
  def genome_pos(self, p, bias=0):
    m = self.cdna_length()
    if p < 0 or p > m: return None
    if p == 0 : return self.end5
    if p == m: return self.end3
    p1=p
    pos=self.start
    l=range(self.blockCount)
    if self.is_reverse():
      l=l[::-1]
    for i in l:
      if self.blockSizes[i]-p1>=bias:
        if self.is_reverse():
          pos+=self.blockStop(i)-p1
        else:
          pos+=p1+self.blockStarts[i]
        return pos
      else:
        p1-=self.blockSizes[i]
    else:
      return None 
  
  #def cdna_type(p):
    
  def type(self, p, genome = True):
    if genome == False:
      pg = self.genome_pos(p)
      return self.type(pg)
    if type(p) != int : return None
    if p > self.stop or p < self.start : return None
    if p < self.cds_start:
      if self.is_reverse(): return "3UTR"
      else: return "5UTR"
    elif  p > self.cds_stop:
      if self.is_reverse(): return "5UTR"
      else: return "3UTR"
    else: return "CDS"

  def is_exon(self,p): #If in exon, Bed12
    if not self.is_contain(p):
      return False
    l=range(self.blockCount)
    if self.is_reverse():
      l=l[::-1]
    p1=p-self.start
    for i in l:
      if self.blockStarts[i]<=p1<=self.blockStop(i):
        return True
    else:
      return False
    #return self.is_contain(p)
  def is_intron(self,p): #If in intron
    #return self.is_contain(p) and not self.is_exon(p)
    if not self.is_contain(p):
      return False
    l=range(1,self.blockCount)
    if self.is_reverse():
      l=l[::-1]
    p1=p-self.start
    for i in l:
      if self.blockStarts[i] > p1 > self.blockStop(i-1):
        return True
    else:
      return False

def bed3Iter(file, bin = False):
  for l in file:
    try:
      yield Bed3(l, bin)
    except ValueError:
      pass


def bed6Iter(file, bin = False):
  for l in file:
    try:
      yield Bed6(l, bin)
    except ValueError:
      pass

def bed12Iter(file, bin = False):
  for l in file:
    try:
      yield Bed12(l, bin)
    except ValueError:
      pass

class refGene(Bed12):
  def __init__(self, x):
    l = x.strip().split('\t')
    self.exonStarts = comTotup(l[9])
    self.exonEnds = comTotup(l[10])
    start = int(l[4])
    blockSizes = [self.exonEnds[i] - self.exonStarts[i] for i in range(len(self.exonStarts))]
    blockStarts = [self.exonStarts[i] - start for i in range(len(self.exonStarts))]
    l2 = [l[2],l[4],l[5],l[1],l[11],l[3],l[6],l[7],'0',l[8],tupTocom(blockSizes),tupTocom(blockStarts)]
    Bed12.__init__(self, l2)
    self.bin = int(l[0])
    self.name2 = l[12]
    self.exonFrames = comTotup(l[15])
    self.cdsStartStat = l[13]
    self.cdsEndStat = l[14]
    
  def __repr__(self): 
    return 'refGene object:\n'+str(self)+'\n'
    
def refGeneIter(file):
  for l in file:
    try:
      yield refGene(l)
    except ValueError:
      pass
    
def refFlatIter(file):
  for l in file:
    lst = l.strip().split()
    starts = comTotup(lst[9])
    stops = comTotup(lst[10])
    st1 = tuple(map(lambda x: x - int(lst[4]), starts))
    sizes = tuple(map(lambda x, y: y - x, starts, stops))
    lstb = [lst[2],lst[4],lst[5],lst[1],lst[0],lst[3],lst[6],lst[7],"0,0,0",lst[8],sizes,st1]
    bed = Bed12(lstb)
    yield bed
    
def shortBed(s, name = ''): #like chr1:1-200:+
  lst = s.strip().split(":")
  l1 = lst[1].split("-")
  if name == "": name = s.strip()
  if len(lst) >= 3: strand = lst[2]
  else: stand = '.'
  return Bed6([lst[0], l1[0], l1[1], name, 0, strand])