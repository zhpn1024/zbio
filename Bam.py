
import pysam

class Bamfile(pysam.Samfile):
    
    def __repr__(self):
        return 'pysam.Samfile '+self.filename

class Bam():#AlignedRead
    
    def __init__(self,read,bamfile):
        self.read=read
        self.ref=bamfile.references
        
    #def __getattr__(self,attr):
        #a=getattr(self,attr)
        #return a
    
    @property
    def chr(self): 
        return self.ref[self.read.tid]
    @property
    def start(self):
        return self.read.pos
    @property
    def stop(self):
        return int(self.read.aend)
    @property
    def id(self): #Bed3 only
        return self.read.qname
    @property
    def score(self): #Bed3 only
        return 0.0
    @property
    def strand(self): #Bed3 only
        if self.read.is_reverse: return '-'
        else: return '+'
    
    def __len__(self): #All bed
        return self.read.alen
    @property
    def length(self):
        return len(self)
    
    def cdna_length(self): 
        return self.read.qlen
    
    def center(self): #Middle point, NEED REVISE!!
        return (self.start+self.stop)/2.0
    
    def is_reverse(self): #All bed
        return self.read.is_reverse
    
    @property
    def end5(self): #5' end, all bed
        if not self.read.is_reverse: return self.start
        else : return self.stop
    @property
    def end3(self): #3' end, all bed
        if not self.read.is_reverse: return self.stop
        else : return self.start
        
    