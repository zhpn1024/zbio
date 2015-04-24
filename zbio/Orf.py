from zbio import Tools

codonSize = 3
cstart = ['ATG']
cstartlike = ['TTG', 'CTG', 'GTG', 'AAG', 'AGG', 'ACG', 'ATT', 'ATA', 'ATC']
cstop = ['TGA', 'TAA', 'TAG']

senseframe = [1, 2, 3]
antiframe = [-1, -2, -3]
frame = [1, 2, 3, -1, -2, -3]

def findORF(seq, strand = '+', altcstart = False) :
  seq = seq.upper().replace('U','T')
  if strand == '+' : fr = senseframe
  elif strand == '-' : 
    fr = antiframe
    antiseq = Tools.rc(seq)
  else: 
    fr = frame
    antiseq = Tools.rc(seq)
  if altcstart: cs = cstart + cstartlike
  else: cs = cstart
  
  length = len(seq)
  for f in fr:
    if f > 0:
      s = seq
      fa = f
    else: 
      s = antiseq
      fa = -f
    orfstart = orfstop = -1
    for i in range(fa-1, length, codonSize):
      try: codon = s[i:i+codonSize]
      except: break
      if codon in cs:
        if orfstart < 0: orfstart = i
      elif codon in cstop:
        orfstop = i + codonSize
        if orfstart >= 0:
          orflen = (i - orfstart) / 3
          if f > 0:
            yield orfstart, orfstop, f, orflen
          else:
            yield length - orfstop, length - orfstart, f, orflen
          orfstart = orfstop = -1
          
          