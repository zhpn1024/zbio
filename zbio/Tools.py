
Nt=['A','C','G','T']

def rc(seq):
	comps = {'A':"T", 'C':"G", 'G':"C", 'T':"A",
		   'B':"V", 'D':"H", 'H':"D", 'K':"M",
		   'M':"K", 'R':"Y", 'V':"B", 'Y':"R",
		   'W':'W', 'N':'N', 'S':'S'}
	return ''.join([comps[x] for x in seq.upper()[::-1]])

#def u2t(seq):
#	return ''.join(x=='U'?'T':x for x in seq.upper())
	
#def t2u(seq):
#	return ''.join(x=='T'?'U':x for x in seq.upper())
	
def overlap(A, B):
	'''
	if A is overlapping with B.
	'''
	if(A.chr != B.chr) : return False
	if (A.stop < B.start) : return False
	if (B.stop < A.start) : return False
	return True

def inside(A, B):
	'''
	if A is inside B.
	'''
	if(A.chr != B.chr) : return False
	if (B.stop < A.stop) : return False
	if (A.start < B.start) : return False
	return True
	
def distance(A,B):
	if A.chr!=B.chr: return None
	if overlap(A,B): return 0
	if A.start<B.start : return B.start-A.stop
	else: return A.start-B.stop
	
def centerin(A,B):
	'''
	If center of A is in B
	'''
	if overlap(A,B):
		return B.start<=A.center()<=B.stop
	else: return False

def end5(A):
	if A.strand=='-':
		return A.stop
	else:
		return A.start

def end3(A):
	if A.strand=='-':
		return A.start
	else:
		return A.stop


def faIter(file):
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

def coverIter(bedIter, weight= lambda x: 1): 
	current = []
	chr = ""
	start = 0
	clen = 0
	for bed in bedIter:
		bed.weight = weight(bed)
		#print bed
		if clen == 0 :
			chr = bed.chr
			start = bed.start
			current.append(bed)
			clen = 1
		elif chr == bed.chr:
			if start == bed.start: pass
			else:
				assert start < bed.start, "Records must be sorted!\n"
				while clen > 0 and current[0].stop <= bed.start:
					cover = 0
					for i in range(clen):
						cover += current[i].weight
					stop = current[0].stop
					yield (chr, start, stop, cover)
					#print clen
					start = stop
					i = 0
					while i < clen and current[i].stop <= stop:
						i += 1
					current[0:i] = []
					clen = len(current)
					#print "clean", i, clen
				if clen > 0 and current[0].stop > bed.start:
					cover = 0
					for i in range(clen):
						cover += current[i].weight
					stop = bed.start
					yield (chr, start, stop, cover)
					#print clen
					start = stop
			i = clen - 1
			while i >= 0 and bed.stop < current[i].stop:
				i -= 1
			current[i+1:i+1] = [bed]
			clen = len(current)
			if clen == 1: start = bed.start
		else:
			while clen > 0 :
				cover = 0
				for i in range(clen):
					cover += current[i].weight
				stop = current[0].stop
				yield (chr, start, stop, cover)
				start = stop
				i = 0
				while i < clen and current[i].stop <= stop:
					i += 1
				#for i in range(clen):
					#if current[i].stop > stop: break
				current[0:i] = []
				clen = len(current)
			chr = bed.chr
			start = bed.start
			current.append(bed)
			clen = 1
	
	while clen > 0 :
		cover = 0
		for i in range(clen):
			cover += current[i].weight
		stop = current[0].stop
		yield (chr, start, stop, cover)
		start = stop
		i = 0
		while i < clen and current[i].stop <= stop:
			i += 1
		#for i in range(clen):
			#if corrent[i].stop > stop: break
		current[0:i] = []
		clen = len(current)
		
def overlapIter(bedIterA, bedIterB, func=overlap, ignoreStrand = True, counts = [0,0]):
	lst = []
	ac = bedIterA.next()
	counts[0] += 1
	Aend = False
	for b in bedIterB:
		#print len(lst)
		counts[1] += 1
		j = -1
		for i in range(len(lst)):
			if lst[i].chr < b.chr:
				j = i
			elif lst[i].chr == b.chr:
				if lst[i].stop < b.start:
					j = i
			#if lst[i].id == 'piR-mmu-10797635': print b
				if func(lst[i], b):
					if ignoreStrand or lst[i].strand == b.strand :
						yield (lst[i], b)
		lst[0:(j+1)] = []
		if ac.chr > b.chr or (ac.chr == b.chr and ac.start > b.stop) :
			continue
		if func(ac, b):
			if ignoreStrand or ac.strand == b.strand :
				yield (ac, b)
		if Aend == False :
			lst.append(ac)
		c = 0
		for a in bedIterA:
			counts[0] += 1
			assert a.chr > ac.chr or a.start >= ac.start, "Records must be sorted!\n"+str(a)
			ac = a
			c += 1
			#lst.append(a)
			#if a.id == 'piR-mmu-10797635': 
			#	print b
			#	return
			if a.chr < b.chr: continue
			elif a.chr == b.chr:
				if a.stop < b.start: continue
				if func(a, b):
					if ignoreStrand or a.strand == b.strand :
						yield (a, b)
				elif a.start > b.stop: break
			else: break
			lst.append(a)
		if c == 0 : Aend = True
		#print len(lst)

def randOverlapIter(bedIterA, bedListB, func=overlap, ignoreStrand = True):
	m = len(bedListB)
	for a in bedIterA:
		i0 = 0
		i1 = m
		while i1 - i0 > 1 :
			i = (i1 - i0) / 2
			if a > bedListB[i] : i0 = i
			else: i1 = i
		if func(a, bedListB[i0]):
			if ignoreStrand or a.strand == bedListB[i0].strand :
				yield (a, bedListB[i0])
		if func(a, bedListB[i1]):
			if ignoreStrand or a.strand == bedListB[i0].strand :
				yield (a, bedListB[i1])

def bed2Seq(seq, bed):
  s = ''
  for e in bed.exons():
    es = seq[e.start:e.stop]
    if e.strand == '-':
      es = rc(es)
    s += es
  return s
