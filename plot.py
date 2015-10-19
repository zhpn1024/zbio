import matplotlib
matplotlib.use('pdf')
from matplotlib.pylab import *

def plottrans(t, ypos, intv, r = [0.1, 0.3], color = 'blue',rid = -0.5):
  plot([t.start,t.stop],[ypos,ypos],color=color)

  x = []
  y = []
  for i in range(t.start+intv/2, t.stop-intv/3, intv):
    x.append(i)
    y.append(ypos)
  arr = '>'
  if t.is_reverse() : arr = "<"
  plot(x,y,'w'+arr)
  x = [[],[]]
  y = [[],[]]
  cds = 0
  exons = t.exons[:]
  exons.sort()
  ts1 = t.thick_start
  ts2 = t.thick_stop
  orf = t(start = ts1, stop = ts2)
  #print ts1,ts2
  for e in exons:
    for es in e-orf : 
      x[0].append(es.start)
      y[0].append(len(es))
    for ei in e.intersect(orf):
      x[1].append(ei.start)
      y[1].append(len(ei))
    #if e.start < ts1 < e.stop :
      #x[0].append(e.start)
      #y[0].append(ts1-e.start)
      #if e.start < ts2 < e.stop :
        #x[1].append(ts1)
        #y[1].append(ts2-ts1)
        #x[0].append(ts2)
        #y[0].append(e.stop-ts2)
      #else :
        #x[1].append(ts1)
        #y[1].append(e.stop-ts1)
        #cds = 1
    #elif e.start < ts2 < e.stop :
      #x[1].append(e.start)
      #y[1].append(ts2-e.start)
      #x[0].append(ts2)
      #y[0].append(e.stop-ts2)
      #cds = 0
    #else :
      #x[cds].append(e.start)
      #y[cds].append(e.stop-e.start)
  #x.append(e.start)
  #y.append(len(e))
  bar(x[0],[r[0]*2]*len(x[0]),width=y[0],bottom=ypos-r[0],edgecolor=color,color=color)
  bar(x[1],[r[1]*2]*len(x[1]),width=y[1],bottom=ypos-r[1],edgecolor=color,color=color)
  text((t.start+t.stop)/2,ypos+rid,t.id)
  #print x[1],y[1]

def save(file):
  savefig(file)
