#!/usr/bin/env python

'''
getSplit.py

convert annotations to genome split sites for pseudo reads sampling.
Designed by Peng Zhang on Apr. 2015
'''

import sys
import getopt
import gzip
from zbio import Bed
from zbio import Gtf

use_message = '''

To be added.
'''

class Params:
  def __init__(self, num_threads = 1):
    self.num_threads = num_threads
  def parse_options(self, argv):
    try:
      opts, args = getopt.getopt(argv[1:], "hp:",["help", "num-threads="])
    except:
      print(use_message)
      exit(1)
    self.args = args
    for option, value in opts:
      if option in ("-p", "--num-threads"):
        #print(option + value)
        self.num_threads = int(value)
      if option in ("-h", "--help"):
        print(use_message)
        exit(1)
    return args
  def check(self):
    if self.num_threads < 1 :
      die("Error: arg to --num-threads must be greater than 0")


def fopen(filename):
  ns = filename.split('.')
  if ns[-1]=='gz' or ns[-1]=='gzip':
    fin = gzip.open(filename,'r')
  else:
    fin = open(filename,'r')
  return fin

class split:
  def __init__(self):
    self.spls = {}
  def add(self, chr, pos, type, strand):
    if chr not in self.spls:
      self.spls[chr] = {}
    if pos not in self.spls[chr]:
      self.spls[chr][pos] = [type, strand]
    else :
      if self.spls[chr][pos][0] < 2 and type != self.spls[chr][pos][0]:
        self.spls[chr][pos][0] = 2
      if len(self.spls[chr][pos][1]) < 2 and strand != self.spls[chr][pos][1]:
        self.spls[chr][pos][1] = '+-'
  def loadJunc(self, fin):
    for l in fin:
      lst = l.strip().split()
      split.add(lst[0], int(lst[1]), 0, lst[3])
      split.add(lst[0], int(lst[2]), 1, lst[3])
  def loadBed12(self, fin):
    for b in Bed.bed12Iter(fin):
      strand = b.strand
      if strand not in ['+', '-']:
        strand = '+-'
      for e in b.exons():
        self.add(e.chr, e.start, 1, strand)
        self.add(e.chr, e.stop, 0, strand)
  def loadGtf(self, fin):pass
  def export(self, fout = sys.stdout):
    chrs = self.spls.keys()
    chrs.sort()
    for chr in chrs:
      ss = self.spls[chr].keys()
      ss.sort()
      for pos in ss:
        stat, strand = self.spls[chr][pos]
        fout.write(chr +'\t'+ str(pos) +'\t'+ str(stat) +'\t'+ strand +'\n')

def getSplit(split = split(), argv = None):
  params = Params()
  if argv is None:
    argv = sys.argv
  args = params.parse_options(argv)
  params.check()
  #print params.num_threads

  if len(args) == 0:
    print(use_message)
    exit(1)
  for arg in args:
    ns = arg.strip().split('.')
    i = -1
    if ns[-1] in ['gz', 'gzip']: i = -2
    if ns[i] == 'junc' :
      split.loadJunc(fopen(arg))
    elif ns[i] == 'bed' :
      split.loadBed12(fopen(arg))
    else:
      print 'Unknown file format: ' + arg
  return split

if __name__ == '__main__':
  split = getSplit()
  split.export()
  

