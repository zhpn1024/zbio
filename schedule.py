'''
Task scheduling
Copyright (c) 2023 Peng Zhang <zhpn1024@163.com>
'''
import os
from os.path import isfile
import time
import sys

class Pipeline:

  def __init__(self, pipe = None):
    if pipe is not None: self.pipe = pipe
    else: self.pipe = []

  def add_step(self, name, command, cores = 1, limit = None):
    self.pipe.append([name, command, cores, limit])

  def __len__(self):
    return len(self.pipe)

  def __getitem__(self, i):
    return self.pipe[i]

  def cores(self, i):
    return self.pipe[i][2]

  def limit(self, i):
    return self.pipe[i][3]

class Job:

  def __init__(self, name, pipeline, batch = ''):
    self.id = name
    self.pipeline = pipeline
    self.batch = batch
    self.progress = [0] * len(self.pipeline)
    self.current = 0
    self.label = ''

  def __len__(self):
    return len(self.pipeline)

  def status(self):
    if self.current >= len(self): return 'Job complete' # self.progress[-1] == 2: return 'Job complete'
    else:
      #for i in range(len(self)):
        #if self.progress[i] != 2: break
      if self.progress[self.current] == 0: return '{} {} waiting'.format(self.current, self.pipeline[self.current][0])
      elif self.progress[self.current] == 1: return '{} {} running'.format(self.current, self.pipeline[self.current][0])
      else: return '{} {} error'.format(self.current, self.pipeline[self.current][0])

  def submit(self, i = None):
    if i is None: i = self.current
    else: self.current = i
    os.popen(self.pipeline[i][1].format(self.id) + ' >> log_{}_{}_{}_running.txt 2>&1 &'.format(self.current, self.pipeline[i][0], self.id)) # background run?
    self.progress[self.current] = 1

  def is_complete(self):
    return self.current >= len(self.progress)

  def is_error(self):
    if self.is_complete(): return False
    return self.progress[self.current] == -1

  def is_running(self):
    if self.is_complete(): return False
    return self.progress[self.current] == 1

  def is_waiting(self):
    if self.is_complete(): return False
    return self.progress[self.current] == 0

  def update(self):
    if self.is_complete(): return 3 # complete
    if self.is_error(): return -1
    if self.progress[self.current] == 0: return 0 # waiting
    out = os.popen('tail -n1 log_{}_{}_{}_running.txt'.format(self.current, self.pipeline[self.current][0], self.id))
    s = out.read()
    if s.startswith('Complete'):
      self.progress[self.current] = 2
      print('{} Step completed {} {} {}_{}'.format(time.ctime(), self.id, self.label, self.current, self.pipeline[self.current][0]))
      out = os.popen('mv log_{0}_{1}_{2}_running.txt log_{0}_{1}_{2}_complete.txt'.format(self.current, self.pipeline[self.current][0], self.id))
      self.current += 1
      if self.is_complete(): return 3
      else:
        self.progress[self.current] = 0
        return 2
      #return 2 # complete
    elif s.startswith('Uncomplete'):
      self.progress[self.current] = -1
      os.system('echo "{} {} {} {} UnComplete `date +%Y/%m/%d--%H:%M`" >>log_error_{}.lst'.format(time.ctime(), self.id, self.current, self.pipeline[self.current][0], self.batch))
      out = os.popen('mv log_{0}_{1}_{2}_running.txt log_{0}_{1}_{2}_error.txt'.format(self.current, self.pipeline[self.current][0], self.id))
      print('{} Step error {} {} {}_{}'.format(time.ctime(), self.id, self.label, self.current, self.pipeline[self.current][0]))
      return -1 # error
    else: return 1 # running

  def using_cores(self):
    if self.progress[self.current] != 1: return 0 # not running
    return self.pipeline.cores(self.current)

  def reset(self, i = None):
    if i is None: i = self.current
    else: self.current = i
    self.progress[self.current] = 0
    print('reset {} current = {}'.format(self.id, i))

class Scheduler:

  def __init__(self, name, max_cores, pipeline, samples, priority = 'sample'):
    self.id = name
    self.max_cores = max_cores
    #self.cores = 0
    self.pipeline = pipeline
    self.jobs = [Job(s, pipeline, name) for s in samples]
    self.steps = len(pipeline)
    self.priority = priority
    n = len(samples)
    for i, j in enumerate(self.jobs):
      j.label += '{}/{}'.format(i+1, n)

  def run(self, wait = 100):
    running = True
    while running:
      running = False
      changed = False

      if isfile('exec_{}.txt'.format(self.id)):
        l = open('exec_{}.txt'.format(self.id)).read()
        print(l)
        try: exec(l)
        except Exception as e: print('Error exec exec_{}.txt! {}'.format(self.id, str(e)))
        os.system('rm exec_{}.txt'.format(self.id))
        changed = True

      reset = {}
      if isfile('reset_{}.txt'.format(self.id)):
        resetfile = open('reset_{}.txt'.format(self.id))
        for l in resetfile:
          lst = l.strip().split()
          try: reset[lst[0]] = int(lst[1])
          except Exception as e: print('Error reset reset_{}.txt! {} {} {}'.format(self.id, lst[0], lst[1], str(e)))
        resetfile.close()
        os.system('rm reset_{}.txt'.format(self.id))
        changed = True

      dstep, drun = {}, {}
      cores = 0
      skip = {}
      for j in self.jobs:
        if j.id in reset:
          j.reset(reset[j.id])
        status = j.update()
        if j.is_complete(): continue
        else: running = True
        if status == 1: cores += j.using_cores()
        elif status == 2:
          changed = True
          skip[j.current-1] = 1 # hold 1 round
        elif status == -1:
          changed = True
        if j.current not in dstep:
          dstep[j.current] = []
          drun[j.current] = 0
        dstep[j.current].append(j)
        if j.is_running(): drun[j.current] += 1

      if self.priority == 'sample':
        for j in self.jobs:
          if not j.is_waiting(): continue
          i = j.current
          if i in skip: continue
          if cores + self.pipeline.cores(i) > self.max_cores: continue
          limit = self.pipeline.limit(i)
          if limit is not None and drun[i] >= limit: continue
          j.submit()
          cores += self.pipeline.cores(i)
          drun[i] += 1
          changed = True
          print('{} Submit {} {} {}_{} {} cores, cores used {}/{}, step limit {}/{}'.format(time.ctime(), j.id, j.label, i, self.pipeline[i][0], self.pipeline.cores(i), cores, self.max_cores, drun[i], limit))
          break

      else:
        for i in range(self.steps):
          if i in skip: continue
          if i not in dstep: continue
          if cores + self.pipeline.cores(i) > self.max_cores: continue
          limit = self.pipeline.limit(i)
          num_running = 0
          for j in dstep[i]:
            if j.is_running(): num_running += 1
          #if limit is not None and num_running >= limit: continue
          for j in dstep[i]:
            if limit is not None and num_running >= limit: break
            if cores + self.pipeline.cores(i) > self.max_cores: break
            if j.is_waiting():
              j.submit()
              cores += self.pipeline.cores(i)
              num_running += 1
              changed = True
              print('{} Submit {} {} {}_{} {} cores, cores used {}/{}, step limit {}/{}'.format(time.ctime(), j.id, j.label, i, self.pipeline[i][0], self.pipeline.cores(i), cores, self.max_cores, num_running, limit))
              #sys.stdout.flush()
              break  # one job submit per step
      if changed: sys.stdout.flush()

      time.sleep(wait)
    print('{} All jobs completed!'.format(time.ctime()))


