#!/usr/bin/env python2

import os,sys,argparse

class MortonSegmentTable:
  """ A table to enumerate the types of Morton segments of various lengths """

  def __init__(self,dim,logfile=None):
    self.dim        = dim
    self.table      = {}
    self.calculated = {}
    self.caseOne    = []
    self.caseWeak   = [None] * dim
    self.logfile    = logfile
    self.typestr    = {-1: "disconnected", 0: "strongly connected"}
    for d in range(0,dim):
      self.typestr[d + 1] = "weakly connected in direction %d" % d
    if not logfile:
      self.logfile  = os.devnull
    return

  def calculate_case_one(self):
    if self.caseOne:
      return
    else:
      self.logfile.write("===\nCalculate child segments of unitary segment\n")
      for i in range (1,2**self.dim + 1):
        if i == 1:
          self.logfile.write("  %d child segments of length 1\n" % 2**self.dim)
          self.caseOne.append((1,2**self.dim,0))
        else:
          self.logfile.write("  Counting child segments of length %d\n" % i)
          disconnected = 0
          connected = 0
          weak = [0] * self.dim
          for j in range (0,2**self.dim - i + 1):
            self.logfile.write("    Considering segment of length %d starting at %d ..." % (i,j))
            first = j;
            last  = j + i - 1;
            diff = last ^ first
            diffDim = int.bit_length(diff) - 1
            shifted = first ^ 2**diffDim
            if shifted < last:
              self.logfile.write(" disconnected\n")
              connected += 1
            elif shifted == last:
              self.logfile.write(" weakly connected in direction %d\n" % diffDim)
              weak[diffDim] += 1
            else:
              self.logfile.write(" strongly connected\n")
              disconnected += 1
          self.logfile.write("  Total disconnected length %d child segments: %d\n" % (i,disconnected))
          self.caseOne.append((i,disconnected,-1))
          self.logfile.write("  Total strongly connected length %d child segments: %d\n" % (i,connected))
          self.caseOne.append((i,connected,0))
          for d in range(0,self.dim):
            self.logfile.write("  Total weakly connected length %d child segments in direction %d: %d\n" % (i,d,weak[d]))
            self.caseOne.append((i,weak[d],d + 1))
      self.logfile.write("===\n")
    return

  def calculate_case_weak(self,way):
    if self.caseWeak[way]:
      return
    else:
      self.caseWeak[way] = []
      self.logfile.write("===\nCalculate child segments of weakly connected segment in direction %d\n" % way)
      counters = [0] * (2**(self.dim + 1) - 1) * 3
      onFace = [0] * 2**(self.dim)
      for i in range(0,2**(self.dim-1)):
        ilow = i & (2**way - 1)
        ihi  = (i^ilow) & (2**self.dim - 1)
        idx = ilow | (ihi << 1)
        for j in range(idx,2**self.dim):
          onFace[j] = i
      for i in range(0,2**self.dim):
        self.logfile.write("  If first cell has %d children in child segment ...\n" % (i + 1))
        for j in range(0,2**self.dim):
          self.logfile.write("    ... and last cell has %d children in child segment ..." % (j + 1))
          childLength = i + j + 2
          iOnFace = onFace[i]
          jOnFace = 2**(self.dim - 1) - 1 - onFace[j]
          if iOnFace < jOnFace:
            self.logfile.write(" no first/last children are neighbors: disconnected\n")
            counters[(childLength - 2) * 3 + 0] += 1
          elif iOnFace == jOnFace and (not i & (2**way)) and (not j & (2**way)):
            self.logfile.write(" only endpoint children are neighbors: weakly connected in direction %d\n" % way)
            counters[(childLength - 2) * 3 + 2] += 1
          else:
            self.logfile.write(" first/last children other than endpoints are neighbors: strongly connected\n")
            counters[(childLength - 2) * 3 + 1] += 1
      for i in range(2,2**(self.dim+1) + 1):
        idx = i - 2
        self.logfile.write("  Total disconnected length %d child segments: %d\n" % (i,counters[idx * 3 + 0]))
        self.caseWeak[way].append((i,counters[idx * 3 + 0],-1))
        self.logfile.write("  Total strongly connected length %d child segments: %d\n" % (i,counters[idx * 3 + 1]))
        self.caseWeak[way].append((i,counters[idx * 3 + 1],0))
        self.logfile.write("  Total weakly connected length %d child segments in direction %d: %d\n" % (i,way,counters[idx * 3 + 2]))
        self.caseWeak[way].append((i,counters[idx * 3 + 2],1+way))
      self.logfile.write("===\n")
    return

  def child_tuple_insert(self,level,length,segtype,count,pad):
    childTuple = (level,length,segtype)
    if count:
      self.logfile.write(pad + "Adding %d segments to (level = %d, length = %d, type = %s)\n" % (count,level,length,self.typestr[segtype]))
    if childTuple in self.table:
      self.table[childTuple] += count
    else:
      self.table[childTuple] = count
    return

  def enumerate_internal(self,level,length,segtype,pad):

    assert length > 0
    assert length <= 2**(self.dim*level)

    nextpad = pad + ' '

    if (level,length,segtype) in self.calculated:
      return self.table[(level,length,segtype)]
    else:
      self.logfile.write(pad + 'Calculating segments on level %d with length %d of type %s ...\n' % (level,length,self.typestr[segtype]))
      if (level,length,segtype) not in self.table:
        self.table[(level,length,segtype)] = 0
      if level == 0:
        if segtype == 0:
          self.table[(level,length,segtype)] = 1
      else:
        # figure out which ancestor segments could add to me
        minNumParents = ((length - 1) / 2**(self.dim)) + 1
        remainder = length % 2**(self.dim)
        if remainder <= 1:
          maxNumParents = length/2**(self.dim) + 1
        else:
          maxNumParents = length/2**(self.dim) + 2
        maxNumParents = min(maxNumParents,2**(self.dim*(level - 1)))
        for numParents in range(minNumParents,maxNumParents+1):
          self.enumerate_internal(level-1,numParents,-1,nextpad)
          self.enumerate_internal(level-1,numParents, 0,nextpad)
          for d in range(0,self.dim):
            self.enumerate_internal(level-1,numParents,1+d,nextpad)

      self.calculated[(level,length,segtype)] = True
      count = self.table[(level,length,segtype)]
      self.logfile.write(pad + '... counted %d segments on level %d with length %d of type %s\n' % (count,level,length,self.typestr[segtype]))
      if count:
        if length == 1:
          self.logfile.write(pad + 'Adding to descendant counts by multiplying unitary segment counts\n')
          self.calculate_case_one()
          for (childLength,childCount,childType) in self.caseOne:
            childAdd = count * childCount;
            self.child_tuple_insert(level + 1,childLength,childType,childAdd,nextpad)
        else:
          if segtype == -1:
            self.logfile.write(pad + 'Adding to disconnected descendant counts\n')
            for i in range(2,2**(self.dim + 1) + 1):
              childLength = (length - 2) * 2**(self.dim) + i
              childAdd    = count * min(i - 1,2**(self.dim + 1) + 1 - i)
              self.child_tuple_insert(level + 1,childLength,-1,childAdd,nextpad)
          elif segtype == 0:
            self.logfile.write(pad + 'Adding to strongly connected descendant counts\n')
            for i in range(2,2**(self.dim + 1) + 1):
              childLength = (length - 2) * 2**(self.dim) + i
              childAdd    = count * min(i - 1,2**(self.dim + 1) + 1 - i)
              self.child_tuple_insert(level + 1,childLength,0,childAdd,nextpad)
          else:
            self.logfile.write(pad + 'Adding to descendant counts by multiplying %s counts\n' % self.typestr[segtype])
            self.calculate_case_weak(segtype - 1)
            for (childLength,childCount,childType) in self.caseWeak[segtype - 1]:
              childAdd = count * childCount
              self.child_tuple_insert(level + 1,childLength + 2**(self.dim) * (length - 2),childType,childAdd,nextpad)
      return count

  def enumerate(self,level,length):

    self.logfile.write("Calculating counts of all types for length %d segments on level %d\n" % (level, length))
    count = [0] * (self.dim + 2)
    for segtype in range(-1,self.dim+1):
      self.logfile.write("***\nEntering recursion to count %s segments\n" % self.typestr[segtype])
      count[segtype + 1] = self.enumerate_internal(level,length,segtype,' ')
      self.logfile.write("***\nExiting recursion to count %s segments, counted %d\n" % (self.typestr[segtype],count[segtype+1]))
    assert sum(count) == 2**(self.dim*level) + 1 - length
    return count

def pairs(s):
  try:
    level, length = map(int, s.split(','))
    if length <= 0:
      raise argparse.ArgumentTypeError("length %s not allowed: must be positive" % length)
    return (level,length)
  except argparse.ArgumentTypeError:
    raise
  except:
    raise argparse.ArgumentTypeError("segments must be level,length")

if __name__ == '__main__':

  parser = argparse.ArgumentParser(description='Enumerate continuous and discontinuous Morton curve segments')
  parser.add_argument('--dimension', '-d', help='spatial dimension', type=int, default=2)
  parser.add_argument('--segments', '-s', help='list of (level, length) segments to calculate',type=pairs,nargs='*',default=[(0,1)])
  parser.add_argument('--verbose', '-v', help='verbose output (optional filename)',nargs='?', type=argparse.FileType('w'), default=os.devnull)
  parser.add_argument('--pgfplots', '-p', help='pgfplots-style output (incompatible with --verbose)',action='store_true')
  parser.add_argument('--random', '-r', help='add (level, N): add N randomly chosen segment lengths on level',type=pairs,nargs='*',default=[])
  args = parser.parse_args()

  for pair in args.random:
    import random

    level = pair[0]
    N = pair[1]
    total = 2**(args.dimension*level)
    for i in range(0,N):
      j = int(2**(random.random()*args.dimension*level))
      while [level,j] in args.segments:
        j = int(2**(random.random()*args.dimension*level))
      args.segments.append([level,j])


  if not args.verbose:
    args.verbose = sys.stdout

  table = MortonSegmentTable(args.dimension,args.verbose)

  if args.pgfplots:
    print "dimension level length continuous discontinuous contfrac"
    print "# " + " ".join(sys.argv)
  
  for (level,length) in args.segments:
    counts = table.enumerate(level,length)
    if args.pgfplots:
      logstring = "%d %d %d %d %d %f" % (table.dim, level, length, sum(counts[1:]), counts[0],float(sum(counts[1:]))/float(sum(counts)))
    else:
      logstring = "%d-D Morton curve with refinement level %d has %d segments of length %d:\n  %d continuous segments and %d discontinuous segments\n  (%f%%,%f%%)" % (table.dim, level, sum(counts), length, sum(counts[1:]), counts[0],100.*float(sum(counts[1:]))/float(sum(counts)),100.*float(counts[0])/float(sum(counts)))
    print logstring
    if table.logfile is not sys.stdout:
      table.logfile.write("%s\n" % logstring)
