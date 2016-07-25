
#!/usr/bin/env python
import random
import sys
import optparse

random.seed()

parser = optparse.OptionParser()
parser.add_option('-f', '--file', dest="filename", action="store", default=None, help="fasta or fastq file to use")
parser.add_option('-r', '--random', dest="random", action="store", default='all', help="number of random samples to return")
parser.add_option('-m', '--min', dest="min", action="store", default=0, type='long', help="minimum length to include in the set")
parser.add_option('-n', '--max', dest="max", action="store", default=None, type='long', help="maximum length to include in the set")

options, args = parser.parse_args()

#name of the input file (fasta or fastq)
#assumes input file is standard fasta/fastq format
fileName = options.filename
#number of sequences to subsample
numSeq = options.random
increment = 0

#if it's a fasta file
if (fileName.find(".fasta") != -1):
  increment = 2
#else if it's a fastq file
elif (fileName.find(".fastq") != -1):
  increment = 4
#quit if neither
else:
  sys.stdout.write("not a fasta/fastq file\n")
  sys.exit()
 
FILE = open(fileName, 'r')
totalReads = list()
index = 0
total = 0
for line in FILE:
  if(index % increment == 0):
    totalReads.append(index/increment)
    total += 1
  if(index % increment == 1):
    line = line.strip()
    l = len(line)
    min = options.min
    max = options.max
    if max != None and ( l < min or l > max):
      totalReads.pop()
  index += 1
FILE.close()
if numSeq != 'all':
  numSeq = int(numSeq)
  if(len(totalReads) < numSeq):
    sys.stdout.write("You only have "+str(len(totalReads))+" reads!\n")
    sys.exit()
 
random.shuffle(totalReads)
if numSeq != 'all':
  totalReads = totalReads[0: numSeq]
else:
  numSeq = len(totalReads)
totalReads.sort()
 
FILE = open(fileName, 'r')
listIndex = 0
 
for i in range(0, total):
  curRead = ""
  for j in range(0, increment):
    curRead += FILE.readline()
  if (i == totalReads[listIndex]):
    sys.stdout.write(curRead)
    listIndex += 1
    if(listIndex == numSeq):
      break
FILE.close()

