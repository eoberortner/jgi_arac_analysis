import argparse
import sys
import os

parser = argparse.ArgumentParser()
parser.add_argument('varfile',type=argparse.FileType('r'))
parser.add_argument('reffile',type=argparse.FileType('r'))
parser.add_argument('--spcutoff',type=int,default=10)
parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),default=sys.stdout)

args = parser.parse_args()

from Bio import SeqIO
from collections import defaultdict

'''
  Setup reference sequence
'''
ref = SeqIO.read(args.reffile,'fasta')
## ref = SeqIO.read("../reference.fasta",'fasta')
ref_positions = [i for i,v in enumerate(str(ref.seq)) if v in 'acgt']
variable_positions = [i for i,v in enumerate(str(ref.seq)) if v not in 'acgt']

cplx = [dict((ch,0) for ch in 'ATCG') for _ in variable_positions]
species = defaultdict(int)

total = 0

for l in args.varfile:
  seq = l.strip('\n').split('\t')[1]
  species[seq] += 1
  for i,ch in enumerate(seq):
    cplx[i][ch] += 1
  total += 1

print >>args.outfile, "Total number of variants:\t%d" % total
print >>args.outfile, "Number of unique species:\t%d" % len(species)
print >>args.outfile, "Percent of unique species:\t%f" % (float(len(species))/total)
spcutoff = args.spcutoff
sorted_species = sorted([t for t in species.iteritems() if t[1]>=spcutoff],key=lambda x:x[1],reverse=True)
print >>args.outfile, "Number of species > %d:\t%d" % (spcutoff,len(sorted_species)) #sum(x[1]>spcutoff for x in sorted_species)
nread_common = sum(s[1] for s in sorted_species)
print >>args.outfile, "Number of reads in species > %d:\t%d" % (spcutoff,nread_common)
print >>args.outfile, "Percent of reads in species > %d:\t%d" % (spcutoff,(float(nread_common)/total))

print >>args.outfile, "\n[---- Per base summary ---]\n"
for i,vp in enumerate(variable_positions):
  print "Complexity at position %d" % vp
  #t = sum(cplx[i][ch] for ch in 'ATCG')
  for ch in 'ATCG':
    print >>args.outfile, '\t%s\t%f' % (ch,float(cplx[i][ch])/total)
