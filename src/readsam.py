import argparse
import sys
import os

parser = argparse.ArgumentParser()
parser.add_argument('bamfile')
parser.add_argument('reffile',type=argparse.FileType('r'))
parser.add_argument('--cutoff',type=int,default=10)
parser.add_argument('--stopnum',type=int)
parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),default=sys.stdout)

args = parser.parse_args()
if not os.path.exists(args.bamfile):
  print >>sys.stderr, "Bamfile %s does not exist" % args.bamfile

from Bio import SeqIO
import pysam
from collections import defaultdict

'''
  Setup reference sequence
'''
ref = SeqIO.read(args.reffile,'fasta')
ref_positions = [i for i,v in enumerate(str(ref.seq)) if v in 'acgt']
variable_positions = [i for i,v in enumerate(str(ref.seq)) if v not in 'acgt']
refstr = str(ref.seq).upper()
reflen = len(ref)

def process_read(r,cutoff):
  ''' Processes read '''
  if r.is_unmapped:
    #--- Read is unmapped  
    return ('NM',r.qname)
  if r.is_reverse:
    #--- Read is reversed
    return ('REV',r.qname)
  if r.alen != reflen:
    #--- Read is incomplete: aligned length on reference does not cover reference 
    return ('INC',r.qname)
  ###
  nins = sum(ap[1] is None for ap in r.aligned_pairs)
  ndel = sum(ap[0] is None for ap in r.aligned_pairs)
  alndict = dict((ap[1],ap[0]) for ap in r.aligned_pairs if ap[1] is not None)
  mismatches = sum(refstr[rp]!=r.query[alndict[rp]] for rp in ref_positions if rp in alndict and alndict[rp] is not None)
  if (nins+ndel+mismatches) > cutoff:
    return ('ERR',r.qname)    
  else:
    varstring = ''.join(r.query[alndict[rp]] for rp in variable_positions if rp in alndict and alndict[rp] is not None)
    if len(varstring)==len(variable_positions):
      return ('PASS',r.qname,varstring)
    else:
      return ('INDEL',r.qname)
  assert False, "Should not reach here!"

##########################################################################################

print >>sys.stderr, "[----- starting readsam.py ----- ]"
'''
Do it
'''
bamfile = pysam.Samfile(args.bamfile,'rb')
fails = defaultdict(list)
npasses = 0
if args.stopnum is not None:
  counter = 0
  for r in bamfile:
    counter += 1
    pr = process_read(r,args.cutoff)
    if 'PASS' in pr[0]:
      npasses += 1
      print >>args.outfile, '%s\t%s' % (pr[1],pr[2])
    else:
      fails[pr[0]].append(pr[1])
    if counter >= args.stopnum:
      print >>sys.stderr, "Complete: counter (%d) exceeded stopnum" % counter
      break
else:  
  for r in bamfile:
    pr = process_read(r,args.cutoff)
    if 'PASS' in pr[0]:
      npasses += 1
      print >>args.outfile, '%s\t%s' % (pr[1],pr[2])
    else:
      fails[pr[0]].append(pr[1])    

treads = npasses + sum(len(qnames) for code,qnames in fails.iteritems())
print >>sys.stderr, "Total reads:\t%d" % treads

print >>sys.stderr, "PASS\t%s" % npasses
for code,qnames in fails.iteritems():
  print >>sys.stderr, '%s\t%d' % (code,len(qnames))

print >>sys.stderr, "[----- readsam.py complete ----- ]"

