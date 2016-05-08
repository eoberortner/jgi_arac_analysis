import argparse
import sys
import os
from collections import defaultdict

## BioPython
from Bio import SeqIO

## pysam for reading SAM/BAM files
import pysam

## for sorting the bins' counts
import operator

parser = argparse.ArgumentParser()
##parser.add_argument('bamfile')
parser.add_argument('reffile',type=argparse.FileType('r'))

args = parser.parse_args()
##if not os.path.exists(args.reffile):
##  print >>sys.stdout, "reference sequences %s does not exist" % args.bamfile


'''
  Setup reference sequence
'''
ref = SeqIO.read(args.reffile,'fasta')
ref_indices = [i for i,v in enumerate(str(ref.seq)) if v in 'acgt']

## positions of mutations
## i.e. 
variable_indices = [i for i,v in enumerate(str(ref.seq)) if v not in 'acgt']

## the reference sequence
refstr = str(ref.seq).upper()
## length of the reference sequence
reflen = len(ref)

def process_read(r):
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
 
    
    bRandom = False

    read_DNA_seq = ''
    counter = 0

    ## ap ... aligned read pair
    ## ap[0] ... position in raw read
    ## ap[1] ... position in reference sequence 
    for ap in r.aligned_pairs:

        if ap[0] == None and ap[1] != None:
            read_DNA_seq += r.seq[ap[1]]
    
            if ap[1] not in variable_indices:
                bRandom = True
    
        elif ap[0] != None and ap[1] != None:
            counter += 1
    
            ref_idx = ap[1]
            read_idx = ap[0]
    
            ## if there's a mutation at a non-indented position, 
            ## then we should store the read in the "random" bin
            if ap[1] not in variable_indices:
                if r.seq[read_idx] != refstr[ref_idx]:
                    bRandom = True 
    
            read_DNA_seq += r.seq[read_idx]

    ## the length of the aligned read must 
    ## match the length of the reference sequence
    assert len(read_DNA_seq) == reflen

    ## translate the read's sequence into AA space
    read_AA_seq = translate_to_AA(read_DNA_seq)

    if not bRandom:
        ## a mutation occurred at an intended position
        return ('COMBINATORIAL', r.qname, read_AA_seq)
  
    ## a mutation occurred at a random position (different 
    ## from the intended positions)
    return ('RANDOM', r.qname, read_AA_seq)


## Codon to AminoAcid map
map = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
       "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
       "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
       "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
       "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
       "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
       "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
       "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
       "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
       "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
       "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
       "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
       "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
       "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
       "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
       "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

def translate_to_AA(str_sequence):
    """
    The translate_to_AA function translates a DNA sequence 
    into a Protein sequence using the global map
    """
    
    codon = '' 
    aa = ''
    for i in range(0, len(str_sequence) - 1):
        codon += str_sequence[i]
        if (i+1)%3 == 0:
            if codon.upper() == 'NNS':
                aa += 'X'
            else:
		if codon.upper() in map:
                	aa += map[codon.upper()]
		else:
			aa += 'N'
            codon = ''

    return aa

 
##########################################################################################

##
## MAIN
##
if __name__ == '__main__':
    
    ## translate the ref. sequence from DNA to Protein sequence
    ref_AA_seq = translate_to_AA(str(ref.seq))

    wildtype = 'PTRHYH'
    
    libraries = ['M4683','M4684', 'M4689', 'M4690','M6607', 'M6608', 'M6609', 'M6610']

    for library in libraries:
        
        print >>sys.stdout, '----- Library: {} -----'.format(library)

        bam_filename = 'aligned/{}.bam'.format(library)
        print >>sys.stdout, 'bam-file: {}'.format(bam_filename)
        
        bamfile = pysam.Samfile(bam_filename,'rb')
        fails = defaultdict(list)
        
        ## counters
        npasses = 0
        nrandoms = 0
        ncombinatorials = 0
        nreads = 0
        
        print >>sys.stdout, "[----- starting binning -----]"

        random_bin = {}
        combinatorial_bin = {}
        
        for r in bamfile:
            pr = process_read(r)
            nreads += 1
        
            if nreads % 10000 == 0:    
                print >>sys.stdout, '{} reads processed'.format(nreads)
        
            if 'RANDOM' in pr[0]:
              nrandoms += 1
        
              if pr[2] not in random_bin:
                  random_bin[pr[2]] = 0 
              
              random_bin[pr[2]] += 1
              
            elif 'COMBINATORIAL' in pr[0]:
              ncombinatorials += 1
              
              if pr[2] not in combinatorial_bin:
                combinatorial_bin[pr[2]] = 0
        
              combinatorial_bin[pr[2]] += 1
        
            else:
              fails[pr[0]].append(pr[1])    
        
        ## print some stats
        print >>sys.stdout, "Total reads:\t%d" % nreads
        print >>sys.stdout, "RANDOMS\t%s" % nrandoms
        print >>sys.stdout, "COMBINATORIALS\t%s" % ncombinatorials
        
        ##
        ## MERGING THE BINS
        ##
        print >>sys.stdout, "[----- merging bins -----]"
        
        merged_bin = {}
        for comb in combinatorial_bin:
            merged_bin[comb] = combinatorial_bin[comb]
        for random in random_bin:   
            if random not in merged_bin:
                merged_bin[random] = random_bin[random] 
            else:
                merged_bin[random] += random_bin[random]
        
        ##
        ## CSV OUTPUT
        ##
        print >>sys.stdout, "[----- generating CSV -----]"

        ## sort the merged bin

        ##sorted_merged = sorted(merged_bin.items(), key=operator.itemgetter(1), reverse=True)
        ##sorted(merged_bin.items(), key=lambda x:x[1])
        ##print sorted_merged
        csv = open('results/{}.csv'.format(library), 'w')
        
        ## write header
        print >>csv, 'Sequence,Count,Random_mutagenesis/Combinatorial, 6-mer'
        
        sorted_list = sorted(merged_bin.items(), key=lambda x:x[1], reverse=True)
        
        ## we only take the first 500 highest counts
        max_entries = 500
        entry_counter = 0
        for entry in sorted_list:

            ## stop after we're reached 500 haplotype counts
            if entry_counter == max_entries:
                break

            haplotype = entry[0]
            count = entry[1]   
            
            ## figure out to what bin the sequence belongs to
            bin = 'R'
            if haplotype in combinatorial_bin:
                bin = 'C'
#                 
            ## construct the 6-mer
            six_mer_AA = ''
            for i in range(0, len(variable_indices)):
                 
                if variable_indices[i] % 3 != 0:
                    continue
                 
                idx = variable_indices[i] / 3
                six_mer_AA += haplotype[idx]

            assert len(six_mer_AA) == len(variable_indices) / 3     
            print >>csv, '{},{},{},{}'.format(haplotype, count, bin, six_mer_AA)
            
            entry_counter += 1
            
        csv.close()
        
        print >>sys.stdout, '--------------------------'


