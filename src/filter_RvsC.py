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

## 
import analysis_RvsC

parser = argparse.ArgumentParser()
##parser.add_argument('bamfile')
parser.add_argument('reffile',type=argparse.FileType('r'))

args = parser.parse_args()
##if not os.path.exists(args.bamfile):
##  print >>sys.stdout, "Bamfile %s does not exist" % args.bamfile


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
                aa += map[codon.upper()]
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
    
    libraries = ['M4683', 'M4684', 'M4689', 'M4690', 'M6607', 'M6608', 'M6609', 'M6610']

    for library in libraries:
        
        print >>sys.stdout, '----- Library: {} -----'.format(library)

        csv_filename = '../results/{}.csv'.format(library)
        
        if not os.path.isfile(csv_filename):
            continue
        
        print >>sys.stdout, 'csv-file: {}'.format(csv_filename)
        
        print >>sys.stdout, "[----- filtering -----]"
        
        ## key ... AA sequence
        ## value ... 
        ##           [count, 'C'/'R']
        bins = {}
        
        with open(csv_filename, 'r') as csvfile:
            
            ## ignore the header row
            next(csvfile)
            
            for line in csvfile:
                
                ## split the line
                elements = line.split(',')
                assert len(elements) == 4

                aa_sequence = elements[0].strip()
                count = int(elements[1].strip())
                type = elements[2].strip()
                haplotype = elements[3].strip()

                bAdd = True                
                if elements[2] == 'R':
                    
                    if haplotype != wildtype:
                        
                        for key in bins:
                            if bins[key][2] == wildtype:
                                ## increment the count of the wild type
                                bins[key][0] += count
                                ## don't add the false positive to the output table
                                bAdd = False

                        ## real random
                        type = 'C'

                if bAdd:
                    bins[aa_sequence] = []
                    bins[aa_sequence].append(int(count))
                    bins[aa_sequence].append(type)     ## R/C
                    bins[aa_sequence].append(haplotype)

        ## serialize the bins to a CSV file
        print >>sys.stdout, "[----- generating filtered CSV -----]"
        
        filtered_csv_file = open('../results/{}.filtered.csv'.format(library), 'w')

        ## write header
        print >>filtered_csv_file, 'Sequence,Count,Random_mutagenesis/Combinatorial,6-mer'

        sorted_bins = sorted(bins.items(), key=lambda x:x[1][0], reverse=True)
        
        for entry in sorted_bins:
  
            print >>filtered_csv_file, '{},{},{},{}'.format(entry[0], entry[1][0], entry[1][1], entry[1][2])
             
        filtered_csv_file.close()
        
        print >>sys.stdout, '--------------------------'

