araC MiSeq analysis
===================

Here is a description of the protocol used for analyzing the araC MiSeq data

##  Setup reference file

The reference sequence is in the file `ref/reference.fasta`
The reference used is shown below. An 'n' is used in the degenerate positions.

```
>Arac_amplimer_degenerate
atggctgaagcgcaaaatgatnnsctgctgccgggatactcgtttaacgcccatctggtg
gcgggtttannsccgattgaggccaacggttatctcgatttttttatcgacnnsccgctg
ggaatgaaaggttatattctcaatctcaccattcgcggtcagggggtggtgaaaaatcag
ggacgagaattcgtctgccgaccgggtgatattttgctgttcccgccaggagagattnns
cacnnsggtcgtcatccggaggctcgcgaatggtatnnscagtgggtttactttcgtccg
cgcgcctactggcatgaatggctt
```

Index the reference using BWA
```bash
module load bwa/0.7.4
cd ref && bwa index reference.fasta && cd ..
```

##  Setup sample information

Information about samples M4713 and M4714 was found on the RQC webpage (rqc.jgi-psf.org).
The table with individual libraries was cut-and-pasted into the files `sample.M4713.txt`
and `sample.M4714.txt`. The same information is in the tables below:

### sample.M4713.txt

Seq Unit | Individual Library Name | Index Name | Index Sequence | Read Count | Read %
---------|-------------------------|------------|----------------|------------|-------
8209.1.92779.AACAGGTTCGC.fastq.gz | M4686 | None | AACAGGTTCGC | 3,115,112 | 5.88 %
8209.1.92779.AAGTGGACTCT.fastq.gz | M4690 | None | AAGTGGACTCT | 3,724,826 | 7.03 %
8209.1.92779.ACGCTATCTGG.fastq.gz | M4689 | None | ACGCTATCTGG | 4,787,140 | 9.03 %
8209.1.92779.ACTCGTCGGTA.fastq.gz | M4683 | None | ACTCGTCGGTA | 2,888,938 | 5.45 %
8209.1.92779.ACTGTCGAAGC.fastq.gz | M4687 | None | ACTGTCGAAGC | 3,182,442 | 6.00 %
8209.1.92779.AGCGTGTCCAA.fastq.gz | M4684 | None | AGCGTGTCCAA | 3,172,064 | 5.98 %
8209.1.92779.CACATTGTGAG.fastq.gz | M4691 | None | CACATTGTGAG | 3,084,418 | 5.82 %
8209.1.92779.CCTGGATATAG.fastq.gz | M4696 | None | CCTGGATATAG | 3,069,394 | 5.79 %
8209.1.92779.CGATGTCGTCA.fastq.gz | M4692 | None | CGATGTCGTCA | 3,556,850 | 6.71 %
8209.1.92779.CTATGTGAACC.fastq.gz | M4694 | None | CTATGTGAACC | 2,957,946 | 5.58 %
8209.1.92779.CTTGCTGAAGA.fastq.gz | M4697 | None | CTTGCTGAAGA | 2,653,172 | 5.01 %
8209.1.92779.GCACGTTCTAA.fastq.gz | M4693 | None | GCACGTTCTAA | 3,310,564 | 6.25 %
8209.1.92779.GTATCCATGCG.fastq.gz | M4685 | None | GTATCCATGCG | 3,332,354 | 6.29 %
8209.1.92779.TAAGGCCTATC.fastq.gz | M4695 | None | TAAGGCCTATC | 2,955,966 | 5.58 %
8209.1.92779.TTAACGGCTGA.fastq.gz | M4688 | None | TTAACGGCTGA | 2,891,484 | 5.46 %

### sample.M4714.txt

Seq Unit | Individual Library Name | Index Name | Index Sequence | Read Count | Read %
---------|-------------------------|------------|----------------|------------|-------
8215.1.92872.AACAGGTTCGC.fastq.gz | M4701 | None | AACAGGTTCGC | 3,286,480 | 6.92 %
8215.1.92872.AAGTGGACTCT.fastq.gz | M4705 | None | AAGTGGACTCT | 3,007,136 | 6.33 %
8215.1.92872.ACGCTATCTGG.fastq.gz | M4704 | None | ACGCTATCTGG | 2,805,252 | 5.91 %
8215.1.92872.ACTCGTCGGTA.fastq.gz | M4698 | None | ACTCGTCGGTA | 2,472,888 | 5.21 %
8215.1.92872.ACTGTCGAAGC.fastq.gz | M4702 | None | ACTGTCGAAGC | 2,759,514 | 5.81 %
8215.1.92872.AGCGTGTCCAA.fastq.gz | M4699 | None | AGCGTGTCCAA | 2,678,128 | 5.64 %
8215.1.92872.CACATTGTGAG.fastq.gz | M4706 | None | CACATTGTGAG | 3,560,984 | 7.50 %
8215.1.92872.CCTGGATATAG.fastq.gz | M4711 | None | CCTGGATATAG | 2,442,816 | 5.14 %
8215.1.92872.CGATGTCGTCA.fastq.gz | M4707 | None | CGATGTCGTCA | 3,239,198 | 6.82 %
8215.1.92872.CTATGTGAACC.fastq.gz | M4709 | None | CTATGTGAACC | 2,693,864 | 5.67 %
8215.1.92872.CTTGCTGAAGA.fastq.gz | M4712 | None | CTTGCTGAAGA | 4,766,514 | 10.03 %
8215.1.92872.GCACGTTCTAA.fastq.gz | M4708 | None | GCACGTTCTAA | 2,478,634 | 5.22 %
8215.1.92872.GTATCCATGCG.fastq.gz | M4700 | None | GTATCCATGCG | 2,627,380 | 5.53 %
8215.1.92872.TAAGGCCTATC.fastq.gz | M4710 | None | TAAGGCCTATC | 2,436,196 | 5.13 %
8215.1.92872.TTAACGGCTGA.fastq.gz | M4703 | None | TTAACGGCTGA | 2,552,368 | 5.37 %


## Create symlinks to gzipped fastq files

Here are a few lines of code that iterates over the sample files, calculates the path
to the gzipped fastq file, and creates a symlink in the `reads` directory

```bash
# This is where the original sequence files are located
seqloc=/global/dna/dm_archive/sdm/illumina/00

mkdir -p reads

tail -n+2 sample.M4713.txt | while read l; do
  filename=$(echo "$l" | cut -f1)
  libname=$(echo "$l" | cut -f2)
  filepath="$seqloc/${filename:0:2}/${filename:2:2}/$filename"
  [[ -e $filepath ]] && ln -s $filepath reads/$libname.fastq.gz || echo "not found $filename"
done

tail -n+2 sample.M4714.txt | while read l; do
  filename=$(echo "$l" | cut -f1)
  libname=$(echo "$l" | cut -f2)
  filepath="$seqloc/${filename:0:2}/${filename:2:2}/$filename"
  [[ -e $filepath ]] && ln -s $filepath reads/$libname.fastq.gz || echo "not found $filename"
done
```


## Align reads and call haplotypes

The script `analysis.sh` performs all the analysis for each sample.
The steps in the analysis are:
- Reads are merged using FLASH.
- Merged reads are aligned to reference using BWA
- The haplotypes are called using a custom script call readsam.py

Final output is in the `aligned` directory
The `*.varseqs.txt` file contains the haplotype found for each read that passed.

```bash
mkdir -p aligned
tail -n+2 sample.M4713.txt | while read l; do
  libname=$(echo "$l" | cut -f2) && qsub -v samp=$libname analysis.sh
done
tail -n+2 sample.M4714.txt | while read l; do
  libname=$(echo "$l" | cut -f2) && qsub -v samp=$libname analysis.sh
done
```

## Generate summary information

The script `countvars.py` can summarize the information from `*.varseqs.txt`
and provide information such as the complexity per position.

```bash
# Precompute the summary for all samples
for f in aligned/*.varseqs.txt; do python countvars.py $f ref/reference.fasta > ${f%%.*}.varsummary.txt; done
```

