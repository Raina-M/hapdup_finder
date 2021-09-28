# HapDup-Finder
This script aims to find the haplotype-specific duplications that caused by heterozygous variations.
Since it is invetible that certain homologous loci exist divergence for most diploid organisms,
assemblers might report the same regions twice, one of which is duplicated sequences for a haploid assembly.
This script will identify such regions and keep only one sequence of every homologous pair, then integrate it
into primary seqeunces, so the other is alternate sequences for corresponding region. In addition,
primary sequences also include those consensuses (homozygous regions collapse sequnces from pairf chromosomes).

```
sh ./hapdup_finder.sh -h
This script identifies haplotype-specific duplications caused by heterozygous variantions.
The names of duplicated sequences will be indicated in the end. The script can be applied
both at the contig and scaffold level

Syntax: sh ./main.sh -f <FA>  [-n|x|t|o]
options:
-s		Specify when sequence file is scaffolds.
-f FILE		Sequence file that is required to inspect haplotype-specific duplicates.
-n INT          Number of chromosome-level scaffolds. Only used when -s is activated.
-x FLOAT	Upper bound of sequence divergence used in minimap2 alignments.
		Allowed input: 0.05, 0.1 (by default) and 0.2. Default: 0.1.
-t INT          Number of threads used in minimap2. Default: 3.
-o DIR		Output directory. All generated files are put in this directory.
		Default: present working directory.
-h		Print this help and exit.
```

## Dependencies
minimap2
python 2.7+
samtools (not required if index of sequences, i.e., `.fai` file, exists)

## Output
If `-s` is activated, which means the input file is at scaffold level ('N' in the sequences), you will get two files naming:
- primary_seqs.fa
- alternate_seqs.fa

If not, i.e., input sequeces are contigs, then you will have two files naming:
- primary_contigs.fa
- alternate_contigs.fa
