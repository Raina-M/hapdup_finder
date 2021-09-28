#!/usr/bin/ python

# # Identify duplicated small scaffolds
#
# #### Author: Meng Zhang (mzhang@molgen.mpg.de)
# #### Date: Feb 19, 2020
#
# ----------------------------------

# python3
import re
import argparse
import pandas as pd
import numpy as np
from collections import Counter


parser = argparse.ArgumentParser()
parser.add_argument('-a','--alninfo', help='File that contains contigs alignment information. String.')
parser.add_argument('-m','--mapq', type=int, help='Mapping quality no greater than [mapq] is regarded as low mapping quality. Integer.')
parser.add_argument('-p','--mispercent', type=float, help='Sequence pair whose mismatch rate no geater than [mispercent] is similar sequences. Float.')
parser.add_argument('-o','--outdir', help="Output directory. String.")

args = parser.parse_args()


'''Read the alignments information of small scaffolds to chromosomes'''

# A reduced PAF file will be used
df = pd.read_csv(args.alninfo, header=None, sep=" ")
df.columns = ['query_seq','query_length','query_start','query_end', 'mapq']        # Name columns, refering to PAF format interpretation


# ########## Step 1: Mapping quality control ##########

df_high_mapq=df[df['mapq']>=args.mapq]


# ########## Step 2: Analysis of aligned small scaffolds with high mapping quality ##########
#
# High mapping quality is not sufficient to say aligned scaffolds are duplicates,
# because minimap2 split and soft-clip sequences in order to obtain as many alignments
# as possible resulting in a regiment of split scaffold fragments.

# Next, the aligned scaffolds with high confidence will be further analyzed in two situations respectively:
# 1) **complete mapping**: the whole piece of a scaffold aligns to a chromosome.
# 2) **split mapping**: scaffolds are split into multiple fragments, and at least one of them has an alignment.

'''Seperate complete mapping and split mapping,
   both with high mapping quality'''

cmplt_map_list, split_map_list =[], []
for seq,occurrence in Counter(df['query_seq']).items():
    if occurrence==1:
        cmplt_map_list.append(seq)


cmplt_map_sfds = set(df_high_mapq['query_seq']).intersection(cmplt_map_list)

# Seperate complete and split mapping dataframe from high mapping quality dataframe
df_cmplt_mapping = df_high_mapq[df_high_mapq['query_seq'].isin(cmplt_map_sfds)]


# Even though a scaffold only maps once to determined chromosomal regions, it does not necessarily mean
# that most, if not all of bases in this scaffold are aligned, because it might be softly clipped
# by minimap2 when doing alignment. In order to select the 'true' duplicated scaffolds,
# we have to make sure the majority portion of a scaffold aligns, which is limited by `args.mispercent`.

'''Clipping fraction control of complete mapping:
   (soft clipping length / scaffold length) greater than args.mispercent is seen as long clipping,
   and scaffolds with long clipping length are not taken as dupliacted sequences.'''

df_wo_clipping = df_cmplt_mapping[(abs(df_cmplt_mapping['query_end']-df_cmplt_mapping['query_start'])/df_cmplt_mapping['query_length'])
                                  >= (1-args.mispercent)]

sfds_wo_clip = set(df_wo_clipping['query_seq']).intersection(cmplt_map_sfds)


# Small scaffolds with long soft clipping are not duplications because the sequence divergences between
# these scaffolds and their targeted regions are over our preset tolerence, which means we consider them
# from different loci in the genome. Those small scaffolds without or with short clipping and only aligning
# to one target region with high mapping quality are duplications/alternate sequences.


# ########## Step 3: Extract duplicated small scaffolds ##########
#
# In previous steps, we have found all haplotype-specific duplicated scaffolds.
# Now, we write out their names to helps with subsequent assembly.

alternate_sfds = list(sfds_wo_clip)

def get_number(txt):
	match = re.findall(r'[0-9]+$',txt)
	return [int(i) for i in match]
	
alternate_sfds.sort(key=get_number)

# Write out the names of all alternate contigs in a text file
f_alt_seqs=open(args.outdir + "alternate_seq_names.txt","w+")
for seq in alternate_sfds:
    f_alt_seqs.write(seq+"\n")

f_alt_seqs.close()

