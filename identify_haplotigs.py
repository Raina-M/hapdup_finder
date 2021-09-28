#!/usr/bin/env python
# coding: utf-8

# # Selection of duplicate contigs from contig assembly
#
# #### Author: Meng Zhang (mzhang@molgen.mpg.de)
# #### Date: Feb 27, 2020
#
# ----------------------------------
#
# This script aims to discern the duplicated contigs in a contig set.
# The duplication occurs when a certain genomic locus has multiple representations,
# e.g., for diploid organisms, one locus can have at most two different nucleotides.
# Thus, the duplicated heterozygous sequences usually exhibit in a contig assembly,
# because contig assemblers fail to collapse honologous contigs due to sequence divergence.
# The duplicated heterozygous sequences will be called haplotigs in the following.

# The duplicated contig pairs are put in **alternate cluster**, which indicates heterozygous sequences,
# and all other contigs are in **primary cluster**. Under the circumstances of purging contigs,
# one further step would be the selection of only one contig from every duplicated pairs in alternate cluster.
# The selected contigs are **primary contigs**, which then will be used for genome assembly.
# The remaining contigs in **alternate cluster** are alternate contigs.
# And surely, all contigs in primary cluster are primary contigs.

# python3
import argparse
import pandas as pd
import numpy as np
from collections import Counter
#import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument('-a','--alninfo', help='File that contains contigs alignment information. String.')
parser.add_argument('-l','--namelst', help='File that contains names of all contig sequences. String.')
parser.add_argument('-m','--mapq', type=int, help='Mapping quality no greater than [mapq] is regarded as low mapping quality. Integer.')
parser.add_argument('-p','--mispercent', type=float, help='Sequence pair whose mismatch rate no geater than [mispercent] is similar sequences. Float.')
parser.add_argument('-o','--outdir', help="Output directory. String.")

args = parser.parse_args()


'''Read contigs-contigs overlap information'''

# A reduced PAF file will be used to identify haplotigs in the contigs assembly.
df = pd.read_csv(args.alninfo, header=None, sep=" ")
df.columns = ['query_contig','query_length','query_start','query_end',
              'target_contig','target_length','target_start','target_end',
              'map_length','mapq']        # Name columns, refering to PAF format interpretation

# Read all contig names
df_contig_names = pd.read_csv(args.namelst, header=None)
NUM_ALL_CTG = len(df_contig_names)


# ########## Step 1: Extract contigs without aligned pairs as primary contigs ##########

'''Seperate unmapped and aligned contigs'''

# extract the names of unmapped contigs
unmapped_contigs = set(df_contig_names.iloc[:,0]) - set(df['query_contig'])
'''
NUM_UMP_CTG = len(unmapped_contigs)                 # number of all unmapped contigs
FRAC_UMP_CTG = round(NUM_UMP_CTG / NUM_ALL_CTG, 4)  # fraction of unmapped contigs
print("# unmapped contigs: {} ({}%)".format(NUM_UMP_CTG, round(FRAC_UMP_CTG*100, 2)))

NUM_ALN_CTG = len(set(df['query_contig']))          # number of all aligned contigs
FRAC_ALN_CTG = round(NUM_ALN_CTG / NUM_ALL_CTG, 4)  # fraction of all aligned contigs

if NUM_ALL_CTG != (NUM_ALN_CTG + NUM_UMP_CTG):
    # Consistency check
    print("Inconsistency of contig numbers! \n")
    print("Please make sure the sum of unmapped and aligned contigs equals to the number of all contigs.")
    quit()

print("# contigs involving in alignment events: {} ({}%)".format(NUM_ALN_CTG, round(FRAC_ALN_CTG*100, 2)))
'''


# ########## Step 2: Aligned contigs with low mapping quality are categorized as primary contigs ##########
#
# Although contigs are able to align to some others, not all of them are real homologous pairs due to false mapping.
# So in this step, alignment pairs with mapping quality lower than `args.mapq` are put in **primary cluster**.

'''Mapping quality control'''

df_high_mapq=df[df['mapq']>=args.mapq]

aln_low_mapq_contigs = set(df['query_contig']) - set(df_high_mapq['query_contig'])
'''
NUM_ALN_HIGH_MPQ_CTG = len(set(df_high_mapq['query_contig'])) # number of aligned contigs with high mapping quality
NUM_ALN_LOW_MPQ_CTG  = len(aln_low_mapq_contigs)              # number of aligned contigs with low mapping quality

if NUM_ALN_CTG != (NUM_ALN_LOW_MPQ_CTG + NUM_ALN_HIGH_MPQ_CTG):
    # Consistency check
    print("Inconsistency of contig numbers!")
    print("Please make sure the sum of aligned contigs with high and low mapping quality equals to the number of all aligned contigs.\n")
    quit()

FRAC_ALN_HIGH_MPQ_CTG = round(NUM_ALN_HIGH_MPQ_CTG / NUM_ALL_CTG, 4)
FRAC_ALN_LOW_MPQ_CTG  = round(NUM_ALN_LOW_MPQ_CTG / NUM_ALL_CTG, 4)

print("Seperate aligned contigs by mapping quality " + str(args.mapq) + ':')
print("  # aligned contigs with high mapping quality: {} ({}%)".format(NUM_ALN_HIGH_MPQ_CTG, round(FRAC_ALN_HIGH_MPQ_CTG*100, 2)))
print("  # aligned contigs with low mapping quality: {} ({}%)".format(NUM_ALN_LOW_MPQ_CTG, round(FRAC_ALN_LOW_MPQ_CTG*100, 2)))
'''


# ########## Step 3: Analysis of aligned contigs with high mapping quality ##########
#
# High mapping quality is not sufficient to say aligned contigs are duplicates,
# because minimap2 split and soft-clip contigs in order to obtain as many alignments
# as possible resulting in a regiment of incomplete contig fragments.
# As yet, what can be explicitly determined group is contigs without aligned pairs
# and those with low mapping qualities, which should be categorized as **primary cluster**.

# Next, the aligned contigs with high confidence will be further analyzed in two situations respectively:
# 1) **complete mapping**: the whole piece of a contig aligns to other contigs
# 2) **split mapping**: contigs are split into multiple fragments, and at least one of them has an alignment

'''Seperate complete mapping and split mapping,
   both with high mapping quality'''

cmplt_map_list, split_map_list =[], []
for tig,occurrence in Counter(df['query_contig']).items():
    if occurrence==1:
        cmplt_map_list.append(tig)
    else:
        split_map_list.append(tig)

cmplt_map_contigs = set(df_high_mapq['query_contig']).intersection(cmplt_map_list)
split_map_contigs = set(df_high_mapq['query_contig']).intersection(split_map_list)

# Seperate complete and split mapping dataframe from high mapping quality dataframe
df_cmplt_mapping = df_high_mapq[df_high_mapq['query_contig'].isin(cmplt_map_contigs)]
df_split_mapping = df_high_mapq[df_high_mapq['query_contig'].isin(split_map_contigs)]

'''
NUM_CMPLT_MAP_CTG = len(cmplt_map_contigs)
NUM_SPLIT_MAP_CTG = len(split_map_contigs)

if NUM_ALN_HIGH_MPQ_CTG != (NUM_CMPLT_MAP_CTG + NUM_SPLIT_MAP_CTG):
    # Consistency check
    print("Inconsistency of contig numbers!")
    print("Please make sure the sum of contigs in complete and split mappings equals to the number of all high-quality aligned contigs.\n")
    quit()

FRAC_CMPLT_MAP_CTG = round(NUM_CMPLT_MAP_CTG / NUM_ALL_CTG, 4)
FRAC_SPLIT_MAP_CTG = round(NUM_SPLIT_MAP_CTG / NUM_ALL_CTG, 4)

print("# contigs in complete mapping: {} ({}%)".format(NUM_CMPLT_MAP_CTG, round(FRAC_CMPLT_MAP_CTG*100, 2)))
print("# contigs in split mapping: {} ({}%)".format(NUM_SPLIT_MAP_CTG, round(FRAC_SPLIT_MAP_CTG*100, 2)))
'''

# Draw a pie chart to present the distribution for each category
'''
class_num = {'class': ['unmapped','low_mapq','high_mapq_complete','high_mapq_split'], 
             'number': [NUM_UMP_CTG, NUM_ALN_LOW_MPQ_CTG, NUM_CMPLT_MAP_CTG, NUM_SPLIT_MAP_CTG]}

df_piechart = pd.DataFrame(data = class_num)

plt.figure(1, figsize=(10,10))

cmap = plt.get_cmap('Spectral')
colors = [cmap(i) for i in np.linspace(0, 1, 8)]

plt.pie(df_piechart['number'], labels=df_piechart['class'], autopct='%1.1f%%', shadow=True, colors=colors)
plt.suptitle('Contigs Alignment Patterns', fontsize=16)

plt.show()
'''


# ### --------------- 3.1 Analysis and categorization of complete mapping ---------------
# First, even though a contig only maps once to another contig, it does not necessarily mean
# that most, if not all of bases in this contig are aligned, because it might be softly clipped
# by minimap2 when doing alignment. In order to select the 'true' duplicate contigs,
# we have to make sure the majority portion of a contig aligns, which is limited by `args.mispercent`.


'''Clipping fraction control of complete mapping:
   (soft clipping length / contig length) greater than args.mispercent is seen as long clipping,
   and contigs with long clipping length are categorized into primary cluster'''

df_long_clippig=df_cmplt_mapping[(abs(df_cmplt_mapping['query_end']-df_cmplt_mapping['query_start'])/df_cmplt_mapping['query_length'])
                                 < (1-args.mispercent)]
df_no_clipping =df_cmplt_mapping[(abs(df_cmplt_mapping['query_end']-df_cmplt_mapping['query_start'])/df_cmplt_mapping['query_length'])
                                 >= (1-args.mispercent)]

long_clip_contigs = set(df_long_clippig['query_contig']).intersection(cmplt_map_contigs)
no_clip_contigs = set(df_no_clipping['query_contig']).intersection(cmplt_map_contigs)

'''
NUM_CTG_LG_CLP = len(long_clip_contigs)  # number of contigs with long soft clipping
NUM_CTG_WO_CLP = len(no_clip_contigs)    # number of contigs without distinct clipping

if NUM_CMPLT_MAP_CTG != (NUM_CTG_LG_CLP + NUM_CTG_WO_CLP):
    # Consistency check
    print("Inconsistency of contig numbers in analysis of complete mapping!")
    quit()

FRAC_CTG_LG_CLP1 = round(NUM_CTG_LG_CLP / NUM_ALL_CTG, 4)
FRAC_CTG_WO_CLP1 = round(NUM_CTG_WO_CLP / NUM_ALL_CTG, 4)
FRAC_CTG_LG_CLP2 = round(NUM_CTG_LG_CLP / NUM_CMPLT_MAP_CTG, 4)
FRAC_CTG_WO_CLP2 = round(NUM_CTG_WO_CLP / NUM_CMPLT_MAP_CTG, 4)

print("when the clipping tolerence is no greater than " + str(args.mispercent)+':')
print(" # contigs with long soft clipping: {} ({}% / {}%)".format(
    NUM_CTG_LG_CLP,round(FRAC_CTG_LG_CLP1*100, 2), round(FRAC_CTG_LG_CLP2*100, 2)))
print(" # contigs without or with short clipping and only aligned once in high mapq: {} ({}% / {}%)".format(
    NUM_CTG_WO_CLP, round(FRAC_CTG_WO_CLP1*100, 2), round(FRAC_CTG_WO_CLP2*100, 2)))
'''

# Contigs with long soft clipping are put into **primary cluster** because the sequence divergences
# between these contigs and their targets are over our preset tolerence, which means we consider them
# as different loci in the genome. Those contigs without or with short clipping and only aligning
# to one target with high mapping quality are extremely similar to at least one of other contigs.
# They are going to be put into **alternate cluster**, and only one contig in each aligned pair will
# be used in the downstream assembly. Others contigs are homologous duplicates for haploid assembly.

# So how to decide which contig in each aligned pair is used to assemble scaffolds, i.e., as a **primary contig**?
# The main idea of the strategy is keeping as many bases as possible. Generally, there are two types of alignments:
# - **a short contig is totally contained by a long contig**: The long contig is surely taken to assemble genome in this scenario.
#   And in our case, all short contigs in such alignment pairs are query sequences, so no contigs with longer targets are categoterized as primary contigs.
#
# - **two contigs are aligned to each other in full length**: This means these two contigs are almost identical, and we keep the one with more bases as primary contig.


df_no_clip=df_cmplt_mapping[df_cmplt_mapping['query_contig'].isin(no_clip_contigs)]

# dataframe of two contigs aligned to each other in full length
df_ = df_no_clip[df_no_clip['query_contig'].isin(np.intersect1d(df_no_clip['query_contig'],
                                                                df_no_clip['target_contig']))]
pairs = [tig for tig, counts in Counter(df_['query_contig'].append(df_['target_contig'])).items() if counts==2]
df_pairs = df_[df_["query_contig"].isin(pairs)]


primary_in_noclip = []
for index, row in df_pairs.iterrows():
    if row['query_length'] >= row['target_length']:
        primary_in_noclip.append(row['query_contig'])
    else:
        primary_in_noclip.append(row['target_contig'])

primary_in_noclip = set(primary_in_noclip)
alternate_in_noclip = no_clip_contigs - primary_in_noclip


# ### --------------- 3.2 Analysis and categorization of split mapping ---------------
#
# Contigs may find confident aligned pairs after splitting into multiple fragments.
# But this does not necessarily mean that those contigs are haplotigs.
# Next, alternate contigs will bw selected through the following four 'filters':
# (N.B. 4 filters are concatenated, i.e., the output of last filter is the input of next filter)
#
# - **Filter 1: occurrence filter**
# Contigs that are mapped more than twice are put into **primary cluster**,
# i.e., the upper limitation for occurrence of split events is 1,
# because the scenarios would become hugely diverse and fiendishly complicated if contigs are split into 3 or more pieces,
# e.g., the existence of repetitive regions in the contig may result in one fragment aligns to another longer fragments from
# the same contig, while other fragments aligned to different contigs, and in this circumstance, this contig is not a haplotig.
# Furthermore, it is rare to see the case that split fragments from one contig align to the same target contig once the number
# of split fragments exceeds 3. And we usually do not regard those contigs mapped to different targets as duplicates.
#
# - **Filter 2: target filter**
# Contigs that are aligned to more than one target contigs are filtered out into **primary cluster**, as is demonstrated above,
#
# - **Filter 3: proximity filter**
# After two previous filters, the survived contigs are split into at most two fragments, and those fragments align to one single target.
# To make sure split fragments are in the vicinity of each other, the gap between those two fragments or the gap between their target regions
# should be no greater than `mapping length*args.mispercent`.
#
# - **Filter 4: clipping filter**
# Contigs whose mis-match fraction is more than `args.mispercent` are put into **primary cluster**,
# because when minimap2 does alignments, some contigs are soft-clipped at sequence ends to decrease sequnece divergence.
# Thus, the clipping sequences should be taken into account to make a better judgement on the similarity of aligned pairs.


"""occurrence filter"""

# contigs splits 2 pieces
contigs_2splits=df['query_contig'].value_counts()[df['query_contig'].value_counts() == 2].index
# contigs occur twice in df_high_mapq
contigs_2highMapq=df_high_mapq['query_contig'].value_counts()[df_high_mapq['query_contig'].value_counts() == 2].index

# contigs that are split into 2 fragments and both aligned with high mapping qualities
df_2splits=df[df['query_contig'].isin(contigs_2highMapq.intersection(contigs_2splits))]


"""target filter + proximity filter + clipping filter"""

contigs_2splits_1target = []

for row in np.arange(0,len(df_2splits),2):
    target_odd = df_2splits['target_contig'][row:row+1].values     # target contig name in odd number rows
    target_even = df_2splits['target_contig'][row+1:row+2].values  # target contig name in even number rows
    map_length_odd = df_2splits['map_length'][row:row+1].values    # length of alignments in odd number rows
    map_length_even = df_2splits['map_length'][row+1:row+2].values # length of alignments in even number rows
    query_length = df_2splits['query_length'][row:row+1].values    # length of query contig
    map_length = map_length_odd + map_length_even

    query_map_locations = [df_2splits['query_start'][row:row+1].values[0],
                           df_2splits['query_end'][row:row+1].values[0],
                           df_2splits['query_start'][row+1:row+2].values[0],
                           df_2splits['query_end'][row+1:row+2].values[0]]
    target_map_locations = [df_2splits['target_start'][row:row+1].values[0],
                           df_2splits['target_end'][row:row+1].values[0],
                           df_2splits['target_start'][row+1:row+2].values[0],
                           df_2splits['target_end'][row+1:row+2].values[0]]

    query_loc1, query_loc2, query_loc3, query_loc4 = sorted(query_map_locations)
    target_loc1, target_loc2, target_loc3, target_loc4 = sorted(target_map_locations)

    target_filter = (target_odd == target_even)
    
    proximity_filter = (query_loc3 - query_loc2 <= map_length*args.mispercent) or (target_loc3 - target_loc2 <= map_length*args.mispercent)
    
    aln_length = query_loc2 - query_loc1 + query_loc4 - query_loc3
    clipping_filter = ((query_loc1 + query_length - query_loc4 ) / query_length < args.mispercent)

    if (target_filter and proximity_filter and clipping_filter):
        contigs_2splits_1target.append(df_2splits['query_contig'][row:row+1].values[0])

'''
NUM_ALT_SLT_CTG = len(contigs_2splits_1target)   # number of alternate contigs in split mapping group
FRAC_ALT_SLT_CTG1 = round(NUM_ALT_SLT_CTG / NUM_ALL_CTG, 4)
FRAC_ALT_SLT_CTG2 = round(NUM_ALT_SLT_CTG / NUM_SPLIT_MAP_CTG, 4)

print("when the divergence tolerence is no greater than " + str(args.mispercent)+':')
print(" # alternate contigs in split mapping group: {} ({}% / {}%)".format(NUM_ALT_SLT_CTG,
        round(FRAC_ALT_SLT_CTG1*100, 2), round(FRAC_ALT_SLT_CTG2*100, 2)))
'''

# ########## Step 4: Extract names of primary and alternate contigs ##########
#
# In previous steps, we have categorized all contigs as either **alternate** or **primary contigs**based on their alignments.
# Now, we will write out their names to helps with subsequent assembly.


alternate_contigs = list(alternate_in_noclip) + contigs_2splits_1target
primary_contigs = list(unmapped_contigs) + list(aln_low_mapq_contigs) + list(long_clip_contigs) + list(primary_in_noclip) + list(split_map_contigs-set(contigs_2splits_1target))

'''
NUM_ALT_CTG=len(alternate_contigs)
NUM_PRM_CTG=len(primary_contigs)

if (NUM_ALT_CTG + NUM_PRM_CTG) != NUM_ALL_CTG:
    # Consistency check
    print("Inconsistency of contig numbers!")
    quit()

print("# alternate contigs:{} ({}%)".format(NUM_ALT_CTG, round(NUM_ALT_CTG/NUM_ALL_CTG*100, 2)))
print("# primary contigs:{} ({}%)".format(NUM_PRM_CTG, round(NUM_PRM_CTG/NUM_ALL_CTG*100, 2)))
'''

# Draw pie chart
'''
contig_num = {'contig': ['primary contig','alternate contig'],
             'number': [NUM_PRM_CTG, NUM_ALT_CTG]}
df_piechart_num = pd.DataFrame(data = contig_num)

contig_length = {'contig': ['primary contig','alternate contig'],
             'length': [2046321771, 51858048]}
df_piechart_len = pd.DataFrame(data = contig_length)

#fig, ax = plt.subplots(figsize=(20, 10), subplot_kw=dict(aspect="equal"))
fig = plt.figure(figsize=(16,8))

cmap = plt.get_cmap('Spectral')
colors = [cmap(i) for i in np.linspace(0.25, 1, 2)]

ax1 = fig.add_subplot(1,2,1)
wedges1, texts1, autotexts1 = ax1.pie(df_piechart_num['number'], #labels=df_piechart['class'],
                                  autopct='%1.2f%%', textprops=dict(color="w"), colors=colors)
ax1.set_title("Contig number")
ax1.legend(wedges1, df_piechart_num['contig'],
          title="Class",
          loc="center left",
          prop={'size': 14},
          bbox_to_anchor=(0.9, -0.3, 0.5, 1))
plt.setp(autotexts1, size=15, weight="bold")

ax2 = fig.add_subplot(1,2,2)
wedges2, texts2, autotexts2 = ax2.pie(df_piechart_len['length'], #labels=df_piechart['class'],
                                  autopct='%1.2f%%', textprops=dict(color="w"), colors=colors)
plt.setp(autotexts2, size=15, weight="bold")
ax2.set_title("Contig length")

plt.show()
'''

alternate_contigs.sort()
primary_contigs.sort()

# Write out the names of all alternate contigs in a text file
f_alt_contigs=open(args.outdir + "alternate_contig_names.txt","w+")
for tig in alternate_contigs:
    f_alt_contigs.write(tig+"\n")

f_alt_contigs.close()


# Write out the names of all primary contigs in a text file
f_prm_contigs=open(args.outdir + "primary_contig_names.txt","w+")
for tig in primary_contigs:
    f_prm_contigs.write(tig+"\n")

f_prm_contigs.close()

