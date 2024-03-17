###############################################################
### Copyright©¸2024. Caris MPI, Inc. All rights reserved.###
###############################################################


#!/usr/bin/env python
# coding: utf-8

# In[14]:


import math
import numpy as np
import pandas as pd
import pickle
import pysam
import os
import string
import sys
import time
import argparse
import pathlib
import subprocess


def reverse(dna):
    return ''.join([base for base in dna[::-1]])
def complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    return ''.join([complement[base] for base in dna])
def rc(dna):
    return reverse(complement(dna))

parser = argparse.ArgumentParser('Generate Linkome Telemetry.')
parser.add_argument('-b', "--bam", type=pathlib.Path, help="Input Fastq filename and path", required=True)
parser.add_argument("--clipFastq", type=pathlib.Path, help="Output Clip Fastq filename and path", required=True)
parser.add_argument("--clipBam", type=pathlib.Path, help="Output Clip BAM filename and path", required=True)
parser.add_argument('-o', "--output", type=pathlib.Path, help="Output Telemetry filename and path", required=True)
parser.add_argument("-r", "--reference", type=pathlib.Path, help="Genome Fasta filename and path", required=True)
args = parser.parse_args()


bam_FNAP     = args.bam.as_posix()
FASTQ_FNAP   = args.clipFastq.as_posix()
clipBAM_FNAP = args.clipBam.as_posix()
output_FNAP  = args.output.as_posix()
hg38path = args.reference.as_posix()

minClipLength = 12

samfile = pysam.AlignmentFile(bam_FNAP, "rb")
dict_list = []
startTime = time.time()
readsProcessed_cnt = 0
fastq_unit_list = []
for read in samfile.fetch():
    if read.cigartuples == None:
        continue

    reference_start = read.reference_start
    orig_reverse = '1' if read.is_reverse else '0'
    if (read.cigartuples[0][0] == 4 or read.cigartuples[-1][0] == 4):
        orig_mapping_quality = read.mapping_quality
        ClipSide = ""
        if (read.cigartuples[0][0] == 4 ):
            offset = read.cigartuples[0][1]
            ClipSide = "LEFT"
        else:
            offset = 0
            ClipSide = "RIGHT"
        for cigartuple in read.cigartuples:
            if cigartuple[0]==4:
                clipLen = cigartuple[1]
                clipPos = reference_start+offset-1
                nonClippedLen = len(read.query_sequence) - clipLen
                if clipLen >= minClipLength:
                    if ClipSide=="LEFT":
                        orig_BP = reference_start
                        fastq_l2 = read.query_sequence[:clipLen]
                        fastq_l3 = '+'
                        fastq_l4 = ''.join(map(lambda x: chr( x+33 ), read.query_qualities[:clipLen]))
                    else:
                        orig_BP = reference_start + nonClippedLen 
                        fastq_l2 = read.query_sequence[-clipLen:]
                        fastq_l3 = '+'
                        fastq_l4 = ''.join(map(lambda x: chr( x+33 ), read.query_qualities[-clipLen:]))                                
                    if orig_reverse=='1':
                        fastq_l2 = rc(fastq_l2)
                        fastq_l4 = reverse(fastq_l4)

                    partner_info = str(read.reference_id) + "|" + read.reference_name + "|" + str(reference_start) + "|" + str(nonClippedLen) + "|" + str(orig_BP) + "|" + ClipSide + "|" + orig_reverse + "|" + str(orig_mapping_quality)
                    fastq_l1 = '@'+read.query_name+"|"+partner_info
                    fastq_unit = fastq_l1 + "|X||X|" + fastq_l2 + "|X||X|" + fastq_l3 + "|X||X|" + fastq_l4 
                    fastq_unit_list.append(fastq_unit)
                break
    readsProcessed_cnt += 1
    if readsProcessed_cnt % 10_000_000 == 0:
        print(readsProcessed_cnt)
    if readsProcessed_cnt > 10_000_000_000:
        break


print("elapsed time: %f" % (time.time() - startTime), 'seconds')
if len(set(fastq_unit_list)) > 0:
    print(len(set(fastq_unit_list)),'clips being dumped to FASTq')
    with open(FASTQ_FNAP, 'w') as f:
        for fastq_unit in set(fastq_unit_list):
            for line in fastq_unit.split("|X||X|"):
                f.write(f"{line}\n")


    # #### BWA

    # In[17]:



    output = subprocess.check_output(f"source /gpfs1/share/tools/sentieon/bin/lic.txt; /gpfs1/share/tools/sentieon/bin/sentieon bwa mem -R '@RG\\tID:NTC\\tSM:NTC\\tPL:ILLUMINA' -t 24 -k 12 -C -M /gpfs1/share/hg38/hg38.fa {FASTQ_FNAP} | /gpfs1/share/tools/sentieon/bin/sentieon util sort -o {clipBAM_FNAP} -t 8 --sam2bam -i -", shell=True)
    print(output.decode('ascii'))



    # #### Process Clip BAM

    # In[18]:


    clipSAMfile = pysam.AlignmentFile(clipBAM_FNAP, "rb")
    mapped_cnt = 0
    unmapped_cnt = 0
    clip_read_dict_list = []
    for read in clipSAMfile:
        clip_len = len(read.query_sequence)
        clip_reverse = '1' if read.is_reverse else '0'
        read_query_name_chunks = read.query_name.split('|')
        orig_read_id = read_query_name_chunks[0]
        orig_Gene = read_query_name_chunks[1]
        orig_CHR = read_query_name_chunks[2]
        orig_reference_start = read_query_name_chunks[3]
        orig_nonClippedLen = read_query_name_chunks[4]
        orig_BP = read_query_name_chunks[5]
        orig_clipSide = read_query_name_chunks[6]
        orig_reverse = read_query_name_chunks[7]
        orig_mapping_quality = read_query_name_chunks[8]
        if read.reference_id >-1:
            clip_CHR = read.reference_name
            clip_reference_start = read.reference_start
            clip_mapping_quality = read.mapping_quality
            if (orig_clipSide=="RIGHT"):
                clip_BP = clip_reference_start
            else:
                clip_BP = clip_reference_start 


            if   (orig_reverse=='0') and (orig_clipSide=="RIGHT") and (clip_reverse=='0'):
                clip_BP = clip_reference_start
            elif (orig_reverse=='0') and (orig_clipSide=="LEFT") and (clip_reverse=='0'):
                clip_BP = clip_reference_start + clip_len
            elif (orig_reverse=='1') and (orig_clipSide=="RIGHT") and (clip_reverse=='1'):
                clip_BP = clip_reference_start
            elif (orig_reverse=='1') and (orig_clipSide=="LEFT") and (clip_reverse=='1'):
                clip_BP = clip_reference_start + clip_len
            elif (orig_reverse=='0') and (orig_clipSide=="RIGHT") and (clip_reverse=='1'):
                clip_BP = clip_reference_start + clip_len
            elif (orig_reverse=='1') and (orig_clipSide=="LEFT") and (clip_reverse=='0'):
                clip_BP = clip_reference_start
            elif (orig_reverse=='1') and (orig_clipSide=="RIGHT") and (clip_reverse=='0'):
                clip_BP = clip_reference_start + clip_len
            elif (orig_reverse=='0') and (orig_clipSide=="LEFT") and (clip_reverse=='1'):
                clip_BP = clip_reference_start


        else:
            clip_CHR = 'Unmapped'
            clip_reference_start = None
            clip_BP = None
            clip_mapping_quality = None
        d = {'orig_read_id':orig_read_id,
             'orig_Gene':orig_Gene,
             'orig_CHR':orig_CHR,
             'orig_reference_start':orig_reference_start,
             'orig_referenceReads':'TODO',
             'orig_nonClippedLen':orig_nonClippedLen,
             'orig_BP':orig_BP,
             'orig_clipSide':orig_clipSide,
             'orig_mapping_quality':orig_mapping_quality,
             'orig_reverse':orig_reverse,
             'clip_CHR':clip_CHR,
             'clip_reference_start':clip_reference_start,
             'clip_BP':clip_BP,
             'clip_len':clip_len,
             'clip_seq':read.query_sequence,
             'clip_reverse': clip_reverse,
             'clip_mapping_quality':clip_mapping_quality,
        }
        clip_read_dict_list.append(d)
    clip_df = pd.DataFrame(clip_read_dict_list)
    clip_df['orig_CHR_BP'] = clip_df['orig_CHR'] + "|" + clip_df['orig_BP']
    clip_df['clip_Gene_source'] = ''
    clip_df.sort_values(by=['orig_CHR_BP', 'clip_len'], ascending=[True, False], inplace=True)
    clip_df.reset_index(inplace=True, drop=True)


# In[19]:


    print(clip_df.shape)
    print(clip_df.loc[(clip_df.clip_mapping_quality>0)].shape)

    rel_df = clip_df.loc[(clip_df.clip_mapping_quality>10) & (clip_df.orig_CHR != clip_df.clip_CHR) & (~clip_df.orig_CHR.str.contains('alt'))& (~clip_df.clip_CHR.str.contains('alt'))& (~clip_df.orig_CHR.str.contains('random'))& (~clip_df.clip_CHR.str.contains('random'))].copy()
    print(rel_df.shape)

    rel_df['clip_BP'] = rel_df['clip_BP'].astype(int)
    rel_df['clip_CHR_BP'] = rel_df['clip_CHR']  + '|' + rel_df['clip_BP'].astype(str)
    rel_df['Fusion'] = rel_df['orig_CHR_BP']  + '::' + rel_df['clip_CHR_BP']
    fusion_df = rel_df.groupby(['Fusion', 'orig_CHR', 'orig_BP', 'clip_CHR', 'clip_BP']).agg({'orig_read_id':'count','orig_nonClippedLen':'median', 'clip_len':'median', 'orig_mapping_quality':'median', 'clip_mapping_quality':'median'}).reset_index()
    print(fusion_df.shape)



    # #### Cleanup and Persist Output

    # In[ ]:


    os.remove(clipBAM_FNAP)
    os.remove(clipBAM_FNAP+'.bai')
    os.remove(FASTQ_FNAP)
else:
    print('No Clips, writing empty dataframe')
    fusion_df = pd.DataFrame(columns=['Fusion', 'orig_CHR', 'orig_BP', 'clip_CHR', 'clip_BP', 'orig_read_id', 'orig_nonClippedLen', 'clip_len', 'orig_mapping_quality', 'clip_mapping_quality'])
fusion_df.to_csv(output_FNAP, index=False)
print("Output written to", output_FNAP)
print("Job Completed")

