###############################################################
### Copyright©2024. Caris MPI, Inc. All rights reserved.###
###############################################################

import pysam
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import logging
import time
from functools import reduce

"""
Define the functions for the script
"""

def my_bool(s):
    return s != 'False'

def get_args_parser():
    parser = argparse.ArgumentParser(description='Extract array set of per-fragment metrics')
    parser.add_argument('-i', '--input-bam-path', type=str, help= 'The path to the condensed bam file',required=True)
    parser.add_argument('-o', '--output-dir', type=str, help= 'The path where the output serialized pkl file will be saved, The default will be using the same directory as the input file')
    parser.add_argument('-ol', '--output-lvl', type=str,
                        help="The amount of information to pkl for each file. Default: normal",default = 'normal')
    parser.add_argument('-en', '--extended-name', type=str,
                        help='save the pickle file with the extended name to account for duplicate samples from different flowcells. default: No',default = 'no')
    parser.add_argument('-c', '--cbr', type=str, help='The path to the CBR BED file', required=True)
    parser.add_argument('-b', '--blacklist', type=str, help='The path to the blacklist CBR BED file', required=True)


    return parser.parse_args()

def get_nucleotide_metrics(aligned_segment,output_level):
    """
    Get all the additional information which we want to collect for array fragment.
    - start_position -> fragment start position
    - end_position -> fragment end position
    - reference_name -> chromosome
    - GC_percentage -> % of GC nts in the sequence
    - query_quality -> mean of nt phred scores for the sequence
    - mapping_quality -> MAPQ
    - length -> fragment length
    - qc_Fail -> did the fragment QC fail
    - duplicate -> Was the fragment array duplicate
    - left_motif -> last 4 nts at the left end of the fragment
    - right_motif -> last 4 nts at the right end of the fragment
    if extended_info is true then we also extract the following:
    - full_sequence -> full sequence string
    - position_per_sequence_nt -> the corresponding reference positions for the full sequence
    - quality_per_position -> the corresponding phred quality per bp for the full sequence
    """
    sequence = aligned_segment.query_sequence
    # query_qualities = aligned_segment.query_qualities
    if output_level == 'normal':
        output_dict = {'start_position': aligned_segment.reference_start, 'end_position': aligned_segment.reference_end-1,  #-1 since the reference_end parameter points to one past the last aligned residue.
                       'reference_name': aligned_segment.reference_name,
                       'GC_percentage': (sequence.count('G') + sequence.count('C')) / aligned_segment.query_length,  #quence)
                       # 'query_quality': np.mean(query_qualities),
                       'mapping_quality': aligned_segment.mapping_quality, 'length': len(aligned_segment.query_sequence),
                       'qc_Fail': aligned_segment.is_qcfail, 'duplicate': aligned_segment.is_duplicate,
                       'left_motif': sequence[0:4], 'right_motif': sequence[-4:],'Mismatch':aligned_segment.get_tag('NM'),
                       'is_read1': aligned_segment.is_read1,'is_read2': aligned_segment.is_read2, 'is_reverse': aligned_segment.is_reverse,'is_secondary':aligned_segment.is_secondary,
                       'mapped_quality':aligned_segment.mapping_quality,'is_proper_pair':aligned_segment.is_proper_pair,'is_supplementary':aligned_segment.is_supplementary,'is_paired':aligned_segment.is_paired,
                       'read_template_length':aligned_segment.template_length,'read_flag':aligned_segment.flag,'cigarstring':aligned_segment.cigarstring}
    elif output_level == 'extended':
        output_dict = {'start_position': aligned_segment.reference_start, 'end_position': aligned_segment.reference_end-1, #-1 since the reference_end parameter points to one past the last aligned residue.
                       'reference_name': aligned_segment.reference_name,
                       'GC_percentage': (sequence.count('G') + sequence.count('C')) / len(sequence),
                       # 'query_quality': np.mean(query_qualities),
                       'mapping_quality': aligned_segment.mapping_quality, 'length': len(aligned_segment.query_sequence),
                       'qc_Fail': aligned_segment.is_qcfail, 'duplicate': aligned_segment.is_duplicate,
                       'left_motif': sequence[0:4], 'right_motif': sequence[-4:],'full_sequence':sequence,
                       'position_per_sequence_nt':aligned_segment.get_reference_positions(full_length=True),'quality_per_position':query_qualities}
    elif output_level == 'fragment_score':
        output_dict = {'start_position': aligned_segment.reference_start,
                       'reference_name': aligned_segment.reference_name,
                       'GC_percentage': (sequence.count('G') + sequence.count('C')) / aligned_segment.query_length,  #quence)
                       # 'query_quality': np.mean(query_qualities),
                       'mapping_quality': aligned_segment.mapping_quality, 'length': len(aligned_segment.query_sequence),
                       'qc_Fail': aligned_segment.is_qcfail, 'duplicate': aligned_segment.is_duplicate,
                       'is_read1': aligned_segment.is_read1,'is_read2': aligned_segment.is_read2, 'is_reverse': aligned_segment.is_reverse,'is_secondary':aligned_segment.is_secondary,
                       'mapped_quality':aligned_segment.mapping_quality,'is_proper_pair':aligned_segment.is_proper_pair,'is_supplementary':aligned_segment.is_supplementary,'is_paired':aligned_segment.is_paired,
                       'read_template_length':aligned_segment.template_length}
    elif output_level == 'basic':
         output_dict = {'start_position': aligned_segment.reference_start,
                        'reference_name': aligned_segment.reference_name,
                        'length': len(aligned_segment.query_sequence),
                        }

    return output_dict



"""
Run the script
"""


if __name__ == '__main__':

    args = get_args_parser()

    args.input_bam_path = Path(args.input_bam_path)
    if not args.input_bam_path.is_file():
        print(args.input_bam_path)
        print(args.input_bam_path.is_file())
        raise 'Input directory is array mandatory variable!'

    if args.output_dir is None: # if no argument is given then default to the input directory
        args.output_dir = args.input_bam_path.parents[0] # get the directory of the input bam file
    # todo - use os.path.exist and/or assign to var name instead of args.
    else: # if argument is provided here then we will create the directory if it doesn't exist
        args.output_dir = Path(args.output_dir)
        args.output_dir.mkdir(parents=True, exist_ok=True)


    # log the start of the run and all the arguments passed in
    logging.info('#### New Fragment Metric Run ####')
    for arg, value in sorted(vars(args).items()):
        logging.info("Argument %s: %r", arg, value)

    start_time = time.perf_counter()
    # # load the bam file
    bam_file = pysam.AlignmentFile(args.input_bam_path)

    # load both the CBR bed file we want to use and the blacklist CBR bed file we will use
    CBR_bed_df = pd.read_csv(
        args.cbr,
        sep='\t', header=None)
    CBR_bed_df.columns = ['CHROM', 'START', 'STOP', 'GENE_EXON_list']
    CBR_bed_df['CBR'] = CBR_bed_df['CHROM'] + '|' + CBR_bed_df['START'].astype(str) + '|' + CBR_bed_df['STOP'].astype(
        str)
    # omit non autosomal chromosomes
    CBR_bed_df = CBR_bed_df[~CBR_bed_df['CHROM'].isin(['chrX', 'chrY'])].reset_index(drop=True)

    blacklist_CBR_bed_df = pd.read_csv(
        args.blacklist,
        sep='\t', header=None)
    blacklist_CBR_bed_df.columns = ['CHROM', 'START', 'STOP', 'GENE_EXON_list']  # load the blacklist file
    blacklist_CBR_bed_df['CBR'] = blacklist_CBR_bed_df['CHROM'] + '|' + blacklist_CBR_bed_df['START'].astype(
        str) + '|' + blacklist_CBR_bed_df['STOP'].astype(
        str)
    # exclude the CBRs we don't want to use
    CBR_bed_df = CBR_bed_df[~CBR_bed_df['CBR'].isin(blacklist_CBR_bed_df['CBR'].to_list())].reset_index(drop=True)

    # loop through all the remaining CBR reads
    dictlist = []
    for index, row in CBR_bed_df.iterrows():
        region_itr = bam_file.fetch(f"{row['CHROM']}", row['START'] - 1, row['STOP'] + 1)
        for r in region_itr:
            if (not r.is_unmapped) and (r.reference_id != -1) and (not r.is_duplicate) and (not r.is_secondary):
                #dictlist.append({"CHROM":r.reference_name, "POS":r.reference_start, 'FRAG_LEN':len(r.query_sequence), 'ALIGNED_LENGTH':r.query_alignment_length})
                dictlist.append(get_nucleotide_metrics(r,args.output_lvl))
    bam_file.close()
    if args.extended_name == 'no':
        out_path = args.output_dir / (args.input_bam_path.name + '.FragLen.pkl.gz')
    elif args.extended_name == 'yes':
        out_path = args.output_dir / ("_".join(args.input_bam_path.parts[-4:-1]) + '_' + args.input_bam_path.stem + '.FragLen.pkl')
    pd.DataFrame(dictlist).to_pickle(out_path,compression={'method': 'gzip', 'compresslevel': 6})

    for arg, value in sorted(vars(args).items()):
        print("Argument %s: %r", arg, value)

    print("Output written to",out_path)
    print("Job Completed")