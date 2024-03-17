###############################################################
### Copyright©¸2024. Caris MPI, Inc. All rights reserved.###
###############################################################

import pandas as pd
import argparse
from pathlib import Path
import os
import pysam
import numpy as np
import scipy
import time



def get_args_parser():
    parser = argparse.ArgumentParser(description='Extract nucleosome features')
    parser.add_argument('-i', '--input-bam-path', type=str, help='The path to the bam file from which we want to extract the trace features', required=True)
    parser.add_argument('-b', '--bed', type=str,
                        help='The path to CBR bed file')
    parser.add_argument('-o', '--output-dir', type=str,
                        help='The path where the aggregated outputs will be saved, if not provided they are put into ~input-bam-path/stack_outputs_v3/')
    return parser.parse_args()


def get_nucleotide_metrics(aligned_segment):

    output_dict = {'start_position': aligned_segment.reference_start,'end_position': aligned_segment.reference_end-1,'length': aligned_segment.alen,'is_read1': aligned_segment.is_read1,'is_read2': aligned_segment.is_read2,
                    'is_proper_pair':aligned_segment.is_proper_pair,'is_supplementary':aligned_segment.is_supplementary,'is_paired':aligned_segment.is_paired,'read_template_length':aligned_segment.template_length}

    return output_dict

def extract_histograms(input_df):
    #     try:

    collapsed_read_location = (~input_df['is_read1']) & (~input_df['is_read2'])
    non_collapsed_paired_reads_location = ((input_df['is_read1']) | (input_df['is_read2'])) & (
        input_df['is_proper_pair']) & (input_df['read_template_length'] > 0)

    input_df['fragment_length'] = 0
    input_df.loc[collapsed_read_location, 'fragment_length'] = input_df.loc[collapsed_read_location, 'length']
    input_df.loc[non_collapsed_paired_reads_location, 'fragment_length'] = input_df.loc[
        non_collapsed_paired_reads_location, 'read_template_length']

    input_df = input_df[collapsed_read_location | non_collapsed_paired_reads_location].reset_index(drop=True)
    return input_df

def sliding_window_max(binary_vector, window_size= 3 ):
    max_values = np.array([np.max(binary_vector[i:i+window_size]) for i in range(0, len(binary_vector), window_size)])
    return max_values


def get_simple_variation(input_df, column_use, window_size=3):
    unique_values = input_df[column_use].unique()
    # Create a binary vector indicating the presence of each value
    binary_vector = np.isin(np.arange(np.min(unique_values), np.max(unique_values) + 1), unique_values).astype(int)
    unique_bins = sliding_window_max(binary_vector, window_size=window_size)
    total_fragments = input_df.shape[0]
    return sum(unique_bins) / total_fragments


# 80 to 250 for frag length
def get_shannon_entropy(input_df, column_use):
    counts_per_category = input_df[column_use].value_counts()
    proportions_per_category = (counts_per_category / input_df.shape[0])
    entropy = scipy.stats.entropy(proportions_per_category, base=2)

    return entropy


if __name__ == '__main__':

    args = get_args_parser()
    # if an output directory is not provided create one inside the bam location
    if args.output_dir is None:
        args.output_dir = Path(args.input_bam_path).parents[0] / 'stack_outputs_v3'


    # print the arguments
    for arg, value in sorted(vars(args).items()):
        print(f"Argument {arg}: {value}")



    # check if the output directory exists, if it doesn't create it
    isExist = os.path.exists(args.output_dir)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(args.output_dir)
        print("The new directory is created!")
    else:
        print ("The directory already exists")

    # load the bed file which has our positional locations
    CBR_bed_df = pd.read_csv(args.bed, sep='\t', header=None)
    CBR_bed_df.columns = ['CHROM', 'START', 'STOP', 'GENE_EXON_list']
    CBR_bed_df['CBR'] = CBR_bed_df['CHROM'] + '|' + CBR_bed_df['START'].astype(str) + '|' + CBR_bed_df['STOP'].astype(
        str)
    # omit non autosomal chromosomes
    CBR_bed_df = CBR_bed_df[~CBR_bed_df['CHROM'].isin(['chrX', 'chrY'])].reset_index(drop=True)

    # Run the analysis

    sam_file = pysam.AlignmentFile(args.input_bam_path, "rb")

    start = time.time()
    full_list = []
    for index, row in CBR_bed_df.iterrows():
        region_itr = sam_file.fetch(f"{row['CHROM']}", row['START'] - 1, row['STOP'] + 1)

        dictlist = [] # for reads
        output_list = [] # for the CBR feature dictionaries
        for r in region_itr:
            if (not r.is_unmapped) and (r.reference_id != -1) and (not r.is_duplicate) and (
            not r.is_secondary):  # and (r.mapping_quality>=30) (r.reference_start>=row['START']) and (r.reference_end<row['STOP']) and
                dictlist.append(get_nucleotide_metrics(r))
        if len(dictlist) >= 10: # assume we have to have at least 10 reads to perform the requested calculations
            tmp_df = pd.DataFrame(dictlist)
            tmp_fixed_df = extract_histograms(tmp_df)
            if tmp_fixed_df.shape[0] < 10:  # have at least 10 reads with a proper length
                continue
            CBR_value = row['CBR']
            # get the simple features
            start_position_ratio_value = get_simple_variation(tmp_fixed_df, 'start_position')
            output_list.append({'val': start_position_ratio_value,
                                'feature': 'start_position_ratio', 'CBR': CBR_value})

            frag_size_ratio_value = get_simple_variation(tmp_fixed_df, 'fragment_length')
            output_list.append({'val': frag_size_ratio_value,
                                'feature': 'frag_size_ratio', 'CBR': CBR_value})

            # get the entropy features
            entropy_position_value = get_shannon_entropy(tmp_fixed_df, 'start_position')
            output_list.append({'val': entropy_position_value,
                                'feature': 'start_position_entropy', 'CBR': CBR_value})

            entropy_fragment_value = get_shannon_entropy(tmp_fixed_df[tmp_fixed_df['fragment_length'].between(80, 250)],
                                                           'fragment_length')  # 80 to 250 for frag length
            output_list.append({'val': entropy_fragment_value,
                                'feature': 'frag_size_entropy', 'CBR': CBR_value})

            # calculate the simple features multiple
            output_list.append({'val': frag_size_ratio_value * start_position_ratio_value,
                                'feature': 'multiplied_ratio', 'CBR': CBR_value})

            # get the number of fragments (to help filter)
            output_list.append({'val': tmp_fixed_df.shape[0],
                                'feature': 'fragment_count', 'CBR': CBR_value})

            # store as a dictionary
            full_list.extend(output_list)


    end = time.time()
    print('full CBR feature extraction took :',end - start, 'seconds')

    # create the tall dataframe
    output_df = pd.DataFrame(full_list)

    # save the dataframe as a pickle
    output_path = str(Path(args.output_dir) / (Path(args.input_bam_path).stem + '_traceome_v4.pkl'))

    output_df.to_pickle(output_path)
