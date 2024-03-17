###############################################################
### Copyright©2024. Caris MPI, Inc. All rights reserved.###
###############################################################

import pandas as pd
import argparse
from pathlib import Path
import os
import pysam
import numpy as np

def get_args_parser():
    parser = argparse.ArgumentParser(description='Extract nucleosome features')
    parser.add_argument('-i', '--input-bam-path', type=str, help='The path to the bam file from which we want to extract the exon features', required=True)
    parser.add_argument('-m', '--meta', type=str, help='The path to the csv file containing exon 1 positions', required=True)
    parser.add_argument('-o', '--output-dir', type=str,
                        help='The path where the aggregated outputs will be saved, if not provided they are put into ~input-bam-path/nucleosome/')
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

if __name__ == '__main__':

    args = get_args_parser()
    # if an output directory is not provided create one inside the bam location
    if args.output_dir is None:
        args.output_dir = Path(args.input_bam_path).parents[0] / 'nucleosome'

    # check if the output directory exists, if it doesn't create it
    isExist = os.path.exists(args.output_dir)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(args.output_dir)
        print("The new directory is created!")
    else:
        print("The directory already exists")





    # load the exon1 positions we will analyze
    gene_exon1_positions = pd.read_csv(args.meta)

    # load the bam file
    samfile = pysam.AlignmentFile(args.input_bam_path, "rb")


    results = dict()
    for index, gene_row in gene_exon1_positions.iterrows():
        if gene_row['strand'] == '-':
            region_start = gene_row['start']
            region_end = gene_row['region_exon1_end']
            region_chromosome = f"chr{gene_row['CHROM']}"
        if gene_row['strand'] == '+':
            region_start = gene_row['region_exon1_start']
            region_end = gene_row['end']
            region_chromosome = f"chr{gene_row['CHROM']}"

        region_itr = samfile.fetch(region_chromosome, region_start, region_end)  # CPA1 -> high in panc

        dictlist = []
        for r in region_itr:
            if (not r.is_unmapped) and (r.reference_id != -1) and (not r.is_duplicate) and (
            not r.is_secondary):  # and (r.mapping_quality>=30)
                dictlist.append(get_nucleotide_metrics(r))
        if len(dictlist) == 0:
            results[gene_row['gene']] = {'single_nucleosome_like': np.nan, 'TF_like': np.nan, 'relevant_read_count': 0}
            continue

        fragments_df = pd.DataFrame(dictlist)
        fragments_df = extract_histograms(fragments_df)
        fragments_df['midpoint'] = (fragments_df['start_position'] + (fragments_df['fragment_length'] / 2).astype(int))

        # add the labels for fragment size
        condition = [
            fragments_df['fragment_length'].between(146, 170),
            fragments_df['fragment_length'].between(60, 100)
        ]
        case = [
            'single_nucleosome',
            'TF_like'
        ]
        fragments_df['fragment_type'] = np.select(condition, case, default='other')

        # get counts on either side of the exon start/stop
        nucleosome_like_before_exon1 = sum((fragments_df['fragment_type'] == 'single_nucleosome') & (
                    fragments_df['midpoint'] < gene_row['region_exon1_end']) & (
                                                       fragments_df['midpoint'] >= gene_row['region_exon1_start']))
        nucleosome_like_in_exon1 = sum(
            (fragments_df['fragment_type'] == 'single_nucleosome') & (fragments_df['midpoint'] >= gene_row['start']) & (
                        fragments_df['midpoint'] < gene_row['end']))

        tf_like_before_exon1 = sum(
            (fragments_df['fragment_type'] == 'TF_like') & (fragments_df['midpoint'] < gene_row['region_exon1_end']) & (
                        fragments_df['midpoint'] >= gene_row['region_exon1_start']))
        tf_like_in_exon1 = sum(
            (fragments_df['fragment_type'] == 'TF_like') & (fragments_df['midpoint'] >= gene_row['start']) & (
                        fragments_df['midpoint'] < gene_row['end']))
        # if nucleosome_like_before_exon1>1:
        #     break
        nucleosome_like_score = nucleosome_like_before_exon1 / (nucleosome_like_in_exon1 + 1)
        tf_like_score = tf_like_before_exon1 / (tf_like_in_exon1 + 1)
        results[gene_row['gene']] = {'single_nucleosome_like': nucleosome_like_score, 'TF_like': tf_like_score,
                                     'relevant_read_count': tf_like_before_exon1 + tf_like_in_exon1 + nucleosome_like_before_exon1 + nucleosome_like_in_exon1}

    results_df = pd.DataFrame(results).T.reset_index()
    results_df['FNAP'] = args.input_bam_path

    # save the files
    factorome_path = str(Path(args.output_dir) / (Path(args.input_bam_path).stem + '_factorome.pkl'))

    results_df.to_pickle(factorome_path)