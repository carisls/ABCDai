###############################################################
### Copyright©∏è2024. Caris MPI, Inc. All rights reserved.###
###############################################################

import pandas as pd
import argparse
from pathlib import Path
from itertools import product
from collections import Counter
from math import log


# TODO: Get rid of any fragment ends which have a soft or hard clip

def get_args_parser():
    parser = argparse.ArgumentParser(description='Extract array set of per-fragment metrics')
    parser.add_argument('-i', '--input-pkl-path', type=str, help='The path to the fragment pkl file', required=True)
    parser.add_argument('-o', '--output-dir', type=str,
                        help='The path where the aggregated output feature sets will be saved')
    parser.add_argument('-en', '--extended-name', type=str,
                        help='save the pickle file with the extended name to account for duplicate samples from different flowcells. default: no',default = 'no')
    parser.add_argument('-nr12', '--no-read12', type=str,
                        help='filter cohort to ONLY look at those fragments which do NOT not have a read1 or 2 flag (assuming that means it was collapsed)',
                        default='no')
    parser.add_argument('-pp', '--proper-pair', type=str,
                        help='filter cohort to ONLY proper pairs',
                        default='no')
    # parser.add_argument('-sf', '--short-frags', type=str,
    #                     help='filter cohort to ONLY frags shorter than 150 bp',
    #                     default='no')
    # parser.add_argument('-lf', '--long-frags', type=str,
    #                     help='filter cohort to ONLY frags longer than 160 bp',
    #                     default='no')
    parser.add_argument('-e', '--ends', type=str,
                        help='which molecule ends do we extract motifs for. default = "both", any other input else will give 5',
                        default='both')
    parser.add_argument('-cr', '--clip_remove', type=str,
                        help=' Should we filter out any reads which have a soft/hard clip at the ends of the molecule. default = "yes", any other input else will NOT exclude clipped reads',
                        default='yes')



    return parser.parse_args()


def reverse(dna):
    return ''.join([base for base in dna[::-1]])

def complement(dna):
    complement_map = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    return ''.join([complement_map[base] for base in dna])

def rc(dna):
    return reverse(complement(dna))

def extract_collapsed_motifs(df,ends = 'both'):
    collapsed_reads_forward_strand = df[~df['is_reverse']]
    collapsed_reads_reverse_strand = df[df['is_reverse']]
    if ends == 'both':
        forward_motifs = collapsed_reads_forward_strand['left_motif'].to_list() + collapsed_reads_forward_strand['right_motif'].to_list()
        reverse_motifs = collapsed_reads_reverse_strand['left_motif'].to_list() + collapsed_reads_reverse_strand['right_motif'].to_list()
        reverse_motifs = [rc(motif) for motif in reverse_motifs]
        return forward_motifs+reverse_motifs
    else: # only return 5' ends
        forward_motifs = collapsed_reads_forward_strand['left_motif'].to_list()
        reverse_motifs = collapsed_reads_reverse_strand['right_motif'].to_list()
        reverse_motifs = [rc(motif) for motif in reverse_motifs]
        return forward_motifs+reverse_motifs


def extract_paired_motifs(df, ends='both'):
    if ends == 'both':
        read_1_watson = df[(df['is_read1']) & (~df['is_reverse'])]
        read_2_watson = df[(df['is_read2']) & (df['is_reverse'])]
        read_1_clark = df[(df['is_read1']) & (df['is_reverse'])]
        read_2_clark = df[(df['is_read2']) & (~df['is_reverse'])]

        read_1_watson_motifs = read_1_watson['left_motif'].to_list()
        read_2_watson_motifs = read_2_watson['right_motif'].to_list()

        read_1_clark_motifs = read_1_clark['right_motif'].to_list()
        read_1_clark_motifs = [rc(motif) for motif in read_1_clark_motifs]

        read_2_clark_motifs = read_2_clark['left_motif'].to_list()
        read_2_clark_motifs = [rc(motif) for motif in read_2_clark_motifs]

        return read_1_clark_motifs + read_2_clark_motifs + read_1_watson_motifs + read_2_watson_motifs

    else:  # only take the 5' motifs
        read_1_watson = df[(df['is_read1']) & (~df['is_reverse'])]
        read_1_clark = df[(df['is_read1']) & (df['is_reverse'])]

        read_1_watson_motifs = read_1_watson['left_motif'].to_list()
        read_1_clark_motifs = read_1_clark['right_motif'].to_list()
        read_1_clark_motifs = [rc(motif) for motif in read_1_clark_motifs]

        return read_1_watson_motifs + read_1_clark_motifs


def create_motif_features(list_all_motifs):
    # What are all the possible combinations
    li = ['T', 'C', 'G', 'A']
    combs = [''.join(comb) for comb in product(li, repeat=len(li))]
    # count the motifs we see in our fragments
    all_counts = Counter(list_all_motifs)

    # filter down into a dictionary for only the motifs we want
    motif_dict = {key: all_counts[key] for key in combs}

    # normalize into probabilities
    motif_dict_normalized = {key: float(motif_dict[key]) / sum(motif_dict.values()) for key in motif_dict}

    # calculate the shannon scaled entropy
    entropy = sum([((-1 * item) * log(item)) / (log(256)) for item in motif_dict_normalized.values()])

    # create the output dataframe
    out_df = pd.DataFrame(motif_dict_normalized, index=[0])
    out_df['entropy'] = entropy
    return out_df

if __name__ == '__main__':
    args = get_args_parser()

    # print out the arguments used
    for arg, value in sorted(vars(args).items()):
        print("Argument %s: %r", arg, value)


    input_path = Path(args.input_pkl_path)
    input_df = pd.read_pickle(input_path)
    # filter out any fragment which has a mapq <30 (original implementation), is a duplicate, or is not a primary alignment
    input_df = input_df[(input_df['mapping_quality']>=30)&(input_df['duplicate']==False)&(input_df['is_secondary']==False)].reset_index(drop=True)

    # check for additional parsing
    if args.no_read12 == 'yes':
        input_df = input_df[(input_df['is_read1'] == False) & (input_df['is_read2'] == False)].reset_index(
            drop=True)
    if args.proper_pair == 'yes':
        input_df = input_df[input_df['is_proper_pair']].reset_index(drop=True)
    # DEPRECIATED:
    # if args.short_frags == 'yes':
    #     input_df = input_df[input_df['length'] < 150].reset_index(drop=True)
    # if args.long_frags == 'yes':
    #     input_df = input_df[input_df['length'] > 160].reset_index(drop=True)


    # UMI consensus bam files have a peculiar structure. If a forward/reverse read for a fragment reached a certain condition (such as inter-mate distance and minimum read  quality) they get merged together into a single 'read'. This fragment is aligned to the watson strand and both ends can be extracted for motifs. These fragments can be found simply by checking if 'is_read1' and 'is_read2' are both FALSE. The other condition is when read1/read2 are NOT fused together, in this case we need to be careful with what end we use for motifs since we don't want to extract the inside fragment end. Thus in this case we use the following logic (ONLY for properly paired flagged reads)
    """
    We always want the following order:
    5' MOTIF1 ---------MOTIF2 3'
    Read1 = True + Reverse = False = Watson R1 = 4 LEFTMOST nts
    Read2 = True + Reverse = True = Watson R1 = 4 RIGHTMOST nts
    Read1 = True + Reverse = True = Crick R1 = RC of the 4 RIGHTMOST nts
    Read2 = True + Reverse = False = Crick R2 = RC of the 4 LEFTMOST nts
    """

    # split into the merged collapsed reads (fragments) and the read1/read2 combos for different processing paths

    # filter to the collapsed reads
    collapsed_reads = input_df[(~input_df['is_read1']) & (~input_df['is_read2'])]

    # filter to the non-collapsed reads
    non_collapsed_paired_reads = input_df[
        ((input_df['is_read1']) | (input_df['is_read2'])) & (input_df['is_proper_pair'])]

    # extract collapsed motifs
    collapsed_reads_motifs = extract_collapsed_motifs(collapsed_reads,ends = args.ends)

    # extract paired motifs
    paired_reads_motifs = extract_paired_motifs(non_collapsed_paired_reads, ends= args.ends)

    # combine
    all_motifs = collapsed_reads_motifs + paired_reads_motifs

    output_df = create_motif_features(all_motifs)
    output_df['FNAP'] = args.input_pkl_path

    if args.extended_name == 'no':
        output_path = Path(args.output_dir) / (input_path.stem + '_end_motif_features.pkl')
    elif args.extended_name == 'yes':
        output_path = Path(args.output_dir) / (
                    "_".join(input_path.parts[-4:-1]) + '_' + input_path.stem + 'end_motif_features.pkl')

    output_df.to_pickle(output_path)
    print("Output written to",output_path)
    print("Job Completed")
