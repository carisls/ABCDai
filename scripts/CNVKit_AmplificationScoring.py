###############################################################
### Copyright©2024. Caris MPI, Inc. All rights reserved.###
###############################################################


import argparse
from pathlib import Path
import pandas as pd
import os
from glob import glob
import numpy as np

def get_args_parser():
    parser = argparse.ArgumentParser(description='run CNVKit Amplification Calling from a directory with results')
    parser.add_argument('-i', '--input-directory', type=str, help='The directory which has all of the required CNVkit output files', required=True)
    parser.add_argument('-o', '--output-path', type=str,
                        help='The directory where the aggregated output will be saved, if not provided they are put into ~input-directory/AmpCall.csv')
    parser.add_argument('--cnr-log-2', type=float, help='The threshold for the copy number ratio log2 value per gene for a positive call. Default: 0.8', default=0.8)
    parser.add_argument('--cns-low-ci', type=float, help='The threshold for the gene segment ratio log2 lower confidence interval for a positive call. Default: 0.75', default=0.75)
    parser.add_argument('--segment-weight', type=float, help='The threshold for the minimum gene segment weight for a positive call. Default: 100', default=100.0)
    parser.add_argument('--gene-std', type=float,
                        help='The threshold for the minimum all gene bin standard deviation for a positive call. A noisiness filter. Default: 1.0',
                        default=1.0)
    parser.add_argument('--all-genes', type=str, help='experimental. Tag if we want to use ALL genes ("yes") rather than the tissue amp default ("no"). Default: no', default="no")

    return parser.parse_args()

# def Single_Call(row,ci_lo_CNS_threshold,segment_weight_CNS_threshold,):
#     if (row['ci_lo_CNS']>ci_lo_CNS_threshold)&(row['segment_weight_CNS']>segment_weight_CNS_threshold)&(row['log2_CNR']>0.9): #(row['ci_lo_CNS']>.75)&(row['segment_weight_CNS']>200)&(row['log2_CNR']>0.5):
#         return "Amplified"
#     else:
#         return 'Not Amplified'

def genometric_features_cnr(genometric_CNR_file,
                            gene_list=['ASXL1', 'CCND1', 'CDK12', 'ERBB2', 'FGFR2', 'KRAS', 'MET', 'MYC', 'MYCN',
                                       'PDGFRB', 'TACC2', 'TACC3']):
    genometric_df = pd.read_csv(genometric_CNR_file, sep='\t')
    avg_gene_CNR = genometric_df['log2'].median()
    std_gene_CNR = genometric_df['log2'].std()
    output_df = genometric_df[genometric_df['gene'].isin(gene_list)].reset_index(drop=True)


    output_df['FNAP'] = genometric_CNR_file
    output_df['avg_gene_CNR'] = avg_gene_CNR
    output_df['std_gene_CNR'] = std_gene_CNR


    return output_df

def genometric_features_cns(genometric_CNS_file,
                            gene_list=['ASXL1', 'CCND1', 'CDK12', 'ERBB2', 'FGFR2', 'KRAS', 'MET', 'MYC', 'MYCN',
                                       'PDGFRB', 'TACC2', 'TACC3']):
    genometric_df = pd.read_csv(genometric_CNS_file, sep='\t')
    output_df = genometric_df[genometric_df['gene'].isin(gene_list)].reset_index(drop=True)
    output_df['FNAP'] = genometric_CNS_file

    # de-duplicate genes with multiple segments (breakpoints)
    #  Right now what I will do is take the segment with the highest log2.
    #  This should hopefully give us better sensitivity moving forward, but it might be worthwhile to use the larger segment or the one
    # which covers more of the gene
    output_df = output_df.sort_values(by='log2').reset_index(drop=True)
    output_df = output_df.drop_duplicates(subset=['gene'], keep='last').reset_index(
        drop=True)

    return output_df

if __name__ == '__main__':


    # initialize all the genes we want to use. Note, this is from an initial parsing of the 96 batch 65 degree temp cohort.
    # 323 in total but only ~317 should match gene name annotations inside cnvkit
    # updated this to 323 through the following conversions:
    #     if gene == 'CARS':
    #         gene = 'CARS1'
    #     elif gene == 'FGFR1OP':
    #         gene = 'CEP43'
    #     elif gene == 'H3F3A':
    #         gene ='H3-3A'
    #     elif gene == 'H3F3B':
    #         gene ='H3-3B'
    #     elif gene == 'RUNx1T1':
    #         gene ='RUNX1T1'
    #     elif gene == 'WISP3':
    #         gene ='CCN6'
    ALL_Genes = ['AKT2', 'ALK', 'ARID1A', 'AURKB', 'CCND1', 'CCND3', 'CCNE1', 'CD274', 'CDK4', 'CDK6', 'CDKN2A', 'CREBBP',
                 'CRKL', 'EGFR', 'EP300', 'ERBB2', 'EZH2', 'FGF10', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'GATA3',
                 'KDR', 'MAP2K1', 'MCL1', 'MDM2', 'MET', 'MYC', 'NF2', 'NFKBIA', 'NTRK1', 'RB1', 'RICTOR', 'ROS1', 'TOP1',
                 'WT1', 'ABL2', 'ACSL3', 'ACSL6', 'AFF1', 'AFF3', 'AKT3', 'APC', 'ARID2', 'ARNT', 'ASXL1', 'ATF1', 'ATIC',
                 'ATM', 'ATP1A1', 'ATR', 'BAP1', 'BARD1', 'BCL11A', 'BCL2L11', 'BCL3', 'BCL6', 'BCL9', 'BLM', 'BMPR1A',
                 'BRAF', 'BRCA1', 'BRCA2', 'BRIP1', 'CACNA1D', 'CALR', 'CAMTA1', 'CARD11', 'CARS1', 'CASP8', 'CBFB',
                 'CBL', 'CCDC6', 'CCND2', 'CD74', 'CD79A', 'CDC73', 'CDH11', 'CDKN1B', 'CDX2', 'CHEK1', 'CHEK2',
                 'CHIC2', 'CIC', 'CLTCL1', 'CNBP', 'CREB3L2', 'CRTC3', 'CSF1R', 'CTCF', 'CTNNA1', 'CTNNB1', 'CYLD',
                 'DAXX', 'DDR2', 'DDX6', 'DEK', 'DICER1', 'EBF1', 'ECT2L', 'ELK4', 'EPHA3', 'ERBB3', 'ERBB4', 'ERCC2',
                 'ERCC3', 'ERG', 'ESR1', 'ETV1', 'ETV5', 'ETV6', 'EWSR1', 'EXT1', 'EXT2', 'FANCA', 'FANCC', 'FANCD2', 'FANCE',
                 'FANCG', 'FANCL', 'FAS', 'FBXW7', 'FCRL4', 'FGF19', 'FGF23', 'CEP43', 'FGFR4', 'FH', 'FHIT', 'FLCN',
                 'FLI1', 'FLT1', 'FLT4', 'FNBP1', 'FOXA1', 'FOXO1', 'FOXP1', 'FUBP1', 'FUS', 'GID4', 'GMPS', 'GNA13', 'GNAQ',
                 'GNAS', 'GRIN2A', 'H3-3A', 'H3-3B', 'HOOK3', 'HSP90AA1', 'HSP90AB1', 'IDH1', 'IDH2', 'IGF1R', 'IKZF1', 'IL7R',
                 'IRF4', 'ITK', 'JAK1', 'JAK2', 'JAK3', 'JAZF1', 'KEAP1', 'KIAA1549', 'KIF5B', 'KIT', 'KLHL6', 'KMT2A', 'KMT2C',
                 'KMT2D', 'KRAS', 'LCK', 'LHFPL6', 'LIFR', 'LPP', 'LRP1B', 'MAF', 'MAML2', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MDM4',
                 'MDS2', 'MEF2B', 'MEN1', 'MITF', 'MLH1', 'MLLT10', 'MLLT3', 'MRE11', 'MSH2', 'MSH6', 'MSI2', 'MTOR', 'MYB',
                 'MYCN', 'MYD88', 'NCOA2', 'NF1', 'NFIB', 'NFKB2', 'NIN', 'NOTCH2', 'NPM1', 'NR4A3', 'NSD1', 'NTRK2', 'NTRK3',
                 'NUP214', 'NUP93', 'NUP98', 'NUTM1', 'PALB2', 'PAX3', 'PAX5', 'PBRM1', 'PBX1', 'PCM1', 'PDCD1', 'PDCD1LG2',
                 'PDGFRA', 'PDGFRB', 'PER1', 'PIK3CA', 'PIK3R1', 'PIM1', 'PMS2', 'POLE', 'POT1', 'POU2AF1', 'PPARG', 'PRCC',
                 'PRDM1', 'PRKAR1A', 'PRRX1', 'PTCH1', 'PTEN', 'PTPN11', 'RAC1', 'RAD50', 'RAF1', 'RET', 'RNF43', 'RPL22', 'RPN1',
                 'RUNX1', 'RUNX1T1', 'SBDS', 'SDC4', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETBP1', 'SETD2', 'SF3B1', 'SLC34A2',
                 'SMAD2', 'SMAD4', 'SMARCB1', 'SMARCE1', 'SMO', 'SNX29', 'SOX10', 'SPECC1', 'SPEN', 'SRGAP3', 'SRSF2', 'SRSF3',
                 'STAT3', 'STAT5B', 'STIL', 'STK11', 'SUFU', 'SUZ12', 'SYK', 'TAF15', 'TCF7L2', 'TET1', 'TFRC', 'TGFBR2',
                 'TNFAIP3', 'TNFRSF14', 'TP53', 'TPM3', 'TPM4', 'TRIM27', 'TRRAP', 'TSC1', 'TSC2', 'TSHR', 'U2AF1',
                 'USP6', 'VTI1A', 'WDCP', 'CCN6', 'WRN', 'WWTR1', 'XPC', 'YWHAE', 'ZNF217', 'ZNF331', 'ZNF384', 'ZNF521', 'CDK8',
                 'AFDN', 'EZR', 'FLT3', 'MLF1', 'NFE2L2', 'PTPRC', 'HMGA2', 'APOBEC3B', 'AURKA', 'AXL', 'CTLA4', 'HGF', 'PDK1',
                 'PIK3R2', 'VEGFA', 'XPO1', 'VEGFB']

    # Remove the set of genes with multiple genomic positions per the refflat file. NOTE, this does NOT include any gene which has fix/alt patches
    # Example: OR4F29 -  https://genome.ucsc.edu/cgi-bin/hgc?hgsid=1363426453_hKgHZoUrWBVHgsdhGFZn2tbqtOp3&db=hg19&c=chr5&l=112173878&r=180837595&o=180794287&t=180795226&g=refGene&i=NM_001005221
    # Has multiple 100% genomic matches:
    # BROWSER | SIZE IDENTITY CHROMOSOME  STRAND    START     END              QUERY      START  END  TOTAL
    # -----------------------------------------------------------------------------------------------------
    # browser |   939  100.0%          5     + 180794288 180795226          NM_001005221     1   939   939
    # browser |   939  100.0%          1     +    367659    368597          NM_001005221     1   939   939
    # browser |   939  100.0%          1     -    621096    622034          NM_001005221     1   939   939

    # the refflat file can be found here: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ under https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz
    # This is the UCSC genome annotation database. The downloaded file can be found here, for reference, on PHX4:  /gpfs2/home/sklimov/ED_Cases_MolDXPipeline/data/copy_number_data_v2/additional_data/refFlat.txt


    remove_genes = ['AKAP17A', 'ASMT', 'ASMTL', 'ASMTL-AS1', 'BMS1P17', 'BMS1P18', 'BMS1P22', 'CD99', 'CD99P1', 'CRLF2', 'CSF2RA', 'DDX11L16', 'DGCR6', 'DHRSX', 'ENPP7P13', 'FAM138A', 'FAM138C', 'FAM138D', 'FAM138F',
                    'GGT3P', 'GTPBP6', 'H19', 'H3P6', 'IL3RA', 'IL9R', 'LINC00102', 'LINC00106', 'LINC00685', 'LINC01219', 'LINC02564', 'LOC100132062', 'LOC100132287', 'LOC100505874', 'LOC100507412', 'LOC102723376',
                    'LOC102723655', 'LOC102723753', 'LOC102724770', 'LOC102724788', 'LOC103021295', 'LOC105379514', 'MIR10396A', 'MIR10396B', 'MIR12136', 'MIR1244-1', 'MIR1244-2', 'MIR1244-3', 'MIR1244-4', 'MIR1268A',
                    'MIR1302-10', 'MIR1302-11', 'MIR1302-2', 'MIR1302-9', 'MIR3118-1', 'MIR3654', 'MIR3675', 'MIR3690', 'MIR4426', 'MIR4444-1', 'MIR4444-2', 'MIR4454', 'MIR4477A', 'MIR4477B', 'MIR548AI', 'MIR548AQ',
                    'MIR548AU', 'MIR548F3', 'MIR548N', 'MIR5692A1', 'MIR5692A2', 'MIR5701-1', 'MIR5701-2', 'MIR5701-3', 'MIR6089', 'MIR663A', 'MIR6724-1', 'MIR6724-2', 'MIR6724-3', 'MIR6724-4', 'MIR675', 'MIR6859-1',
                    'MIR6859-2', 'MIR6859-3', 'MIR6859-4', 'MIR8069-1', 'MIR8069-2', 'MRPL23', 'NF1P2', 'OR4F16', 'OR4F29', 'OR4F3', 'OR4M2', 'P2RY8', 'PGM5P3-AS1', 'PLCXD1', 'POTEH-AS1', 'PPP2R3B', 'PRODH', 'RNA18SN1',
                    'RNA18SN2', 'RNA18SN3', 'RNA18SN4', 'RNA18SN5', 'RNA28SN1', 'RNA28SN2', 'RNA28SN5', 'RNA45SN1', 'RNA45SN2', 'RNA45SN5', 'RNA5-8SN1', 'RNA5-8SN2', 'RNA5-8SN3', 'RNA5-8SN4', 'RNA5-8SN5', 'RNU1-1',
                    'RNU1-2', 'RNU1-3', 'RNU1-4', 'RNU6-1', 'RNU6-2', 'RNU6-7', 'RNU6-8', 'RNU6-9', 'RNVU1-18', 'RPL21', 'SHOX', 'SLC25A6', 'SNORA105A', 'SNORA105B', 'SNORA59A', 'SNORA59B', 'SNORD131', 'SNORD141A',
                    'SNORD141B', 'SPRY3', 'TP53TG3', 'TP53TG3B', 'TP53TG3C', 'TP53TG3D', 'TP53TG3E', 'TP53TG3F', 'VAMP7', 'WASIR1', 'XGY2', 'ZBED1']


    args = get_args_parser()

    # if an output directory is not provided create one inside the bam location
    if args.output_path is None:
        output_path = Path(str(Path(args.input_directory) / 'AmpCall.csv'))
    else:
        # to remove the last forward slash to be able to fit inside the directory creation
        output_path = Path(args.output_path)
        # print the arguments
    for arg, value in sorted(vars(args).items()):
        print(f"Argument {arg}: {value}")

    # validate that the output directory path exists - if not, create it
    isExist = os.path.exists(output_path.parents[0])
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(output_path.parents[0])
        print("The new directory is created!")
    else:
        print("The directory already exists")

    # Validate that the CNVkit feature files exist
    args.input_directory = str(Path(args.input_directory)) # convert this to a string directory so that the last backslash is gone
    CNR_file_path = glob(args.input_directory+'/*.sorted_genemetrics_cnr.txt')
    assert len(CNR_file_path) == 1
    CNR_file_path = CNR_file_path[0]

    CNS_file_path = glob(args.input_directory+'/*.sorted_genemetrics_cns.txt')
    assert len(CNS_file_path) == 1
    CNS_file_path = CNS_file_path[0]

    if args.all_genes == 'yes':
        # If yes replace the above, tissue amp default gene list, with all the genes annotated by CNVkit
        get_all_genes_df = pd.read_csv(CNR_file_path,sep = '\t')
        ALL_Genes = list(get_all_genes_df['gene'].unique())

    # Remove the duplicated gene positions here
    print (f'Raw Gene list size: {len(ALL_Genes)}')
    ALL_Genes = [gene for gene in ALL_Genes if not any(remove_gene in gene for remove_gene in remove_genes)]
    print(f'Processed Gene list size: {len(ALL_Genes)}')

    results_cnr_df = genometric_features_cnr(CNR_file_path, ALL_Genes)
    results_cns_df = genometric_features_cns(CNS_file_path, ALL_Genes)

    # combine the two results

    # rename first
    columns = results_cnr_df.columns.to_list()
    for column in columns:
        if (column != 'gene')&(column != 'chromosome')&(column != 'start')&(column != 'end'):
            results_cnr_df = results_cnr_df.rename(columns={column: f'{column}_CNR'})


    columns = results_cns_df.columns.to_list()
    for column in columns:
        if column != 'gene':
            results_cns_df = results_cns_df.rename(columns={column: f'{column}_CNS'})

    # merge them
    CNVkit_Combined_df = results_cnr_df[
        ['chromosome','start','end','log2_CNR', 'depth_CNR', 'weight_CNR', 'probes_CNR', 'avg_gene_CNR_CNR', 'std_gene_CNR_CNR',
         'gene','FNAP_CNR']].merge(results_cns_df[
                                  ['log2_CNS', 'depth_CNS', 'weight_CNS', 'ci_hi_CNS', 'ci_lo_CNS', 'probes_CNS',
                                   'segment_weight_CNS', 'segment_probes_CNS', 'gene','FNAP_CNS']],
                              how='left', on=[ 'gene'])

    # Note, the 'chromosome','start','end' (and in turn the CNR gene aggregation through genemetrics) have a bug which causes them to extend the end forward 1 bin.
    # This is due to the masters branch (as of 10/18/2023. Stable branch is still using iloc) switching row slicing from iloc to loc but NOT removing the +1 (many PRs requiring this fix)
    # to the end position (this was required when .iloc was used at it is EXCLUSIVE but using .loc is INCLUSIVE). This occurs in the `.by_gene()` function of the CopyNumArray class.
    # Should be an easy fix but likely has minimal impact (will only 'add' 1 bin to the gene level CNR calculations)
    # You can see a complaint, coming from this issue, filed here: https://github.com/etal/cnvkit/issues/838

    # Assert that there is no duplicates, otherwise throw an error
    assert CNVkit_Combined_df['gene'].value_counts().max() == 1


    # add a column for calling.

    positive_calls = (CNVkit_Combined_df['ci_lo_CNS']>=args.cns_low_ci)&\
                     (CNVkit_Combined_df['segment_weight_CNS']>=args.segment_weight)&\
                     (CNVkit_Combined_df['log2_CNR']>=args.cnr_log_2)& \
                     (CNVkit_Combined_df['std_gene_CNR_CNR'] < args.gene_std) # recommended by Jim

    CNVkit_Combined_df['LP_Amplification'] = np.where(positive_calls,'Amplified','Not Amplified')



    # save the results
    CNVkit_Combined_df.to_csv(output_path,index=False)
    print("Job Completed")
