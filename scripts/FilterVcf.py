###############################################################
### Copyright©¸2024. Caris MPI, Inc. All rights reserved.###
###############################################################


import pandas as pd
import argparse
import pathlib
import os
from pathlib import Path

parser = argparse.ArgumentParser(description='Filter VCF file stroing in perfanemnt folder while full VCF and other intermediate files can be archived')
parser.add_argument('--input', dest='inFile', type=pathlib.Path, help='Input VCF File', required = True)
parser.add_argument('--output', dest='outFile', type=pathlib.Path, help='Output CSV File', required = False)
parser.add_argument('--cprv', dest='cprvFile', type=pathlib.Path, help='CSV file with CPRV column', required = False)
parser.add_argument('--exclusion', dest='exclusionFile', type=pathlib.Path, help='CSV file with CPRV column for CPRVs to exclude', required = False)
parser.add_argument('--force', dest='force', action='store_true', help='Force Overwrite', required = False)
parser.add_argument('--verbose', dest='verbose', action='store_true', help='Verbose output', required = False)


args = parser.parse_args()
VCF_FNAP = args.inFile.as_posix()

if args.outFile is None:
    outputPath = VCF_FNAP + '.uber_filtered.csv'
else:
    outputPath = args.outFile.as_posix()

if os.path.exists(outputPath) and not args.force:
    print("Output file exists. Use --force to overwrite. Exiting...")
    exit(-1)

CPRV_to_include_set = {}
if args.cprvFile is not None:
    CPRV_to_include_set = set(pd.read_csv(args.cprvFile.as_posix()).CPRV.unique())

CPRV_to_exclude_set = {}
if args.exclusionFile is not None:
    CPRV_to_exclude_set = set(pd.read_csv(args.exclusionFile.as_posix()).CPRV.unique())

def get_skip_rows(file):
    with open(file) as f:
        counter = 0
        line = f.readline()
        while line.startswith('#'):
            line = f.readline()
            counter += 1
        return(counter-1)
    
def addAlterationType(df):
    df[['CPRV_CHROM', 'CPRV_POS', 'CPRV_REF', 'CPRV_ALT']] = df['CPRV'].str.split('|', expand=True)

    df['RLen'] = df.CPRV_REF.str.len()
    df['ALen'] = df.CPRV_ALT.str.len()
    df['AlterationType'] = ''
    df['AlterationSubType'] = ''
    df['LenDiff'] = abs(df['RLen'] - df['ALen'])
    df['Modulo3'] = df['LenDiff'] % 3
    df.loc[(df.RLen==1) & (df.ALen==1), 'AlterationType'] = 'SNV'
    df.loc[(df.RLen==1) & (df.ALen==1), 'AlterationSubType'] = 'SNV'
    df.loc[(df.RLen>1) | (df.ALen>1), 'AlterationType'] = 'INDEL'
    df.loc[(df.AlterationType=='INDEL') & (df['Modulo3'] != 0), 'AlterationSubType'] = 'FS'
    df.loc[(df.AlterationType=='INDEL') & (df['Modulo3'] == 0), 'AlterationSubType'] = 'INDEL'

    df.drop(columns=['CPRV_CHROM', 'CPRV_POS', 'CPRV_REF', 'CPRV_ALT', 'RLen', 'ALen', 'LenDiff', 'Modulo3'], inplace=True)
    return df

                           
os.makedirs(Path(outputPath).parents[0], exist_ok=True)

                           
try:
    df_ = pd.read_csv(VCF_FNAP, sep='\t', skiprows=get_skip_rows(VCF_FNAP), low_memory=False)
    print(df_.shape[0], 'rows in raw VCF')
    FORMAT_INFO_col = df_.columns[-1]
    dict_list = []
    df_.rename(columns={'#CHROM':'CHROM',FORMAT_INFO_col:'FORMAT_INFO'}, inplace=True)
    idx = 0
    for d in df_.to_dict(orient='records'):
        idx += 1
        if idx % 100_000 == 0:
            print(idx, end='\r')
        #d = {'CHROM':row['#CHROM'],'POS':row['POS'], 'ID':row['ID'], 'REF':row['REF'], 'ALT':row['ALT'],'FILTER':row['FILTER'] ,'FORMAT':row['FORMAT']  ,'FORMAT_INFO':row[FORMAT_INFO_col]   }
        INFO = d.pop('INFO')
        info_list = INFO.split(';')
        for info in info_list:
            try:
                info_key = info.split('=')[0]
                if info_key not in d:
                    d[info_key] = info.split('=')[1]
                else:
                    d['INFO_' + info_key] = info.split('=')[1]
            except:
                if info not in d:
                    d[info] = info
                else:
                    d['INFO_' + info] = info

        FORMAT = d['FORMAT']
        FORMAT_INFO = d['FORMAT_INFO']
        ANN = d['ANN']
        format_dict = dict(zip(FORMAT.split(':'), FORMAT_INFO.split(':')))
        d.update(format_dict)
        ANN_col_list = ['ALLELE', 'EFFECT', 'IMPACT', 'GENE', 'GENEID', 'FEATURE', 'FEATUREID', 'BIOTYPE', 'RANK', 'DC', 'PC', 'ANN_INFO']
        ANN_dict = dict(zip(ANN_col_list, ANN.split('|', maxsplit=11)))
        d.update(ANN_dict)
        dict_list.append(d)
    df = pd.DataFrame(dict_list)
    if df.shape[0] > 0:
        df['CPRV'] = df['CHROM'] + '|' +  df['POS'].astype(str) + '|' +  df['REF'] + '|' +  df['ALT']
        df['TLOD'] = df['TLOD'].astype(float)
        df['MPOS'] = df['MPOS'].astype(int)

       
    if df.shape[0] > 0:
        print('Splitting comma separated feature columns')
        df[['AD_R', 'AD_A']] = df['AD'].str.split(',', expand=True)
        df[['MBQ_R', 'MBQ_A']] = df['MBQ'].str.split(',', expand=True)
        df[['MMQ_R', 'MMQ_A']] = df['MMQ'].str.split(',', expand=True)
        df[['MFRL_R', 'MFRL_A']] = df['MFRL'].str.split(',', expand=True)
        df[['F2R1_R', 'F2R1_A']] = df['F2R1'].str.split(',', expand=True)
        df[['F1R2_R', 'F1R2_A']] = df['F1R2'].str.split(',', expand=True)
        df[['SB_RefForward', 'SB_RefReverse', 'SB_AltForward', 'SB_AltReverse']] = df['SB'].str.split(',', expand=True)
        df = addAlterationType(df)
        print('Converting new feature columns to int')
        for col in [ 'AD_R', 'AD_A', 'MBQ_R', 'MBQ_A', 'MMQ_R', 'MMQ_A', 'MFRL_R', 'MFRL_A', 'F2R1_R', 'F2R1_A', 'F1R2_R', 'F1R2_A', 'SB_RefForward', 'SB_RefReverse', 'SB_AltForward', 'SB_AltReverse']:
            print(col)
            df[col] = df[col].astype(int)
        print('Generating Ratio Columns')
        df['MBQ_Ratio'] = df['MBQ_A'] / df['MBQ_R']
        df['MMQ_Ratio'] = df['MMQ_A'] / df['MMQ_R']
        df['VCF_FNAP'] = VCF_FNAP

        df = df.loc[(df.CPRV.isin(CPRV_to_include_set)) | ((df.TLOD>=2) & (df.MPOS>=6) & (df.MMQ_A >= 20) & (df.MBQ_A >= 30) & (df.MBQ_Ratio >= 0.5))]

        df = df[~df.CPRV.isin(CPRV_to_exclude_set)]
        df.reset_index(inplace=True, drop=True)

except Exception as ERROR:
    df = pd.DataFrame([{'VCF_FNAP':VCF_FNAP, 'ERROR':ERROR}])
    print(ERROR)
if df.shape[0] == 0:
    df = pd.DataFrame(columns='CHROM,POS,ID,REF,ALT,QUAL,FILTER,FORMAT,FORMAT_INFO,AS_FilterStatus,AS_SB_TABLE,DP,ECNT,GERMQ,MBQ,MFRL,MMQ,MPOS,POPAF,TLOD,ANN,GT,AD,AF,F1R2,F2R1,SB,ALLELE,EFFECT,IMPACT,GENE,GENEID,FEATURE,FEATUREID,BIOTYPE,RANK,DC,PC,ANN_INFO,PGT,PID,PS,LOF,RPA,RU,STR,STRQ,NMD,CPRV,AD_R,AD_A,MBQ_R,MBQ_A,MMQ_R,MMQ_A,MFRL_R,MFRL_A,F2R1_R,F2R1_A,F1R2_R,F1R2_A,SB_RefForward,SB_RefReverse,SB_AltForward,SB_AltReverse,AlterationType,AlterationSubType,MBQ_Ratio,MMQ_Ratio,VCF_FNAP'.split(','))

df.to_csv(outputPath, index=False)
print("Output Written to",outputPath)
print("Job Completed")