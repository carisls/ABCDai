###############################################################
### Copyright©2024. Caris MPI, Inc. All rights reserved.###
###############################################################

import pandas as pd
import pysam
from collections import Counter
import numpy as np
from scipy import stats
import argparse
import pathlib


import re
def getFcLnSample(f):
    fc, ln, sample_s, sample = None, None, None, None
    parts = f.split('/')
    for i,p in enumerate(parts):
        if re.match('\d{6}_[A-Z]\d{5}_\d{4}_[A-Z0-9-]+',p) != None or re.match('\d{8}_[A-Z]+\d{5}_\d{4}_[A-Z0-9-]+',p) != None:
            fc = p
            if len(parts) > i+1:
                ln = parts[i+1]
                if len(parts) > i+2:
                    sample_s = parts[i+2]
                    sample = sample_s.split('_S')[0]
            break
    return fc, ln, sample_s, sample



parser = argparse.ArgumentParser('Generate SmartBin Fragment data.')
parser.add_argument('-b',"--bam", type=pathlib.Path, help="Path to BAM file", required=True)
parser.add_argument('-s',"--smartBins", type=pathlib.Path, help="Path to SmartBin definition File", required=True)
parser.add_argument("-o", "--out", type=pathlib.Path, help="Path to desired output csv file", required=True)
parser.add_argument('-i', "--blacklist", type=pathlib.Path, help='Path to file with regions to be excluded', required=False, default=None)
args = parser.parse_args()


bamPath = args.bam.as_posix()
smartBinDef = args.smartBins.as_posix()
blacklistPath = args.blacklist.as_posix()

smartBin_df = pd.read_csv(smartBinDef)
if blacklistPath is not None:
    blacklist_df = pd.read_csv(blacklistPath)
    blacklist_df.sort_values(['CHROM','START']).reset_index(drop=True, inplace=True)
else:
    blacklist_df = pd.DataFrame(columns=['CHROM','START','END','CBR_NAME'])

dictList = []
fc, ln, sample_s, sample = getFcLnSample(bamPath)
bam = pysam.AlignmentFile(bamPath)
for binId in smartBin_df.BinIndex.unique():
    bin_df = smartBin_df[smartBin_df.BinIndex == binId]
    reads = []
    for idx,row in bin_df.iterrows():
        chrom = row['Chrom']
        start = row['StartPos']
        end = row['EndPos']

        if end == -1:
            end = np.inf
        blacklist_df.sort_values(['CHROM','START']).reset_index(drop=True)

        chromBlkList_df = blacklist_df[(blacklist_df.CHROM == chrom) & ((blacklist_df.START.between(start, end) | (blacklist_df.END.between(start,end)) ))]
        intervals = []
        curStart = start
        blkListIdx = 0
        if chromBlkList_df.shape[0] > 0 and chromBlkList_df.START.iloc[0] < start:
            curStart = chromBlkList_df.END.iloc[0]
            blkListIdx = 1
        while (curStart < end) and blkListIdx < chromBlkList_df.shape[0]:
            curEnd = min(end, chromBlkList_df.START.iloc[blkListIdx])
            intervals.append((curStart, curEnd))
            curStart = chromBlkList_df.END.iloc[blkListIdx]
            blkListIdx+=1
        if (curStart < end):
            if end == np.inf:
                end = -1
            intervals.append((curStart, end))

        for s,e in intervals:
            if e == -1:
                reads += bam.fetch(chrom, s)
            else:
                reads += bam.fetch(chrom, s, e+1)

    readLengths=[]
    for r in reads:
        fragLen = len(r.query_sequence)
        if fragLen != 155 and fragLen < 250 and fragLen >= 50:
            readLengths.append(fragLen)
    if len(readLengths) > 10:
        mn = min(readLengths)
        ksResult = stats.kstest(readLengths,stats.norm.cdf)
        wsResult = stats.shapiro(readLengths)
        adResult = stats.anderson(readLengths)
        cntr = Counter(readLengths)
        dictList.append({"Case":sample, 'Flowcell':fc, 'Sample':sample_s, 'Bin':binId, 'ReadCount':len(readLengths), 'Median': np.median(readLengths), 'Mean': round(np.mean(readLengths),2), 'Mode':cntr.most_common(1)[0][0]
                        ,"KS_Stat":ksResult.statistic, "KS_pValue":ksResult.pvalue, 'WS_Stat':wsResult.statistic, 'WS_pValue':wsResult.pvalue,'AD_Statistic':adResult.statistic})
    else:
        dictList.append({"Case":sample, 'Flowcell':fc, 'Sample':sample_s, 'Bin':binId, 'ReadCount':0})
bam.close()

out_df = pd.DataFrame(dictList)

out_df.to_csv(args.out.as_posix(), index=False)
print("output written to", args.out.as_posix())
print("Job Completed")

