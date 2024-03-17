###############################################################
### Copyright©¸2024. Caris MPI, Inc. All rights reserved.###
###############################################################

import argparse
import pathlib
import os
import pandas as pd


parser = argparse.ArgumentParser(description='Trim LoFreq VCF files to exclude reads below Variant read count threshold, where variant reads is calculated as DP*AF')
parser.add_argument('--input', dest='inFile', type=pathlib.Path, help='Input VCF File', required = True)
parser.add_argument('--output', dest='outFile', type=pathlib.Path, help='Output VCF File', required = True)
parser.add_argument('--reads', dest='varReadsThreshold', type=float, help='Exclude Variants with VarReads less then this value', required = False, default=1)
parser.add_argument('--sb', dest='sbLT', type=float, help='Exclude Variants with SB greater then or equal to this value', required = False, default = 1000)
parser.add_argument('--hrun', dest='hrunLT', type=float, help='Exclude Variants with HRUN greater then or equal to this value', required = False, default = 2)
parser.add_argument('--dp', dest='dpMin', type=float, help='Exclude Variants with DP less then this value', required = False, default = 3)
parser.add_argument('--cprv', dest='CPRV_Meta', type=pathlib.Path, help='CSV file with CPRVs for the PLP variants', required = False)
parser.add_argument('--force', dest='force', action='store_true', help='Force Overwrite', required = False)
parser.add_argument('--verbose', dest='verbose', action='store_true', help='Write excluded lines to stdout', required = False)

args = parser.parse_args()

print(args)

varReadsMin = args.varReadsThreshold
sbLT = args.sbLT
hrunLT = args.hrunLT
dpMin = args.dpMin


if os.path.exists(args.outFile) and not args.force:
    print("Output file exists. Use --force to overwrite. Exiting...")
    exit(-1)

plpCPRVSet = set(pd.read_csv(args.CPRV_Meta).CPRV.unique()) if args.CPRV_Meta != None else set()
    
exclusionCount = 0
writtenCount = 0
commentCount = 0
with open(args.inFile,'rt') as infile, open(args.outFile,'wt') as outfile:
    for l in infile:
        if l.startswith('#'):
            outfile.write(l)
            commentCount += 1
        else:
            #print(l)
            chrom,pos,snpid,ref,alt,qual,fltr,info = l.strip().split('\t')
            cprv = chrom+'|'+pos+'|'+ref+'|'+alt
            infoparts = info.split(';')
            infodict = {i.split('=')[0]:i.split('=')[1] for i in infoparts if '=' in i}
            varreads = int(infodict['DP'])*float(infodict['AF'])
            sb = float(infodict['SB'])
            if 'HRUN' in infodict:
                hrun = int(infodict['HRUN'])
            else:
                hrun = None
            dp = int(infodict['DP'])            
            if (varreads >= varReadsMin and sb < sbLT and (hrun is None or hrun < hrunLT) and dp >= dpMin) or (cprv in plpCPRVSet):
                outfile.write(l)
                writtenCount += 1
            else:
                exclusionCount +=1
                if args.verbose:
                    print('Excluding', l)

                    
print(exclusionCount,'Variants excluded')
print(commentCount,'Comment lines and',writtenCount,'Variants written to',args.outFile)