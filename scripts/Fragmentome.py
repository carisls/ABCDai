###############################################################
### Copyright©∏è2024. Caris MPI, Inc. All rights reserved.###
###############################################################
import os
import sys
import pysam
import subprocess
import yaml
from datetime import datetime
import shutil
import pandas as pd
import argparse
import pathlib
import traceback

launchCmd = " ".join(sys.argv)
app=sys.argv[0]
binDir=os.path.dirname(app)
if(binDir == ""):  binDir = "."

parser = argparse.ArgumentParser('Generate Fragmentome features from BAM file using LoFreq and Sentieon TNHaplotyper2 for variant calling')
parser.add_argument("--bam", type=pathlib.Path, help="Input BAM filename and path", required=True)
parser.add_argument("-f", "--outputFolder", type=pathlib.Path, help="Output Folder", required=True)
parser.add_argument("-c", "--config", type=pathlib.Path, help="Configuration File", required=True)
parser.add_argument('--overwrite', action='store_true')
args = parser.parse_args()

sorted_bam = args.bam.as_posix()
outputFolder = args.outputFolder.as_posix()
configPath = args.config.as_posix()
overwrite = args.overwrite

if outputFolder[-1] == '/':
    outputFolder = outputFolder[:-1]


config = yaml.full_load(open(configPath,'rt'))
pythonCmd = config['pythonCmd']
smartbinmeta = config['smartbinmeta']
blacklistcsv = config['fragExclusionRegions']


os.makedirs(outputFolder, exist_ok=True)

def run_command(cmd, label):
    try:
        print(cmd)
        out = subprocess.check_output(cmd, shell=True)
        out = out.decode('ascii')
        print(out)
    except subprocess.CalledProcessError as e:
        print(label+" ERROR: ",e)
        print(e.output)

        
def makePivots(input_csv):
    input_df = pd.read_csv(input_csv)
    input_df['M4'] = input_df[['Mean','Median','Mode']].mean(axis=1)
    pvt_df = input_df.pivot(index=['Sample'],columns='Bin',values=['Median','Mean','Mode','KS_Stat','KS_pValue','WS_Stat','WS_pValue','AD_Statistic','M4'])
    pvt_df.columns = pvt_df.columns.get_level_values(0)+'_'+pvt_df.columns.get_level_values(1).astype(str)
    pvt_df.reset_index(inplace=True)
    return pvt_df

bamFn = args.bam.name
sample = bamFn.split('.')[0].split('_S')[0]

outfile = outputFolder +'/'+ sample+'.smartbin.csv'

if not os.path.exists(outfile) or (overwrite):
    cmd = f'{pythonCmd} scripts/GetSmartBinFragmentDataV3.py -b {sorted_bam} -s {smartbinmeta} -i {blacklistcsv} -o {outfile}'
    run_command(cmd, 'Fragmentome')
    
    
fragmentome_df = makePivots(outfile)
fragmentome_df.to_csv(f'{outputFolder}/{sample}.fragmentome.csv',index=False)
    
print("Job Completed")
