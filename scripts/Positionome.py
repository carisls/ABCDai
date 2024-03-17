###############################################################
### Copyright©2024. Caris MPI, Inc. All rights reserved.###
###############################################################

import os
import sys
import subprocess
import yaml


import argparse
import pathlib

launchCmd = " ".join(sys.argv)
app=sys.argv[0]
binDir=os.path.dirname(app)
if binDir == "":  binDir = "."

parser = argparse.ArgumentParser('Generate Transcriptome features from BAM file')
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
cbr_bed = config['cbr_bed']
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



config = yaml.full_load(open(configPath,'rt'))
pythonCmd = config['pythonCmd']
exon1_meta = config['exon1_meta']
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


output_path = str(pathlib.Path(args.outputFolder) / (pathlib.Path(args.bam).stem + '_factorome.pkl'))

if not os.path.exists(output_path) or overwrite:
    cmd = f'{pythonCmd} scripts/run_nucleosome.py -i {sorted_bam} -m {exon1_meta} -o {outputFolder}'
    run_command(cmd, 'Positionome')


print("Job Completed")
