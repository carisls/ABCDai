###############################################################
### Copyright©2024. Caris MPI, Inc. All rights reserved.###
###############################################################

import argparse
from pathlib import Path
from time import perf_counter
import subprocess
import yaml

def get_args_parser():
    parser = argparse.ArgumentParser(description='run CNVKit on a tumor sample')
    parser.add_argument('-i', '--input-bam-path', type=str, help='The path to the bam file from which we want to extract depth statistics', required=True)
    parser.add_argument("-c", "--config", type=str, help="Configuration File", required=True)
    parser.add_argument('-ir', '--reference', type=str,
                        help='The path to the reference .cnn file to apply', required=True)
    parser.add_argument('-nj', '--n-jobs', type=int,
                        help='number of jobs. default = 1', default=1)
    parser.add_argument('-o', '--output-dir', type=str,
                        help='The path where the aggregated outputs will be saved, if not provided they are put into ~input-bam-path/cnvkit_outputs/')

    return parser.parse_args()


def run_cnvkit_command(bam_path,reference_path,output_directory,n_jobs):
    start = perf_counter()
    cnv_kit_bash_command = f"{cnvkitpath} batch {bam_path} -r {reference_path} -p {n_jobs} -d {output_directory} --diagram --scatter --segment-method hmm"
    subprocess.call(cnv_kit_bash_command, shell=True)
    end = perf_counter()
    print(end-start)

def run_cnvkit_metrics_command(bam_name,output_directory):
    start = perf_counter()
    cnr_file_path = f"{output_directory}/{bam_name}.cnr"
    cns_file_path = f"{output_directory}/{bam_name}.cns"
    metrics_output_path = f"{output_directory}/{bam_name}_metrics.txt"
    print(f"CNR metrics path: {metrics_output_path}")
    cnv_kit_bash_command = f"{cnvkitpath} metrics {cnr_file_path} -s {cns_file_path}> {metrics_output_path}"
    subprocess.call(cnv_kit_bash_command, shell=True)
    end = perf_counter()
    print(end - start)

def run_cnvkit_segmetrics_command(bam_name,output_directory):
    start = perf_counter()
    cnr_file_path = f"{output_directory}/{bam_name}.cnr"
    cns_file_path = f"{output_directory}/{bam_name}.cns"
    metrics_output_path = f"{output_directory}/{bam_name}_segmetrics.txt"
    print(f"CNR segmetrics path: {metrics_output_path}")
    cnv_kit_bash_command = f"{cnvkitpath} segmetrics {cnr_file_path} -s {cns_file_path} --mean --std --mad --iqr --bivar --mse --ci --pi -o {metrics_output_path}"
    subprocess.call(cnv_kit_bash_command, shell=True)
    end = perf_counter()
    print(end - start)

def run_cnvkit_genometrics_cnr_command(bam_name,output_directory):
    start = perf_counter()
    cnr_file_path = f"{output_directory}/{bam_name}.cnr"
    metrics_output_path = f"{output_directory}/{bam_name}_genemetrics_cnr.txt"
    print(f"CNR genometrics path: {metrics_output_path}")
    cnv_kit_bash_command = f"{cnvkitpath} genemetrics {cnr_file_path} -t 0 -m -0 > {metrics_output_path}"
    subprocess.call(cnv_kit_bash_command, shell=True)
    end = perf_counter()
    print(end - start)


def run_cnvkit_genometrics_cns_command(bam_name, output_directory):
    start = perf_counter()
    cnr_file_path = f"{output_directory}/{bam_name}.cnr"
    cns_file_path = f"{output_directory}/{bam_name}.cns"
    metrics_output_path = f"{output_directory}/{bam_name}_genemetrics_cns.txt"
    print(f"CNS genometrics path: {metrics_output_path}")
    cnv_kit_bash_command = f"{cnvkitpath} genemetrics {cnr_file_path} -s {cns_file_path} -t 0 -m -0 > {metrics_output_path}"
    subprocess.call(cnv_kit_bash_command, shell=True)
    end = perf_counter()
    print(end - start)

if __name__ == '__main__':

    args = get_args_parser()
    
    configPath = args.config
    config = yaml.full_load(open(configPath,'rt'))
    cnvkitpath = config['cnvkitPath']
    
    # if an output directory is not provided create one inside the bam location
    if args.output_dir is None:
        args.output_dir = str(Path(args.input_bam_path).parents[0] / 'depth_outputs')
    else:
        # to remove the last forward slash to be able to fit inside the directory creation
        args.output_dir = str(Path(args.output_dir))
        # print the arguments
    for arg, value in sorted(vars(args).items()):
        print(f"Argument {arg}: {value}")

    # run CNVkit on the new sample. It will create the output folder if it doesn't exist
    run_cnvkit_command(args.input_bam_path, args.reference, args.output_dir, args.n_jobs)

    # now we can run the rest of the tabular reports inside the same output directory
    sample_name = Path(args.input_bam_path).stem
    run_cnvkit_metrics_command(sample_name, args.output_dir)
    run_cnvkit_segmetrics_command(sample_name, args.output_dir)
    run_cnvkit_genometrics_cnr_command(sample_name, args.output_dir)
    run_cnvkit_genometrics_cns_command(sample_name, args.output_dir)

    print ('Job Completed')
