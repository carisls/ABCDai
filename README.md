Copyright©2024. Caris MPI, Inc. All rights reserved.

Code used for "Validation of an AI-enabled exome/transcriptome liquid biopsy platform for early detection, MRD, and therapy selection for solid tumors"   
For additional information, please contact [Dr.David Spetzler](mailto:dspetzler@carisls.com).


# Running Pillar Extraction Code
The 8 pillars each have their individual feature extraction `.py` scripts which require the following arguments:
- A path to the config file
- A path to the .bam file
- A path to the output directory where extracted data will be stored
This will output the raw feature files and allow developers to filter/ETL accordingly.


The `.py` files for each pillar can be found in the following location:
- scripts/Motifome.py
- scripts/Copyome.py
- scripts/Positionome.py
- scripts/Transcriptionome.py
- scripts/Fragmentome.py
- scripts/Fusionome.py
- scripts/DNA_Linkome.py
- scripts/Mutationome.py

# Config file meta data 
The config file points to required meta data to process the code. This metadata is dependent on individual lab/sequencing considerations including WES vs. WGS, bait boosting strategy, variant universe, specific exlusions, etc. To faciliate running the code we elaborate on the required format for the files which the config path points towards: <br>
- **pythonCmd** <br>
Path to the python interpreter *('~/miniconda3/envs/ud_pipeline/bin/python')*.<br>
- **sentiEnv** <br>
The shell command to create enviornment variables required for sentieon *('export SENTIEON_LICENSE=IP:PORT')* . The script can be modified to use the Broads Mutect2 (https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2)<br>
- **sentiPath** <br>
Path to the sentieon directory. The script can be modified to use the Broads Mutect2 (https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2)<br>
- **cnvkitPath** <br>
Path to the main cnvkit python script. *('~/miniconda3/envs/ud_pipeline/bin/cnvkit.py')*.<br>
- **threads** <br>
Number of threads to use for `TNHaplotyper2`<br>
- **hg38_RefSeq_file** <br>
Path to the Hg38 FASTA file i.e.,(*http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz*)<br>
- **gatk** <br> 
Path to the Broad Genome Analysis Toolkit (GATK)  *(https://gatk.broadinstitute.org/hc/en-us)*<br>
- **snpeff** <br>
Path to the SnpEff toolbox *(https://pcingola.github.io/SnpEff/)* <br>
- **java** <br>
Path to the Java Virtual Machine (JVM) executable <br>
- **targets_padded**  <br>
Path to a `.bed` file which has the genomic intervals used for mutation calling (example below):  <br>
>The header row is for explanatory purposes
>
| Chromosome | Start   | End     |
|------------|---------|---------|
| chr1       | 10000   | 20000   |
| chr2       | 15000   | 25000   |
| chr3       | 20000   | 30000   |

- **seq_dict**  <br>
Path to the sequence dictionary for your reference sequence. This can be generated using:
  - java -jar picard.jar CreateSequenceDictionary

- **cprvUniverse** <br>
Path to a `.csv` file which contains the universe of all mutations seen in your dataset (example below):<br>

| gene | proteinChange | ngs_chromosome | ngs_position | ref_allele | variant_allele | CPRV                 | occurrences |
|------|---------------|----------------|--------------|------------|----------------|----------------------|-------------|
| EGFR | S229C         | chr7           | 55152602     | A          | T              | chr7&#124;55152602&#124;A&#124;T | 1000        |
| BRAF | V600K         | chr7           | 140753335    | CAC        | CTT            | chr7&#124;140753335&#124;CAC&#124;CTT | 999         |
| TP53 | K382fs        | chr17          | 7669644      | GT         | G              | chr17&#124;7669644&#124;GT&#124;G | 998         |

- **exon1_meta** <br>
Path to a `.csv` file that is an exon1 filtered MANE Select gff file with an additional pair of columns (`region_exon1_start` and `region_exon1_end`) which represent the baiting you have avaliable for that exon (example below):<br>

| source     | ID                    | Parent               | gbkey | tag          | seqid       | TI             | TI_Root       | TI_suffix    | type | gene  | product                                                     | exonNumber | CHROM | start  | end    | score | strand | phase | Ensembl          | GeneID | GenBank        | HGNC  | MIM   | IMGT/GENE-DB | length | region_exon1_start | region_exon1_end |
|------------|-----------------------|----------------------|-------|--------------|-------------|----------------|---------------|--------------|------|-------|-------------------------------------------------------------|------------|-------|--------|--------|-------|--------|-------|------------------|--------|----------------|-------|-------|--------------|--------|--------------------|------------------|
| BestRefSeq | exon-NM_001005484.2-1 | rna-NM_001005484.2   | mRNA  | MANE Select  | NC_000001.11 | NM_001005484.2 | NM_001005484 | 2            | exon | OR4F5 | olfactory receptor family 4 subfamily F member 5             | 1          | 1     | 65419  | 65433  | .     | +      | .     | ENST00000641515.2 | 79501  | NM_001005484.2 | 14825 |       |              | 14     | 65419              | 65433            |
| BestRefSeq | exon-NM_001005221.2-1 | rna-NM_001005221.2   | mRNA  | MANE Select  | NC_000001.11 | NM_001005221.2 | NM_001005221 | 2            | exon | OR4F29| olfactory receptor family 4 subfamily F member 29           | 1          | 1     | 450740 | 451678 | .     | -      | .     | ENST00000426406.4 | 729759 | NM_001005221.2 | 31275 |       |              | 938    | 450740             | 451678           |
| BestRefSeq | exon-NM_001005277.1-1 | rna-NM_001005277.1   | mRNA  | MANE Select  | NC_000001.11 | NM_001005277.1 | NM_001005277 | 1            | exon | OR4F16| olfactory receptor family 4 subfamily F member 16           | 1          | 1     | 685716 | 686654 | .     | -      | .     | ENST00000332831.5 | 81399  | NM_001005277.1 | 15079 |       |              | 938    | 685716             | 686654           |
| BestRefSeq | exon-NM_015658.4-1   | rna-NM_015658.4      | mRNA  | MANE Select  | NC_000001.11 | NM_015658.4    | NM_015658     | 4            | exon | NOC2L | NOC2 like nucleolar associated transcriptional repressor | 1          | 1     | 959215 | 959256 | .     | -      | .     | ENST00000327044.7 | 26155  | NM_015658.4    | 24517 | 610770 |              | 41     | 959215             | 959256           |

- **cnvkitReference** <br>
Path to a CNVkit Copy number reference profile `.cnn` file built with your normal samples <br>
- **cbr_bed** <br>
Path to a contiguous baited region `.bed` file which contains your baited regions of the genome. *note that the 4th column contains a list of gene_exons which are overlapped by the range* <br>
>The header row is for explanatory purposes
>
| Chromosome | Start    | End      | Gene_Exon List       |
|------------|----------|----------|----------------------|
| chr1       | 10000    | 20000    | ['OR4F5_1', 'OR4F5_2'] |
| chr2       | 10000    | 20000    | ['TNP1_1', 'TNP1_2']   |
| chr10      | 10000 | 20000 | ['COX15_9']            |

- **contamination_bed** <br>
Path to a `.bed` file which has genomic regions that should be excluded from analysis <br>

- **smartbinmeta** <br>
Path to a `.csv` file which contains the size of bins used for the Fragment Pillar (example below):<br>

| BinIndex | Chrom | StartPos | EndPos   |
|----------|-------|----------|----------|
| 1        | chr1  | 0        | 1313839  |
| 2        | chr1  | 1313840  | 2558334  |
| 3        | chr1  | 2558335  | 6339481  |
| 4        | chr1  | 6339482  | 10408056 |
| 5        | chr1  | 10408057 | 11232432 |

- **fragExclusionRegions** <br>
Path to a `.csv` file which contains the contiguous baited regions (in the same format as cbr_bed `.bed`) the  which we want to exclude from the Fragment Pillar aggregation<br>
- **Common_CH**<br>
Path to a `.csv` file which contains common CHIP mutations (example below):<br>
| CPRV               |
|--------------------|
| chr10|100425|TA|T  |
| chr10|10027425|T|A |
| chr10|100499032|GT|G|
| chr10|100988900|A|G |
| chr10|101008657|GCCT|G |
| chr10|101009374|A|G |

- **binMeta**<br>
Path to a `.pkl` file which contains a dictionary of CHROMOSOME|GENE keys and bin values used for copy number genomic binning (example below):<br>
```
{'chr1|OR4F5': 'chr1_1',
 'chr1|OR4F29': 'chr1_1',
 'chr1|OR4F16': 'chr1_1',
 'chr1|SAMD11': 'chr1_1',
 'chr1|NOC2L': 'chr1_1',
 'chr1|CDK11B': 'chr1_2',
 'chr1|SLC35E2B': 'chr1_2',
 'chr1|CDK11A': 'chr1_2',
 'chr1|NADK': 'chr1_2',
 'chr1|GNB1': 'chr1_2',
 'chr1|CALML6': 'chr1_2',
 'chr1|AJAP1': 'chr1_3',
 'chr1|NPHP4': 'chr1_3',
 'chr1|KCNAB2': 'chr1_3',...}
```
