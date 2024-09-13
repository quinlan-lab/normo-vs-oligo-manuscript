#!/bin/bash

import os
import glob

projects = ["19610R", "19959R", "20070R", "20963R"]

def make_sample_sheet(project):
    
    sample_sheet = f"/scratch/ucgd/lustre-labs/quinlan/u1240855/fastq_to_duplex_bam/sample_sheets/{project}_sample_sheet.csv"
    merged_fastq_dir = "/scratch/ucgd/lustre-labs/quinlan/u1240855/fastq_to_duplex_bam/output/merged_fastqs"
    structure = "8M1S+T"

    with open(sample_sheet, 'w') as f:
        f.write("sample,fastq_1,fastq_2,read_structure\n")

    for filename in glob.glob(os.path.join(merged_fastq_dir, project, "*_R1.merged.fastq.gz")):
        sample = os.path.basename(filename).replace("_R1.merged.fastq.gz", "")

        print(sample)

        R2 = filename.replace("_R1", "_R2")
        with open(sample_sheet, 'a') as f:
            f.write(f"{sample},{filename},{R2},\"{structure} {structure}\"\n")

for project in projects:
    make_sample_sheet(project)