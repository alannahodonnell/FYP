#!/bin/bash

# List of SRA accessions
sra_files=(
    SRR6175537
    SRR6175538
    SRR6175539
    SRR6175540
    SRR6175541
    SRR6175542
    SRR6175543
    SRR6175544
)

# Loop through each SRA accession and download the FASTQ files
for sra in "${sra_files[@]}"
do
    echo "Downloading ${sra}..."
    fastq-dump --split-3 --gzip $sra
done
