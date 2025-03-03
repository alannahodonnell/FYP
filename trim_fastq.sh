#!/bin/bash

# Path to Trimmomatic JAR file
TRIMMOMATIC_JAR="trimmomatic"
# Path to adapter file
ADAPTER_FILE="$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-SE.fa"

# Loop through all FASTQ files in the current directory
for file in *.fastq; do
    # Extract the base name (SRR6175537, SRR6175538, SRR6175539, SRR6175540, SRR6175541, SRR6175542, SRR6175543, SRR6175544.fastq)
    BASENAME=$(basename "$file" .fastq)
    # Define the output file name
    OUTPUT_FILE="${BASENAME}.trimmed.fastq"

    # Run Trimmomatic
    echo "Trimming $file..."
    "$TRIMMOMATIC_JAR" SE -phred33 \
        "$file" "$OUTPUT_FILE" \
        ILLUMINACLIP:"$ADAPTER_FILE":2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    echo "Finished trimming $file. Output: $OUTPUT_FILE"
done

echo "Trimming completed for all FASTQ files in the current directory."
