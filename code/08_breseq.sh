#!/bin/bash

# Set reference GenBank file
REFERENCE="/home/nlzoh.si/ursmik1/projects/C_saudiense/02_prokka/H005_118.fasta.prokka.output/PROKKA_H005_118.fasta.gbk"

# Output directory for breseq results
OUTPUT_DIR="/home/nlzoh.si/ursmik1/projects/C_saudiense/03_breseq_results"
qc_path="/home/storage/finished_projects/UrsaM/C_saudiense"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all *_1.fastq.gz files to find pairs
for READ1 in "${qc_path}/2trimmomatic/*R1_P1.trim.fastq.gz"; do
    # Ensure READ1 exists
    if [ ! -f "$READ1" ]; then
        echo "File $READ1 not found. Skipping..."
        continue
    fi

    # Derive the corresponding READ2 file
    READ2="${READ1/R1_P1.trim.fastq.gz/R2_P1.trim.fastq.gz}"

    # Ensure READ2 exists
    if [ ! -f "$READ2" ]; then
        echo "File $READ2 not found. Skipping..."
        continue
    fi

    # Get the sample name by stripping the file extension (adjust as needed)
    SAMPLE_NAME=$(basename "$READ1" _R1_P1.trim.fastq.gz)

    # Debugging: Show which files are being processed
    echo "Processing $SAMPLE_NAME with $READ1 and $READ2"

    # Run breseq for each pair
    breseq -r "$REFERENCE" "$READ1" "$READ2" -o "$OUTPUT_DIR/$SAMPLE_NAME"
done

echo "All samples processed. Results are in the $OUTPUT_DIR directory."

# Fetch .gd files into a separate folder 
cd $OUTPUT_DIR
mkdir gd_files

destination=~/projects/C_saudiense/03_breseq_results/gd_files

for dir in */ ; do
  if [ -d "$dir/output" ]; then
    folder_name=$(basename "$dir")
    if [ -f "$dir/output/output.gd" ]; then
      cp "$dir/output/output.gd" "$destination/${folder_name}.gd"
    fi
  fi
done

# Compare mutattions in genomes 
gdtools COMPARE -o ../compare.html -r "$REFERENCE" `ls *.gd`

gdtools COUNT -o ../count.csv -r "$REFERENCE" `ls *.gd`

# Only mutations present in all 
gdtools UNION -o ../unique.gd  `ls *.gd`
gdtools COUNT -o ../unique_count.csv -r "$REFERENCE" ../unique.gd

# For here R script breseq.R 
