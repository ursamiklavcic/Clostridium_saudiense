#!/bin/bash

# Concatenate all files in each directory and name the output file after the directory name
for f in H00*; do cat "$f"/* > "${f}.cat.fastq.gz"; done


####################################### PLOT #######################################

# Run NanoPlot on each concatenated file
for f in *.fastq.gz; do     
    NanoPlot -t 8 --fastq "$f" --no_static --no_supplementary --tsv_stats --outdir "${f%.cat.fastq.gz}.plot"; 
done

# Rename NanoStats.txt files to reflect their folder names
for p in */NanoStats.txt; do 
    [ -f "$p" ] && mv -i "$p" "${p%/*}/${p%/*}.NanoStats.txt"; 
done

# Create a 'stats' directory and copy all renamed NanoStats files there
mkdir -p stats
cp */*NanoStats* stats/

# Add filenames as headers and merge all NanoStats into merged.tsv without redundancy
awk -i inplace -v ORS='\r\n' 'FNR==1{print FILENAME}1' stats/*.NanoStats.txt
paste stats/*.NanoStats.txt > stats/merged.tsv 

#################################### FILTER ########################################

# Filter sequences by quality and length, save output in .filtered.fastq.gz files
for f in *cat.fastq.gz; do
    gunzip -c "$f" | NanoFilt -q 10 -l 3000 | gzip > "${f%.cat*}.filtered.fastq.gz";
done

################################# PLOT AGAIN ########################################

# Run NanoPlot again on filtered files
for f in *.filtered.fastq.gz; do     
    NanoPlot -t 8 --fastq "$f" --no_static --no_supplementary --tsv_stats --outdir "${f%.fastq.gz}.plot"; 
done

# Rename new NanoStats files to include folder names and move them to a new stats directory
for p in *filtered*/NanoStats.txt; do 
    [ -f "$p" ] && mv -i "$p" "${p%/*}/${p%/*}.NanoStats.txt"; 
done

# Create 'stats2' directory and copy all new NanoStats files there
mkdir -p stats2
cp *filtered*/*NanoStats* stats2/

# Add filenames as headers and merge all filtered NanoStats into merged2.tsv
awk -i inplace -v ORS='\r\n' 'FNR==1{print FILENAME}1' stats2/*.NanoStats.txt
paste stats2/*.NanoStats.txt > stats2/merged2.tsv 

###################### ORGANIZE OUTPUT ######################

# Create directories and move files for organization
mkdir -p 1plots 2sequences 3filtered_sequences
mv *.plot 1plots/
mv *cat.fastq* 2sequences/
mv *filtered* 3filtered_sequences/
