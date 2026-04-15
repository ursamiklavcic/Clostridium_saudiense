# This is everything I did in this folder: 


# Prepare bam files for anvio 
### qc_wgs.sh script used
#!/bin/bash
#SBATCH --job-name=qc_bam_only
#SBATCH --cpus-per-task=20
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=100M
#SBATCH --time=20:00:00
#SBATCH --output=qcjob_%j.out
#SBATCH --error=qcjob_%j.err

qc_path="/home/storage/finished_projects/UrsaM/C_saudiense"

source /home/nlzoh.si/larbez1/miniconda3/bin/activate /home/nlzoh.si/larbez1/miniconda3/envs/wgs_qc

# BOWTIE2
bowtie2-build --threads 20 01_FASTA/reference.fasta 03_BAM/reference_index

cd 03_BAM/

for file in "${qc_path}"/2trimmomatic/*R1_P1.trim.fastq.gz; do
  fq_name=$(basename "$file" _R1_P1.trim.fastq.gz)
  sample_id=$(echo "$fq_name" | cut -d'_' -f1,2)

  bowtie2 --fr --end-to-end --threads 16 --sensitive \
    -x reference_index \
    -1 "$file" \
    -2 "${file%R1_P1.trim.fastq.gz}R2_P1.trim.fastq.gz" \
    -U "${file%R1_P1.trim.fastq.gz}R1_R2_P3.trim.fastq.gz" \
    -S ./${fq_name}.sam 
done

# SAMTOOLS
# Convert SAM to BAM
for file in ./*.sam; do 
  name=$(basename "$file" .sam)
  id=$(echo "$name" | cut -d'_' -f1,2)
  samtools view -S -b "${name}.bam > ./${id}.bam
  samtools sort ${file} -o ${id}.bam
  samtools index ${id}.bam
done 

## metaG samples H005, H007, H009 
#cp MH001 MH002 MH003 MH004 MH005 MH006 MH007 MH008 MH009 MH010 MH011 MH012 to 01_BAM

for id in MH001 MH002 MH003 MH004 MH005 MH006 MH007 MH008 MH009 MH010 MH011 MH012; do
    echo "Working on ${id} ..."
    bowtie2 --fr --end-to-end --threads 20 --sensitive -x reference_index -1 ${id}.1.fq.gz -2 ${id}.2.fq.gz -S ${id}.sam
    echo "Moving to samtools on ${id}..."
    samtools view -F 4 -bS ${id}.sam -o ${id}-RAW.bam
    samtools sort ${id}-RAW.bam -o ${id}.bam
    samtools index ${id}.bam
done

rm *-RAW.bam

conda deactivate 


# Coverage 
source /home/nlzoh.si/larbez1/miniconda3/bin/activate /home/nlzoh.si/larbez1/miniconda3/envs/cge_env

for file in 01_BAM/*.bam; do
  ime=$(basename "$file" .bam)
  coverage_file=${ime}_"cov_stat.txt"
  genomeCoverageBed -ibam 01_BAM/${ime}.bam -d > "$coverage_file"
done

mkdir -p 01_COVERAGE
mv *_cov_stat.txt 01_COVERAGE/
cd 01_COVERAGE 

reference_length=3707884

output_file="./alignment_summary_exclude.txt"
echo -e "Sample\tAligned_Length\t%\tStd" > "$output_file"
for file in ./*cov_stat.txt; do
    sample=$(basename "$file" _cov_stat.txt)
    aligned_length=$(awk '$3 > 0 {count++} END {print count}' "$file")
    std=$(awk '{sum+=$3; sumsq+=$3*$3} END {print sqrt(sumsq/NR - (sum/NR)**2)}' "$file")
    percentage_aligned=$(awk -v aligned="$aligned_length" -v ref="$reference_length" 'BEGIN {printf "%.4f", (aligned / ref) * 100}')
    echo -e "$sample\t$aligned_length\t$percentage_aligned%\t$std" >> "$output_file"
done

paste ./*cov_stat.txt | awk '{
    pos=$2; all_positive=1;
    for (i=3; i<=NF; i+=3) if ($i <= 0) all_positive=0;
    if (all_positive) print pos; }' > ./shared_positions.txt

shared_count=$(wc -l < ./shared_positions.txt)
percentage=$(awk "BEGIN {print ($shared_count / $reference_length * 100)}")

echo -e "\nSummary:" >> "$output_file"
echo -e "Total Reference Length:\t$reference_length" >> "$output_file"
echo -e "Positions Aligned in All Samples:\t$shared_count" >> "$output_file"
echo -e "Percentage of All-Aligned Positions:\t$percentage%" >> "$output_file"


cd ../ 

## ANVI'O

# using the anvio-8
# Firstly following this tutorial: https://merenlab.org/tutorials/vibrio-jasicida-pangenome/
# and https://anvio.org/tutorials/pangenome-graphs/
# make sure to run  - & anvi-setup-ncbi-cogs if its your first time using anvio-8 or anvio of any kind 

# ISOLATE GENOMES 

# Reformat fasta files, so that anvi'o likes them 
path='~/projects/C_saudiense/01_FASTA'

while IFS=$'\t' read -r sample path; do
    output="/home/nlzoh.si/ursmik1/projects/C_saudiense/01_FASTA/${sample}.fasta"
    anvi-script-reformat-fasta "$path" -l 1000 -o "$output"
done < fastas.txt
# removed 30-50% of all contigs, but up to 1% off all nucleotides for all samples 

# Create fasta.txt 
echo -e 'name\tpath' > fasta.txt

for filename in 01_FASTA/*.fasta; do
    echo -e $(basename $filename | cut -d. -f1)'\t'$(realpath $filename) >> fasta.txt
done

# Create contigs-db files = anvio's fasta files! 
mkdir 02_CONTIGS/

for f in 01_FASTA/*.fasta; do
    echo
    echo "Working on $f ..."
    echo

    # Extract just the base filename, e.g., "H005_124"
    base=$(basename "$f" .fasta)

    anvi-gen-contigs-database -f "$f" \
                              -o "/home/nlzoh.si/ursmik1/projects/C_saudiense/contigs_db/C_saudiense_${base}.db" \
                              --num-threads 50 \
                              -n "C_saudiense_${base}"
done

# Annotate contigs databases 
for g in 02_CONTIGS/*.db
do
    anvi-run-hmms -c $g --num-threads 50
    anvi-run-ncbi-cogs -c $g --num-threads 50
    anvi-scan-trnas -c $g --num-threads 50
    anvi-run-scg-taxonomy -c $g --num-threads 50
done

anvi-display-contigs-stats 02_CONTIGS/*db
anvi-script-gen-genomes-file --input-dir ./02_CONTIGS/ -o external-genomes.txt

# Investigate contamination 
anvi-estimate-genome-completeness -e external-genomes.txt
# Genome H005_126 contains 2 ganeomes, based on reduncdancy score and total lenght! REMOVED FROM ANALYSIS! and from external-genomes.txt

mkdir 03_PAN/

# make one CONTIGS.db named *-GENOMES.db
anvi-gen-genomes-storage -e external-genomes.txt -o 03_PAN/C_saudiense-GENOMES.db
# compute pangenome
anvi-pan-genome -g 03_PAN/C_saudiense-GENOMES.db \
                --project-name C_saudiense \
                --num-threads 24
                
# ANI SCORE 
anvi-compute-genome-similarity --external-genomes external-genomes.txt \
                               --program pyANI \
                               --output-dir 04_ANI \
                               --num-threads 24 \
                               --pan-db 03_PAN/C_saudiense/C_saudiense-PAN.db
                               
# DISPLAY and save: 
anvi-display-pan -p 03_PAN/C_saudiense/C_saudiense-PAN.db -g 03_PAN/C_saudiense-GENOMES.db

 
# Defined CORE GENOME for this species (by hand in interactive mode) 
# Get aligned sequences out
anvi-get-sequences-for-gene-clusters -g 03_PAN/C_saudiense-GENOMES.db -p 03_PAN/C_saudiense/C_saudiense-PAN.db -b CORE_GENOME --max-num-genes-from-each-genome 1 --min-geometric-homogeneity 1 --max-functional-homogeneity 0.95 --concatenate-gene-clusters -o 03_PAN/concatenated_genes.fasta

# Create phylogenetic tree
anvi-gen-phylogenomic-tree -f 03_PAN/concatenated_genes.fasta -o 03_PAN/phylogenetic_tree_CORE_GENOME

# provide the path to an empty PROFILE.db and create the manual interactive mode 
anvi-interactive -p 03_PAN/PROFILE-tree.db -t 03_PAN/phylogenetic_tree_CORE_GENOME --manual

 # https://merenlab.org/2016/11/08/pangenomics-v2//#calculating-rarefaction-curves-and-heaps-law-parameters
 # This is not availabile yet in v8, but will be in v9, so for the PhD thesis I could do it but not right now! 

# Functional analysis of the genomes 
# layers_additional_data.txt availabile in Excel in the folder metadata_HPC.xlsx 

anvi-import-misc-data layers_additional_data.txt -p 03_PAN/C_saudiense/C_saudiense-PAN.db --target-data-table layers

anvi-display-pan -g 03_PAN/C_saudiense-GENOMES.db -p 03_PAN/C_saudiense/C_saudiense-PAN.db

# FUNCTIONAL ENRICHMENT AB and time 
anvi-compute-functional-enrichment-in-pan -p 03_PAN/C_saudiense/C_saudisense-PAN.db
                                          -g  03_PAN/C_saudiense-GENOMES.db \
                                          --category ANTIBIOTIC_THERAPY \
                                          --annotation-source COG20_FUNCTION \
                                          -o 03_PAN/enriched-functionsAB.txt
# Filtered file for this (based on q_values, in R and saved as enriched_functionsAB_filtered.csv 
# the same was done for TIME (that is why it is a categorical value (this program does not work with numeric values! 
anvi-compute-functional-enrichment-in-pan -p 03_PAN/C_saudiense/C_saudiense-PAN.db  
                                          -g 03_PAN/C_saudiense-GENOMES.db --category TIME 
                                          --annotation-source COG20_FUNCTION 
                                          -o 03_PAN/enriched_functions_TIME.txt
                                          
# At a leter data 15.9.2025; also performed functional enrichment analysis for strain groups, based on ANI score. 


# SNIPPY 
# In the folder snippy/snippy.sh script to obtain SNPs

# Interactive session 
srun --time=1-24:00 --mem=50G --cpus-per-task=20 --job-name=snippy --pty /bin/bash
conda activate snippy 

cd snippy/

snippy-multi reads.tab --ref 01_FASTA/reference.fasta --cpus 20 > runme.sh
less runme.sh

sh ./runme.sh

# FastTree
conda activate fasttree

FastTree -nt -gtr -gamma snippy_2/core.full.aln > snippy/core.tree

# plot_snippy.R 
# results in out/snippy 

# SNIPPY with Gubbins (removal of polymorphic regions!) 

# delete reference out of the aligment file 
sed '/^>reference$/,/^>/ {/^>reference$/d; /^>/!d}' ./core.full.aln > ./core_no_ref.full.aln
snippy-clean_full_aln core_no_ref.full.aln > clean.full.aln

rm core_no_ref.full.aln

# https://github.com/nickjcroucher/gubbins

mkdir -p gubbins/
cd gubbins

source /home/nlzoh.si/larbez1/miniconda3/bin/activate /home/nlzoh.si/larbez1/miniconda3/envs/snp_env

run_gubbins.py -c 20 --prefix gubb ../clean.full.aln

conda activate snippy 

snp-sites -c -o ./Csaudiense_clean.core.fasta gubbins/gubb.filtered_polymorphic_sites.fasta
snp-dists ./Csaudiense_clean.core.fasta > snippy_gubb_dists.tab

# plot_snippy.R 
# results in out/snippy line 226 onwards 

#plot_snps_pertime_group.R


# Prokka and Roary for pangenome 
# srun 

source /home/nejc/miniconda3/bin/activate /home/nejc/miniconda3/envs/prokka

cd ~/projects/C_saudiense/01_FASTA/

for k in *.fasta; do prokka $k --outdir ../02_prokka/"$k".prokka.output --prefix PROKKA_$k; echo $k; done 

conda deactivate

source /home/nlzoh.si/leomar1/miniconda3/bin/activate /home/nlzoh.si/leomar1/miniconda3/envs/roary

cd 02_prokka
mkdir GFF

cp ./*/*.gff GFF

roary -n -e -p 8 ./GFF/*.gff -f roary

conda deactivate

# Load into phandango with metadata etc. 




# ONLY METAGENOMES! # prefix to folders 09_


## StrainPhlan 
/home/stoage/finished_projects/UrsaM/SQM_Microbiota/4_1_strainphlan.sh

# Download pyphlan github repo git clone into ~/bin
# Extract pairwise phylogenetic distances (normalized by th pyhlogenetic distance) into .tsv file 
~/bin/pyphlan/tree_pairwisedists.py -n output/t__SGB6178_mutation_rates/RAxML_bestTree.t__SGB6178.StrainPhlAn4.tre ./output/t__SGB6178_mutation_rates/t__SGB6178_nGD.tsv

# Setting species-specific strain identity thresholds (this was done before, we do not need to repeat it) 0.03

# but R code, to be sure what our threshold for CLostridium saudiense is! 
~/projects/C_saudiense/00_CODE/determine_threshold.R 

# Detection of strain sharing events 
dir="/volumes/homehpc/storage/finished_projects/UrsaM/SQM_Microbiota/strainphlan"

strain_transmission.py --tree "${dir}/output/t__SGB6178_mutation_rates/RAxML_bestTree.t__SGB6178.StrainPhlAn4.tre" --metadata /home/storage/finished_projects/UrsaM/SQM_Microbiota/strainphlan/metadata.tsv --output_dir strainphlan/output/strain_transmission/ --threshold 0.03


## StrainEr 
~/projects/C_saudiense/strainer/README.txt 

# Run, but I think my database for reference is too small so there were no SNPs that would be fixed in my strains genomes
# I would have to employ a larger database 

## LongStrain
~/projects/C_saudiense/LongStrain/README.txt

# Everthing ran as it should after a few corrections by Boyan Zhao. 
# It outputs only the dominant and alternative strain not the number of strains! 

# ConStrains
# DOES NOT WORK, as it is old (2015), python=2.7 but the biggest problem was metaphlan2 dependancy as I can not find metaphlan2 DB.. 

# MIDAS 
# Clostridium saudiense does not have enough coverage in most samples (only in SH) to get allele trajectories that make any sense! 
# But will continue with this as I have the results at least! 


## 
# Additionaly sequeneced genomes H005_73, H005_98, H005_118, H007_176, H009_221 with Nanopore (long-reads) 
# Filtering and hybrid assembly in folder /home/storage/UrsaM/C_saudiense/long_read_v10/ 
# Copied .fasta files to 01_FASTA_hybrid

# Run DRAM (insted of prokka) 00_CODE/DRAM.sb

# Run pyANI
# Run prokka, roary, phandango 
# Run snippy, gubbins, mega12
 


