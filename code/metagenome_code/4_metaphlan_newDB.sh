#!/bin/bash
#SBATCH --job-name=metaphlan
#SBATCH --nodes=1
#SBATCH --nodelist=heracpu02.nlzoh.si
#SBATCH --cpus-per-task=50
#SBATCH --mem=100G
#SBATCH --time=5-24:00:00
#SBATCH --output=result_%j.log        # Standard output and error log

eval "$(conda shell.bash hook)"
source activate metaphlan

path="/volumes/homehpc/storage/finished_projects/UrsaM/SQM_Microbiota/metaphlan"

for f in "${path}"/*bowtie2.bz2; do
    sample_name=$(basename "${f}" .bowtie2.bz2)
    profile_out="${path}/profiled2/profiled_${sample_name}.txt"
    metaphlan "$f" --nproc 50 --input_type mapout --db_dir /home/storage/DB/metaphlan4.2.2 -o "${profile_out}"
done

cd $path 
mkdir -p ./gtdb_profiled

for f in profiled2/*; do
  sample_name=$(basename "$f")
  sgb_to_gtdb_profile.py -i "$f" -o "gtdb_profiled/gtdb_profiled_${sample_name}.txt"
done

merge_metaphlan_tables.py ./gtdb_profiled/gtdb*.txt > ./gtdb_profiled/gtdb_merged_abundance_table.txt

