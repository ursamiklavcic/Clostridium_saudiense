#!/bin/bash

# Year and runID provided by the user
year="$1"             # provide year of the run
runID="$2"            # provide sequencing run ID
expected_count="$3"   # provide expected number of fastq reads

# Generate the full path
raw_path="/volumes/homehpc/storage/raw_data/NextSeq/${year}/${runID}"
qc_path="/volumes/homehpc/storage/qc_reads/NextSeq/${year}/${runID}"

mkdir "$qc_path"

# Check if the folder exists
if [ ! -d "$raw_path" ]; then
    echo "Folder not found: $raw_path"
    exit 1
fi

# Check if the number of files stored matches the number you expect
actual_count=$(ls -1 ${raw_path}/*/*fastq.gz | wc -l)

if [ "$actual_count" -ne "$expected_count" ]; then
  echo "Error: The number of files in folder ($actual_count) does not match the expected number ($expected_count)."
  exit 1
fi


mkdir "${raw_path}/renamed"
cp "${raw_path}"/*/*fastq.gz -t "/${raw_path}/renamed" && rename 's/WGS-//' "${raw_path}/renamed"/*

chmod -R 777 $raw_path
chmod -R 777 $qc_path

#SLURM job script
job_script=$(cat <<EOF
#!/bin/bash
#SBATCH --job-name=qc
#SBATCH --cpus-per-task=24
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=100M
#SBATCH --time=42:00:00
#SBATCH --nodes=1
#SBATCH --output=qcjob_%j.out   # Use unique output filename based on job ID
#SBATCH --error=qcjob_%j.err    # Use unique error filename based on job ID
#SBATCH --nodelist=heracpu02.nlzoh.si

eval "\$(conda shell.bash hook)"
source activate wgs_qc

raw_path="$raw_path"
qc_path="$qc_path"

cd \$qc_path

#rm -r \$qc_path/*

mkdir \$qc_path/1fastqc_pretrim
mkdir \$qc_path/2trimmomatic
mkdir \$qc_path/3fastqc_posttrim
mkdir \$qc_path/4spades_results
mkdir \$qc_path/5quast_results
mkdir \$qc_path/6kraken2
mkdir \$qc_path/7bowtie2
mkdir \$qc_path/8samtools


#FASTQC BEFORE TRIM
echo ''
echo 'RUNNING FASTQC BEFORE TRIM'
echo ''

for file in \$raw_path/renamed/*gz; do
    while [ "\$(jobs -p | wc -l)" -ge "\$SLURM_NTASKS" ]; do
        sleep 30
    done
    srun --ntasks=1 --cpus-per-task=12 --exclusive fastqc -o "$qc_path"/1fastqc_pretrim "\$file" & 
done
wait

# Wait for all FastQC jobs to finish
wait

echo ''
echo 'DONE RUNNING FASTQC'
echo ''


echo ''
echo 'RUNNING MULTIQC'
echo ''
multiqc \$qc_path/1fastqc_pretrim/ -k tsv -o ./1fastqc_pretrim
rm \$qc_path/1fastqc_pretrim/*.zip

#TRIMMOMATIC
echo ''
echo 'RUNNING TRIMMOMATIC'
echo ''
for file in "$raw_path"/renamed/*R1_001.fastq.gz; do
    while [ "\$(jobs -p | wc -l)" -ge "\$SLURM_NTASKS" ]; do
        sleep 30
    done
    filename=\$(basename "\$file")   # Corrected variable expansion
    echo "\$filename"                 # Print the filename for debugging
    srun --ntasks=1 --cpus-per-task=16 trimmomatic PE -threads 48 -phred33 "\$file" "\${file%R1_001.fastq.gz}R2_001.fastq.gz" "\$qc_path"/"\${filename%R1_001.fastq.gz}R1_P1.trim.fastq.gz" "\$qc_path"/"\${filename%R1_001.fastq.gz}R1_P2.trim.fastq.gz" "\$qc_path"/"\${filename%R1_001.fastq.gz}R2_P1.trim.fastq.gz" "\$qc_path"/"\${filename%R1_001.fastq.gz}R2_P2.trim.fastq.gz" ILLUMINACLIP:/home/storage/DB/NEBPE-PE.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:40 &
done
wait

echo ''
echo 'DONE WITH RUNNING TRIMMOMATIC'
echo ''

rm ${qc_path}/*P2.trim.fastq.gz
mv ${qc_path}/*.trim.fastq.gz ${qc_path}/2trimmomatic/.

#FASTQC AFTER TRIM
for file in ${qc_path}/2trimmomatic/*trim.fastq.gz; do
        while [ "\$(jobs -p | wc -l)" -ge "\$SLURM_NTASKS" ]; do
                sleep 30
        done
        srun --ntasks=1 --cpus-per-task=24 fastqc --exclusive -o "${qc_path}"/3fastqc_posttrim "\$file" &
done
wait


multiqc ${qc_path}/3fastqc_posttrim/ -k tsv -o ${qc_path}/3fastqc_posttrim
rm \$qc_path/3fastqc_posttrim/*.zip


#SPADES
echo ''
echo 'RUNNING SPADES'
echo ''

for file in ${qc_path}/2trimmomatic/*R1_P1.trim.fastq.gz;  do
	while [ "\$(jobs -p | wc -l)" -ge "\$SLURM_NTASKS" ]; do
        	sleep 30
	done
	f=\$(basename "\$file")
	srun --ntasks=1 --cpus-per-task=16 spades.py -1 "\${file}" -2 "\${file%R1_P1.trim.fastq.gz}R2_P1.trim.fastq.gz" --careful --threads 38 --cov-cutoff auto -o "\${qc_path}"/4spades_results/"\${f%R1_P1.trim.fastq.gz}output" &
done
wait

echo ''
echo 'DONE WITH RUNNING SPADES'
echo ''

#RENAIMNG CONTIG.FASTA
for f in ${qc_path}/4spades_results/*/contigs.fasta; do
        newname=\$(basename "\${f%/contigs.fasta}.fasta")
        cp \$f \${qc_path}/4spades_results/"\$newname"
done

rm -r \${qc_path}/4spades_results/*_output


#QUAST
echo ''
echo 'RUNNING QUAST'
echo ''
quast ${qc_path}/4spades_results/*.fasta
mv ${qc_path}/quast_results ${qc_path}/5quast_results


echo ''
echo 'DONE WITH RUNNING QUAST'
echo ''

#KRAKEN
echo ''
echo 'RUNNING KRAKEN2'
echo ''
for file in ${qc_path}/2trimmomatic/*R1_P1.trim.fastq.gz; do
        while [ "\$(jobs -p | wc -l)" -ge "\$SLURM_NTASKS" ]; do
                sleep 30
        done
	f=\$(basename "\$file")
	srun --ntasks=1 --cpus-per-task=16 kraken2 --threads 16 --db /home/storage/DB/minikraken2_v3_16GB_20240409/ --report "\${f%R1_P1.trim.fastq.gz}report" --gzip-compressed --paired "\$file" "\${file%R1_P1.trim.fastq.gz}R2_P1.trim.fastq.gz" > "\${f%_R1_P1.trim.fastq.gz}.kraken" &
done
wait

echo ''
echo 'DONE WITH RUNNING KRAKEN'
echo ''

#KRONA
echo ''
echo 'RUNNING KRONA'
echo ''
for file in ${qc_path}/*.kraken; do
        while [ "\$(jobs -p | wc -l)" -ge "\$SLURM_NTASKS" ]; do
                sleep 30
        done
        cat "\$file" | cut -f 2,3 > "\${file%.kraken}.krona";
        srun --ntasks=1 --cpus-per-task=16 ktImportTaxonomy "\${file%.kraken}.krona" -o "\${file##*/}.html" -tax /home/storage/DB/krona/taxonomy &
done
wait

rm ${qc_path}/*kraken;
rm ${qc_path}/*krona;
rm -r ${qc_path}/*html.files; 

mv ${qc_path}/*html ${qc_path}/6kraken2;
mv ${qc_path}/*_report ${qc_path}/6kraken2;

echo ''
echo 'DONE WITH RUNNING KRONA'
echo ''

#BOWTIE2
echo ''
echo 'RUNNING BOWTIE2'
echo ''

for file in ${qc_path}/4spades_results/*_output.fasta; do
        while [ "\$(jobs -p | wc -l)" -ge "\$SLURM_NTASKS" ]; do
                sleep 30
        done
        srun  --ntasks=1 --cpus-per-task=16 bowtie2-build -f \$file "\${file%_output.fasta}"
        mv ${qc_path}/4spades_results/*.bt2 ${qc_path}/7bowtie2
done
wait
for file in ${qc_path}/2trimmomatic/*R1_P1.trim.fastq.gz; do
	 while [ "\$(jobs -p | wc -l)" -ge "\$SLURM_NTASKS" ]; do
                sleep 30
        done;
	sample=\$(basename "\$file")
 
  echo "WORKING ON \$sample"
  srun --ntasks=1 --cpus-per-task=24 bowtie2 --fr --end-to-end --threads 24 -x ${qc_path}/7bowtie2/"\${sample%_R1_P1.trim.fastq.gz}" 2>"\${sample%R1_P1.trim.fastq.gz}stats.txt" -sensitive -1 "\${file}" -2 "\${file%R1_P1.trim.fastq.gz}R2_P1.trim.fastq.gz" -S "\${sample%_R1_P1.trim.fastq.gz}.sam" &
done
wait

mv ${qc_path}/*_stats.txt ${qc_path}/7bowtie2
mv ${qc_path}/*.sam ${qc_path}/7bowtie2

echo ''
echo 'DONE WITH RUNNING BOWTIE2'
echo ''

SAMTOOLS
echo ''
echo 'RUNNING SAMTOOLS'
echo ''

for file in ${qc_path}/7bowtie2/*.sam; do
        while [ "\$(jobs -p | wc -l)" -ge "\$SLURM_NTASKS" ]; do
                sleep 30
        done
        srun --ntasks=1 --cpus-per-task=16 samtools view -S -b \${file} > "\${file##*/}".bam &
             
done
wait

for file in ${qc_path}/7bowtie2/*.sam; do
        while [ "\$(jobs -p | wc -l)" -ge "\$SLURM_NTASKS" ]; do
                sleep 30
        done
        srun --ntasks=1 --cpus-per-task=16 samtools sort "\${file##*/}".bam -o "\${file##*/}".sorted.bam &
              
done
wait

for file in ${qc_path}/7bowtie2/*.sam; do
        while [ "\$(jobs -p | wc -l)" -ge "\$SLURM_NTASKS" ]; do
                sleep 30
        done
        srun --ntasks=1 --cpus-per-task=16 samtools depth "\${file##*/}".sorted.bam | awk '{sum+=\$3} END { print "Average =",sum/NR}' > "\${file##*/}.txt" & 
        
done
wait

echo ''
echo 'DONE WITH RUNNING SAMTOOLS'
echo ''

rm ${qc_path}/7bowtie2/*.sam;
rm ${qc_path}/7bowtie2/*.bt2;
rm ${qc_path}/*.bam;
mv ${qc_path}/*sam.txt ${qc_path}/8samtools;

chmod -R 777 ${qc_path}


#Joining output files from SAMTOOLS and BOWTIE2
output_file="\${qc_path}/8samtools/joined_samtools.txt"

# Header for the output file
echo -e "Filename\tCoverage" > "\$output_file"

# Iterate over each file
for file in ${qc_path}/8samtools/*.txt; do
    # Extract filename without extension
    filename=\$(basename "\$file" .sam.txt)
    # Extract coverage information
    coverage=\$(awk '/Average/ {print \$3}' "\$file")
    # Append to the output file
    echo -e "\$filename\t\$coverage" >> "\$output_file"
done

output_file="\${qc_path}/7bowtie2/joined_bowtie2.txt"

# Header for the output file
echo -e "Filename\tAlignment" > "\$output_file"

# Iterate over each file
for file in ${qc_path}/7bowtie2/*stats.txt; do
    # Extract filename without extension
    filename=\$(basename "\$file" _stats.txt)
    # Extract coverage information
    alignment=\$(sed -n '\$p' \$file)
    # Append to the output file
    echo -e "\$filename\t\$alignment" >> "\$output_file"
done

chmod -R 777 ${qc_path}

echo ''
echo 'DONE WITH RUNNING QC-ANALYSIS. ENJOY THE END PRODUCT AND BE THANKFUL.'
echo ''


conda deactivate
EOF
)

# Create a temporary SLURM job script
temp_job_script="running_qc_analysis.sb"
echo "$job_script" > "$temp_job_script"

# Submit the SLURM job
sbatch "$temp_job_script"

