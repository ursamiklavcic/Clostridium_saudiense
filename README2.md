\# 18.11.2025

\# I have 5 hybrid assemblies H005\_73, Hoo5\_98, H005\_118, H007\_176 and H009\_221 with v10 chemistry

\# Others are short-read only assemblies



\## ANI score with pyANI-plus anim

srun

mkdir /home/nlzoh.si/ursmik1/projects/C\_saudiense/002\_ANI/



pyani-plus anim ~/projects/C\_saudiense/001\_all\_fasta/ --database ~/projects/C\_saudiense/002\_ANI/csaudiense.db --create-db --name 'Clostridium saudiense'



pyani-plus anib ~/projects/C\_saudiense/001\_all\_fasta/ --database ~/projects/C\_saudiense/002\_ANI/csaudiense.db --name 'Clostridium saudiense'







\## SNIPPY + GUBBINS + MEGA 12

\# Interactive session

srun --time=1-24:00 --mem=50G --cpus-per-task=20 --job-name=snippy --pty /bin/bash



conda activate snippy

mkdir 003\_snippy

cd 003\_snippy/



snippy-multi reads.tab --ref ../001\_all\_fasta/H005\_98\_hybrid.fasta --cpus 20 > runme.sh

less runme.sh



sh ./runme.sh



mkdir -p gubbins/

cd gubbins

conda deactivate



source /home/nlzoh.si/larbez1/miniconda3/bin/activate /home/nlzoh.si/larbez1/miniconda3/envs/snp\_env



run\_gubbins.py -c 20 --prefix gubb ../core.full.aln



\# snp-sites -c -o ./Csaudiense\_clean.core.fasta gubbins/gubb.filtered\_polymorphic\_sites.fasta

\# snp-dists ./Csaudiense\_clean.core.fasta > snippy\_gubb\_dists.tab



conda deactivate



\# plot\_snippy.R

\# results in out/snippy







\## PROKKA + ROARY + PHANDANGO



source /home/nejc/miniconda3/bin/activate /home/nejc/miniconda3/envs/prokka



cd ~/projects/C\_saudiense/001\_all\_fasta



for k in \*.fasta; do prokka $k --outdir ../004\_prokka/"$k".prokka.output --prefix PROKKA\_$k; echo $k; done



conda deactivate



source /home/nlzoh.si/leomar1/miniconda3/bin/activate /home/nlzoh.si/leomar1/miniconda3/envs/roary



cd 004\_prokka

mkdir GFF

cp ./\*/\*.gff GFF



for file in ./PROKKA\_\*.fasta.gff; do

&nbsp; newname="${file/PROKKA\_/}"

&nbsp; mv "$file" "$newname"

done



roary -n -e -p 18 ./GFF/\*.gff -f roary



conda deactivate



\# Load into phandango with metadata etc.

# 

\## DRAM anotation 

005\_DRAM/dram.sb



