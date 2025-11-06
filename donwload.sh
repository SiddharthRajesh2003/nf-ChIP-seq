#!/bin/bash

#SBATCH -J DNA-seq
#SBATCH -p general
#SBATCH -o sra.txt
#SBATCH -e sra.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@example.com
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --time=10:00:00
#SBATCH --mem=10GB
#SBATCH -A r00750

base=/N/project/Krolab/Siddharth/Personal/ChIP-seq
cd $base

module load sra-toolkit
mkdir $base/sra_files/






mkdir $base/reference
wget https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -C reference/
gunzip reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa human_ref.fa

