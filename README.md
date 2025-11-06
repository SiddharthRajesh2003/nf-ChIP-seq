Download the SRA Files and convert to FASTQ with fasterq-dump
```bash
module avail sra-toolkit
module load sra-toolkit

mkdir -p sra_files
mkdir -p fastq/raw_reads
mkdir -p fastq/trimmed_reads

prefetch prefetch SRR31058337 SRR31058336 SRR31058335 SRR31058334 SRR31058333 SRR31058332 -O sra_files/
```