# Download the SRA Files using the below script or use it as a job script

```bash
base=your/HPC/working/directory
srr="SRR21948648 SRR21948649 SRR21948650 SRR21948651 SRR21948652 SRR21948653 SRR21948654 SRR21948655 SRR21948656 SRR21948657 SRR21948658 SRR21948659 SRR21948660 SRR21948661 SRR21948662 SRR21948663 SRR21948664 SRR21948665 SRR21948666 SRR21948667 SRR21948668 SRR21948669 SRR21948670 SRR21948671 SRR21948672 SRR21948673 SRR21948674 SRR21948675 SRR21948676 SRR21948677"
cd $base

mkdir -p $base/fastq/raw_reads
module load sra-toolkit

for s in $srr; do
        prefetch $s -O sra_files
        fasterq-dump sra_files/$s -O fastq/raw_reads
done
```
# Download these reference genomes and other necessary files before running the pipeline. B

**But these references are specific to these samples used in the pipeline, always check the reference that the paper used for alignment before choosing your alignment**

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz -P reference/
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.knownGene.gtf.gz -P reference/
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes -P reference/
wget https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt -P reference/     # Motif DB for Homer/MEME-Chip

# https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg19-blacklist.v2.bed.gz
# Download the above file to your local and upload to HPC using scp
# scp your/local/working/directory/hg19-blacklist.v2.bed.gz your.username@gmail.com:/your/HPC/working/directory

gunzip reference/*.gz
```

