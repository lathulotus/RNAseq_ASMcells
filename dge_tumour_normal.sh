#Quality Control
salloc -N1 -t1:00:00 --ntasks=8
module load java fastqc
fastqc *.fastq.gz

____________________________________________________________
#Align the data with HISAT2
salloc -N1 -t1:00:00 --ntasks=8
module load hisat2

#Reference Genome & Annotations (GRCh37)
#http://ftp.ensembl.org/pub/grch37/release-113/fasta/homo_sapiens/dna/
wget http://ftp.ensembl.org/pub/grch37/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa.gz
wget http://ftp.ensembl.org/pub/grch37/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna_sm.chromosome.16.fa.gz
wget http://ftp.ensembl.org/pub/grch37/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz

gunzip *.gz

#Indexing
#https://daehwankimlab.github.io/hisat2/download/
wget https://genome-idx.s3.amazonaws.com/hisat/grch37_snptran.tar.gz

gunzip *.gz | tar -xvf #unzip & untar

#Alignments
hisat2 -p 8 --rg-id=N052611_Dex --rg SM:N052611_Dex --rg PL:ILLUMINA -x /rna_ref_index/grch37_snp_tran/genome_snp_tran --rna-strandness RF -1 /rawSamples/N052611_Dex_r1.fastq.gz -2 /rawSamples/N052611_Dex_r2.fastq.gz -S ./N052611_Dex.sam

nano sbatch_hisat.sh
#SBATCH CODE-----------------------------
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=4:00:00
#SBATCH --job-name A781_hisat2_alignments
module load hisat2
hisat2 -p 8 --rg-id=N052611_Dex --rg SM:N052611_Dex --rg PL:ILLUMINA -x $AHOME/rna_ref_index/grch37_snp_tran/genome_snp_tran --rna-strandness RF -1 $AHOME/rawSamples/N052611_Dex_r1.fastq.gz -2 $AHOME/rawSamples/N052611_Dex_r2.fastq.gz -S $AHOME/hisat2_alignments/N052611_Dex.sam
hisat2 -p 8 --rg-id=N052611_untreated --rg SM:N052611_untreated --rg PL:ILLUMINA -x $AHOME/rna_ref_index/grch37_snp_tran/genome_snp_tran --rna-strandness RF -1 $AHOME/rawSamples/N052611_untreated_r1.fastq.gz -2 $AHOME/rawSamples/N052611_untreated_r2.fastq.gz -S $AHOME/hisat2_alignments/N052611_untreated.sam
hisat2 -p 8 --rg-id=N080611_Dex --rg SM:N080611_Dex --rg PL:ILLUMINA -x $AHOME/rna_ref_index/grch37_snp_tran/genome_snp_tran --rna-strandness RF -1 $AHOME/rawSamples/N080611_Dex_r1.fastq.gz -2 $AHOME/rawSamples/N080611_Dex_r2.fastq.gz -S $AHOME/hisat2_alignments/N080611_Dex.sam
hisat2 -p 8 --rg-id=N080611_untreated --rg SM:N080611_untreated --rg PL:ILLUMINA -x $AHOME/rna_ref_index/grch37_snp_tran/genome_snp_tran --rna-strandness RF -1 $AHOME/rawSamples/N080611_untreated_r1.fastq.gz -2 $AHOME/rawSamples/N080611_untreated_r2.fastq.gz -S $AHOME/hisat2_alignments/N080611_untreated.sam
hisat2 -p 8 --rg-id=N61311_Dex --rg SM:N61311_Dex --rg PL:ILLUMINA -x $AHOME/rna_ref_index/grch37_snp_tran/genome_snp_tran --rna-strandness RF -1 $AHOME/rawSamples/N61311_Dex_r1.fastq.gz -2 $AHOME/rawSamples/N61311_Dex_r2.fastq.gz -S $AHOME/hisat2_alignments/N61311_Dex.sam
hisat2 -p 8 --rg-id=N61311_untreated --rg SM:N61311_untreated --rg PL:ILLUMINA -x $AHOME/rna_ref_index/grch37_snp_tran/genome_snp_tran --rna-strandness RF -1 $AHOME/rawSamples/N61311_untreated_r1.fastq.gz -2 $AHOME/rawSamples/N61311_untreated_r2.fastq.gz -S $AHOME/hisat2_alignments/N61311_untreated.sam
#SBATCH CODE-----------------------------
sbatch sbatch_hisat.sh

____________________________________________________________
#count raw reads with htseq-count
salloc -N1 -t1:00:00 --ntasks=8
module load samtools

module load python/3
virtualenv --no-download $HOME/htseq
source $HOME/htseq/bin/activate
pip install --noindex --upgrade pip
pip install --no-index htseq

#SAM to BAM
samtools sort -@ 8 -n -o ./N052611_Dex.bam $HOME/Assignment7/hisat2_alignments/N052611_Dex.sam
samtools sort -@ 8 -n -o ./N052611_untreated.bam $HOME/Assignment7/hisat2_alignments/N052611_untreated.sam
samtools sort -@ 8 -n -o ./N080611_Dex.bam $HOME/Assignment7/hisat2_alignments/N080611_Dex.sam
samtools sort -@ 8 -n -o ./N080611_untreated.bam $HOME/Assignment7/hisat2_alignments/N080611_untreated.sam
samtools sort -@ 8 -n -o ./N61311_Dex.bam $HOME/Assignment7/hisat2_alignments/N61311_Dex.sam
samtools sort -@ 8 -n -o ./N61311_untreated.bam $HOME/Assignment7/hisat2_alignments/N61311_untreated.sam

#HTSeq Count
htseq-count --format bam --order name --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $HOME/Assignment7/hisat2_alignments/N052611_Dex.bam $HOME/Assignment7/ref37/Homo_sapiens.GRCh37.87.gtf > $HOME/Assignment7/htseq_count/N052611_Dex.tsv

nano sbatch_htseq.sh
#SBATCH CODE--------------------------------
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=4:00:00
#SBATCH --job-name A781_htseq_count
# DEFINE NECESSARY VARIABLES
source $HOME/htseq/bin/activate
htseq-count --format bam --order name --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $HOME/Assignment7/hisat_bams/N052611_Dex.bam $HOME/Assignment7/ref37/Homo_sapiens.GRCh37.87.gtf > $HOME/Assignment7/htseq_count/N052611_Dex.tsv
htseq-count --format bam --order name --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $HOME/Assignment7/hisat_bams/N052611_untreated.bam $HOME/Assignment7/ref37/Homo_sapiens.GRCh37.87.gtf > $HOME/Assignment7/htseq_count/N052611_untreated.tsv
htseq-count --format bam --order name --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $HOME/Assignment7/hisat_bams/N080611_Dex.bam $HOME/Assignment7/ref37/Homo_sapiens.GRCh37.87.gtf > $HOME/Assignment7/htseq_count/N080611_Dex.tsv
htseq-count --format bam --order name --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $HOME/Assignment7/hisat_bams/N080611_untreated.bam $HOME/Assignment7/ref37/Homo_sapiens.GRCh37.87.gtf > $HOME/Assignment7/htseq_count/N080611_untreated
htseq-count --format bam --order name --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $HOME/Assignment7/hisat_bams/N61311_Dex.bam $HOME/Assignment7/ref37/Homo_sapiens.GRCh37.87.gtf > $HOME/Assignment7/htseq_count/N61311_Dex.tsv
htseq-count --format bam --order name --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $HOME/Assignment7/hisat_bams/N61311_untreated.bam $HOME/Assignment7/ref37/Homo_sapiens.GRCh37.87.gtf > $HOME/Assignment7/htseq_count/N61311_untreated.tsv
#SBATCH CODE--------------------------------
sbatch sbatch_htseq.sh

#Merge for EdgeR
#join N052611_untreated.tsv N052611_Dex.tsv | join - N080611_untreated.tsv | join - N080611_Dex.tsv | join - N61311_untreated.tsv | join - N61311_Dex.tsv > allCounts.tsv
#echo "GeneID N052611_untreated N052611_Dex N080611_untreated N080611_Dex N61311_untreated N61311_Dex" > header.txt
#cat header.txt allCounts.tsv > allCountsFinal.tsv

join N052611_untreated.tsv N080611_untreated.tsv | join - N61311_untreated.tsv | join - N052611_Dex.tsv | join - N080611_Dex.tsv | join - N61311_Dex.tsv > allCounts.tsv
echo "GeneID N052611_untreated N080611_untreated N61311_untreated N052611_Dex N080611_Dex N61311_Dex" > header.txt
cat header.txt allCounts.tsv > allCountsFinal.tsv

____________________________________________________________
#Differential expression with edgeR
# In R/Rmd file

____________________________________________________________
#Loop through 6 samples for HISAT2

for SAMPLE in ${SAMPLENAMES[@]}
do
  SAMPLENAME=$(echo $SAMPLE | sed 's/_r1.fastq.gz//')
  echo "Running HISAT2 for $SAMPLENAME"
  R1="${SAMPLENAME}_r1.fastq.gz"
  R2="${SAMPLENAME}_r2.fastq.gz"
  OUT="${SAMPLENAME}.sam"
  hisat2 -p 8 --rg-id=$SAMPLENAME --rg SM:$SAMPLENAME --rg PL:ILLUMINA -x $REF --rna-strandness RF -1 $SAMPLES/$R1 -2 $SAMPLES/$R2 -S $HISATOUT/$OUT
done

