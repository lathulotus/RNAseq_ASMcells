## RNAseq Data Processing Pipeline
Run RNAseq data processing pipeline to ascertain differentially expressed transcripts. Data from [Himes et al., PLOS one, 2014](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0099625#s2).

### Background
Primary human airway smooth muscle (ASM) left untreated or treated with dexamethasone (for inflammatory disease treatment). RNA isolated and prepped using llumina TruSeq RNA Library Prep Kit v2. Sequencing of 75 bp paired-end reads was performed with an Illumina HiSeq 2000 instrument. The first 12 bases of all reads were trimmed by the authors with the FASTX Toolkit due to sequence bias in the initial 12 bases on the 5ʹ end of reads.

### Aim
- Use HISAT2 to obtain quality control metrics and assess quality of the samples
- Use HTSeq to perform counting and normalization
- Use edgeR to perform differential gene expression and visualize data
- Create multidimensional scaling (MDS) plot to analyze findings
- Ascertain top 10 most differentially gene expression
- Run MISO pipeline to assess alternatively expressed isoforms