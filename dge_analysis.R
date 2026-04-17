# Differential Gene Expression
library(edgeR)


#Load Data
geneSymbolMap=read.table("data/ENSG_ID2Name.txt")
colnames(geneSymbolMap)=c("geneName", "geneSymbol") #Add column names to gene name/symbol table

rawCounts=read.table("data/allcounts_final.tsv", header=TRUE, stringsAsFactors=TRUE, row.names=1) #raw counts masterlist for 6 samples


#Edit Raw Counts File
rawCounts2=head(rawCounts, -5) #exclude last 5 rows containing information about reads/pairs

quant=apply(rawCounts2,1,quantile,0.75) #upper quartile of raw counts
keep=which((quant >= 25) == 1) #filter to keep read counts above 25%
rawCountsFilter=rawCounts2[keep,]


#DGE List Object
classLabels=factor(c( rep("untreated",3), rep("treated",3) )) #add untreated vs treated labels to samples

genesData=rownames(rawCountsFilter)
geneNames=geneSymbolMap$geneSymbol[match(genesData, geneSymbolMap$geneName)] #replace ensembl gene ID with gene names

dgeObj=DGEList(counts=rawCountsFilter, genes=geneNames, group=classLabels) #DGE object with filtered counts, gene names, and grouping labels
dgeObj=calcNormFactors(dgeObj) #normalize DGE object


#Q4: MDS Plot
plotMDS(dgeObj, main='MDS Plot of Gene Expression')


#Differentially Expressed Genes
dgeObj=estimateCommonDisp(dgeObj, verbose=TRUE) #common dispersion before exact test
dgeObj=estimateTagwiseDisp(dgeObj) #tagwise dispersion before exact test

dgeET=exactTest(dgeObj, pair=c('untreated','treated'))
dgeTop10=topTags(dgeET, n=10, sort.by='logFC') #top 10 by absolute log-fold change
dgeTop10


#Q6: Cut-Off
dgeTop=topTags(dgeET, sort.by='logFC')
dgeTopCut=dgeTop$table[abs(dgeTop$table$logFC) >= 1, ] #absolute log2 fold change cut-off of >=1
dim(dgeTopCut)
dgeTopCut

dgeObj=decideTests(dgeET, adjust.method = "BH", p.value = 0.05) #BH statistical threshold of p-value=0.05
summary(dgeObj)
dgeETSig = dgeET$table[as.logical(dgeObj[, 1]),] #filtered for significant genes
dim(dgeETSig)

dgeMatrix = data.frame("Gene"= rownames(dgeETSig), #reformatting into dataframe with relevant information
            "Gene_Name"= geneSymbolMap$geneSymbol[match(rownames(dgeETSig), geneSymbolMap$geneName)],
            "Log10_Pvalue"= log10(dgeETSig$PValue),
            "Log_fold_change"= dgeETSig$logFC)

head(dgeMatrix)
dim(dgeMatrix) #double check number of genes (or rows)

write.csv(dgeMatrix, "data/dge_matrix.csv")
#citation('edgeR')
#packageVersion('edgeR')