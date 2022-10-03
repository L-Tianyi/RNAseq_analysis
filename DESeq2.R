###installation of DESeq2
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")

#install tximeta
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("tximeta")



####built a sample table with sample id, conditions, and file path
#location of quantification file: /home/tianyi/salmon/quants2

sample <- data.frame(names=sampleid, condition=conditions)

filepath <- rep(0,nrow(sample))
name <- rep(0,nrow(sample))
for (i in 1:length(filepath)) {
    name[i] <- sample[i,1]
    filepath[i] <- paste0("path/to/your/quantificationfiles",name[i],"_quant/quant.sf")
}

sample$files <- filepath



###import in transcriptome-quantification data with tximeta/tximport
library(tximeta)

se <- tximeta(sample)


###SummarizedExperiment output
############################################################################################################
#The SumarizedExperiment object consists of three components: assay (matrix of counts), rowRanges (genomic #
#ranges), colData (sample informations).                                                                   #
#We can investigate this SummarizedExperiment object by looking at the matrices in the assays slot, the    #
#phenotypic data about the samples in colData slot, and the data about the genes in the rowRanges slot.    #
############################################################################################################
suppressPackageStartupMessages(library(SummarizedExperiment))
colData(se)  #check the columns
assayNames(se)
rowRanges(se)
seqinfo(se) #show the metadata, make sre we have appropriate genome information



###Summarization to gene-level
gse <- summarizeToGene(se)
rowRanges(gse)
head(assay(gse), 3)#to view your results


###to relevel the conditions
library("magrittr")
#e.g:
gse$condition <- factor(gse$condition, levels = c("ASOControl","Control","ASO"))
gse$condition
#  [1] ASO        ASOControl Control    ASO        ASOControl Control
#  [7] ASO        ASOControl Control    ASO        ASOControl Control
# Levels: ASOControl Control ASO


###start analyzing
library("DESeq2")
dds <- DESeqDataSet(gse, design = ~ condition)  #generate DESeqDataSet object
##########################for paired data#######################################
# Yes, you should use a multi-factor design which includes the sample 
# information as a term in the design formula. This will account for 
# differences between the samples while estimating the effect due to 
# the condition. The condition of interest should go at the end of 
# the design formula, e.g. ~ subject + condition.
dds <- DESeqDataSet(gse, design = ~ cell + condition)
################################################################################

dds <- DESeq(dds) 

##Extract results
res <- results(dds)
# we could have specified the coefficient or contrast we want to build a results table for, 
#using the following  commands:                           
res1 <- results(dds, contrast=c("condition","ASO","Control"))
res1


####Add gene-symbol column to the data frame
BiocManager::install("biomaRt")
library(biomaRt)

#Step1: Identifying the database you need
listEnsembl()
ensembl <- useEnsembl(biomart = "genes")
#Step 2: Choosing a dataset
searchDatasets(mart = ensembl, pattern = "hsapiens")
#                 dataset              description    version
#80 hsapiens_gene_ensembl Human genes (GRCh38.p13) GRCh38.p13
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl) 
#↑To use a dataset we can update our Mart object using the function useDataset()

##Or just in on step:
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 105)

##build a biomaRT query
searchAttributes(mart = ensembl, pattern = "symbol")
#                name                description         page
#65       hgnc_symbol                HGNC symbol feature_page
#98 uniprot_gn_symbol UniProtKB Gene Name symbol feature_page
searchFilters(mart = ensembl, pattern = "ensembl.*id")
#                            name
#51               ensembl_gene_id
#52       ensembl_gene_id_version
#53         ensembl_transcript_id
#54 ensembl_transcript_id_version
#55            ensembl_peptide_id
#56    ensembl_peptide_id_version
#57               ensembl_exon_id
#                                                      description
#51                       Gene stable ID(s) [e.g. ENSG00000000003]
#52       Gene stable ID(s) with version [e.g. ENSG00000000003.15]
#53                 Transcript stable ID(s) [e.g. ENST00000000233]
#54 Transcript stable ID(s) with version [e.g. ENST00000000233.10]
#55                    Protein stable ID(s) [e.g. ENSP00000000233]
#56     Protein stable ID(s) with version [e.g. ENSP00000000233.5]
#57                              Exon ID(s) [e.g. ENSE00000000003]

#build an array of sample names
geneid <- rownames(res1)
####################THIS step is IMPORTANT#########################
geneidtrim = sapply(strsplit(geneid, ".", fixed=T), function(x) x[1])
###################################################################

gene_symbol<- getBM(attributes = c("ensembl_gene_id","hgnc_symbol", "description"),
                  filters = "ensembl_gene_id", 
                  values = geneidtrim, 
                  mart = ensembl)
                  
#The BiomaRt out put may have different number of columns to your deqseq reults so need to check and solve it accordingly
#Check duplication
gene_symbol[duplicated(gene_symbol$ensembl_gene_id),]

#Check missingness
geneidtrim[!geneidtrim %in% gene_symbol$ensembl_gene_id]

#merge the duplicated rows by data.table:
#Step1. group gene symbols by ensemble ids
#step2. Convert the element in each list of the gene_symbol into a single string, 
#       then change the gene_symbol column into a charachter vector and 
dt <- as.data.table(gene_symbol) #convert to a data table
dt <- dt[ , .(hgnc_symbol = list(hgnc_symbol), description = unique(description)), by = ensembl_gene_id][ ,
      hgnc_symbol := unlist(lapply(hgnc_symbol, function(x) { paste(unlist(x),collapse = ";")}))]

resdt <- as.data.table(res)
ensembl_gene_id <- geneidtrim
resdt <- cbind(ensembl_gene_id,resdt)

anno <- dt[ensembl_gene_id, on = .(ensembl_gene_id = ensembl_gene_id)] 

resdt <- resdt[anno, on = .(ensembl_gene_id = ensembl_gene_id)]


###Exporting results to CSV files
write.csv(as.data.frame(resdt), 
          file="res1.csv",row.names = FALSE)
         





