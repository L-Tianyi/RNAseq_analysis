###DEXseq document: https://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html#requirements-on-gtf-files

###install Dexseq
BiocManager::install("DEXSeq")

###find the directory of .py files in DEXseq package
pythonScriptsDir = system.file( "python_scripts", package="DEXSeq" )
list.files(pythonScriptsDir)
# [1] "dexseq_count.py"              "dexseq_prepare_annotation.py"
system.file( "python_scripts", package="DEXSeq", mustWork=TRUE )

##move to DEXseq.sh, run python script for：
#(by "dexseq_prepare_annotation.py") creation of .gff 
#(by "dexseq_count.py" ) and count the number of reads 


##when we get the count files, back to R


###reading data into R
setwd("/working/directory/DEXseq")
inDir = "/working/directory/DEXseq/counts"
countFiles = list.files(inDir, pattern="_clean.txt$", full.names=TRUE)
basename(countFiles)
flattenedFile = "/path/to/.gff"


###Now, construct an DEXSeqDataSet object from this data
library( "DEXSeq" )

##if you have more than two conditions, comparisons are made for each pair one by one. Just made a subset of the paths to the sample files that you want to compare.

##Create a sample table: 
library(data.table)
sample1 <- data.frame(row.names=sampleid, condition=conditions)
#sample names as row names and the conditions that vary between samples

countFile_subset = countFiles[1:7]

dxd1 = DEXSeqDataSetFromHTSeq(
   countFile_subset,
   sampleData=sample1,
   design= ~ sample + exon + condition:exon,
   flattenedfile=flattenedFile )

###Normalization
dxd1 = estimateSizeFactors( dxd1 )

###Dispersion estimation
dxd1 = estimateDispersions( dxd1 )
plotDispEsts( dxd1 )


###Testing for differential exon usage
#The function testForDEU() performs these tests for each exon in each gene.
dxd1 = testForDEU( dxd1 )
#Estimate fold changes
dxd1 = estimateExonFoldChanges( dxd1, fitExpToVar="condition")


###call results
#summarize the results without showing the values from intermediate steps
#The result is a DEXSeqResults object, which is a subclass of a DataFrame object.
dxr1 = DEXSeqResults( dxd1 )
dxr1

#The description of each of the column of the object DEXSeqResults can be 
#found in the metadata columns.
mcols(dxr1)$description


####Just For paired-sample Analysis###################################
#####add additional variables to the test#############################
formulaFullModel    =  ~ sample + exon + cell:exon + condition:exon
formulaReducedModel =  ~ sample + exon + cell:exon
# you need to pass the formulae to both estimateDispersions and testForDEU
dxd = estimateDispersions( dxd, formula = formulaFullModel )
dxd = testForDEU( dxd, 
    reducedModel = formulaReducedModel, 
        fullModel = formulaFullModel )
######################################################################


###Visualization
plotDEXSeq( dxr1, "ENSG00000000001", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

#you can make plots for the top significant genes
dxr1ordr <- dxr1[order(dxr1$padj),]
head(dxr1ordr)
topgenes <- dxr1ordr$groupID[1:10]
for (i in 1:10){
   pdf( file = paste0(topgenes[i],".pdf"))
   plotDEXSeq( dxr1, topgenes[i], legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
   dev.off()
}


