####installation
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
BiocManager::install("IsoformSwitchAnalyzeR")


#####Import data from Salmon via tximeta
library(data.table)
library(IsoformSwitchAnalyzeR)

sample <- fread("/path/to/my/samplelist/SampleList.csv")#sample list with names, files, and condition
aSwitchList <- importSalmonData(
    salmonFileDataFrame = sample, #as created above
    showProgress=FALSE,           # For nicer printout in vignette
)

summary(aSwitchList)


#################This below doesn't work for me so I did each step separately (shown after this block)#######################
aSwitchList <- isoformSwitchAnalysisPart1(
    switchAnalyzeRlist = aSwitchList,
    pathToGTF= '/your/path/to/gtf', #should be v33, the name is wrong but the gtf is release 33
    pathToOutput = '/your/output/directory',
    outputSequences = TRUE, # change to TRUE when analyzing your own data 
    prepareForWebServers = FALSE # change to TRUE if you will use webservers for external sequence analysis
    ) 

extractSwitchSummary( aSwitchList ) #The number of switching features
###############################################################################################################################


####Filtering
SwitchListFiltered <- preFilter(
    switchAnalyzeRlist = aSwitchList,
    geneExpressionCutoff = 1,
    isoformExpressionCutoff = 0,
    removeSingleIsoformGenes = TRUE
)

summary(SwitchListFiltered)


####Identifying Isoform Switches
#Testing for Isoform Switches via DEXSeq (The SLOW step!!)
aSwitchListAnalyzed <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = SwitchListFiltered,
    reduceToSwitchingGenes=TRUE
)

######################################for paired data#######################
# The ability to correct for confounding factors rely on that these 
# are annotated in the design matrix (as supplied when creating the 
# switchAnalyzeRlist) and stored in the switchAnalyzeRlist (accessed 
# via exampleSwitchList$designMatrix ). The design matrix is a 
# data.frame with two mandatory columns: sampleID and condition but 
# enables incorporation of additional cofactors in the experimental 
# setup (such as batch effects) - simply add these to the design 
# matrix as additional columns. If additional co-factors have been 
# added isoformSwitchTestDEXSeq will take that into account if 
# correctForConfoundingFactors=TRUE #(default). 
# The co-factor corrected effect sizes (IF and dIF values) can be 
# added to the switchAnalyzeRlist overwriting the old dIF and IF 
# values as controlled by the overwriteIFvalues argument (default is TRUE).
cell <- rep(c("000","111","222","333"),each = 3)
SwitchListFiltered$designMatrix$cell = cell
aSwitchListAnalyzed <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = SwitchListFiltered,
    reduceToSwitchingGenes=TRUE,
    correctForConfoundingFactors=TRUE
)
##############################################################################

extractSwitchSummary(aSwitchListAnalyzed)


####Analyse open reading frames
"orfAnalysis" %in% names( aSwitchListAnalyzed )
#[1] FALSE
##Analyzing Known and Novel Isoforms
aSwitchListAnalyzed <- addORFfromGTF( aSwitchListAnalyzed , pathToGTF = '/home/tianyi/IsoformSwitch/gencode.v39.annotation.gtf.gz') 
#Step 1 of 2: importing GTF (this may take a while)...
#Step 2 of 2: Adding ORF...
#    Added ORF info (incl info about isoforms annotated as not having an ORF) to 3036 isoforms.
#        This correspond to 100% of isoforms in the switchAnalyzeRlist.
#            Which includes 100% of isoforms from annotated genes (novel isoforms not counted) in the switchAnalyzeRlist.
#Done.
"orfAnalysis" %in% names( aSwitchListAnalyzed )
#[1] TRUE

aSwitchListAnalyzed <- analyzeNovelIsoformORF(
switchAnalyzeRlist = aSwitchListAnalyzed,
analysisAllIsoformsWithoutORF = TRUE
)


####Extracting Nucleotide and Amino Acid Sequences
aSwitchListAnalyzed <- extractSequence(
    aSwitchListAnalyzed, 
    pathToOutput = '/output/path/',
    writeToFile= TRUE 
)


####predict alternative splicing, Another type of annotation
aSwitchListAnalyzed <- analyzeAlternativeSplicing(
    switchAnalyzeRlist = aSwitchListAnalyzed,
    quiet=TRUE
)

table( aSwitchListAnalyzed$AlternativeSplicingAnalysis$IR ) #the no. of Intron retention for example
table( aSwitchListAnalyzed$AlternativeSplicingAnalysis$ES ) #exon skipping
#   0    1    2    3    4    5    6
#1672  961  301   73   19    9    1
#Meaning 961 isoforms contain a single exon skipping (ES) and 301 isoform contain two ES.
table( aSwitchListAnalyzed$AlternativeSplicingAnalysis$MEE ) #mutually-exclusive exons
table( aSwitchListAnalyzed$AlternativeSplicingAnalysis$MES ) #multi-exon skipping
table( aSwitchListAnalyzed$AlternativeSplicingAnalysis$A5 )  #alternative 5' end
table( aSwitchListAnalyzed$AlternativeSplicingAnalysis$A3 )  #alternativ 3'end
table( aSwitchListAnalyzed$AlternativeSplicingAnalysis$ATSS ) #alternative transcription start sites 
table( aSwitchListAnalyzed$AlternativeSplicingAnalysis$ATTS ) #alternative transcription termination sites 



####predicting switch consequence
# the consequences highlighted in the text above
consequencesOfInterest <- c(
    'intron_retention',
    'NMD_status',
    'ORF_seq_similarity',
    'tss', 
    'tts',
    'last_exon',
    'isoform_seq_similarity',
    'ORF_genomic',
    '5_utr_seq_similarity',
    '3_utr_seq_similarity'
    )

aSwitchListAnalyzed <- analyzeSwitchConsequences(
    aSwitchListAnalyzed,
    consequencesToAnalyze = consequencesOfInterest, 
    dIFcutoff = 0.1, #default (0.1)
    showProgress=FALSE
)

extractSwitchSummary(aSwitchListAnalyzed, dIFcutoff = 0.1, filterForConsequences = TRUE)


####analyze Individual isoform switching
aSwitchListAnalyzedSubset1 <- subsetSwitchAnalyzeRlist(
    aSwitchListAnalyzed, 
    aSwitchListAnalyzed$isoformFeatures$condition_1 == 'control'
)


### Extract top switching genes (by q-value)
geneSwitch1 <-extractTopSwitches(
            aSwitchListAnalyzedSubset1, 
            filterForConsequences = FALSE, 
            n = Inf, 
            sortByQvals = TRUE
)


### Extract top 2 switching genes (by dIF values)
extractTopSwitches(
    aSwitchListAnalyzedSubset, 
    filterForConsequences = TRUE, 
    n = 2, 
    sortByQvals = FALSE
)






