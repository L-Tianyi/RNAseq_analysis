# RNAseq_analysis

The bam files are generated with STAR alignment tool. 
The number of reads per gene is counted by Salmon tool. 

Deseq2 and IsoformSaitchAnalyzeR take the output of Salmon as an imput,
Whereas DEXSeq and majiq take the bamfiles from STARalignment as an imput. 

Deseq2 analyzes differential genen expression level.
IsoformSwitchAnalyzeR, DEXSeq, and majiq are used to test alternative splicing event among the conditions.

