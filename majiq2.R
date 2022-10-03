#####intallation#######
conda install htslib
conda install gcc
conda install python=3.8

conda install -c conda-forge gxx

python3 -m pip install git+https://bitbucket.org/biociphers/majiq_academic.git
# pip install git+https://bitbucket.org/biociphers/majiq_academic.git#majiq
# conda install -c http://majiq.biociphers.org/download/channel majiq


# download htslib 1.13 archive to current directory
curl -OL https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2
tar -xf htslib-1.13.tar.bz2  # extract archive
cd htslib-1.13  # change into htslib source directory
# configure, make, and install htslib to ~/install/htslib-1.13
./configure --prefix=/your/working/directory/install/htslib-1.13
make
make install

#Install MAJIQ using pip. MAJIQ setup needs to know where HTSLib is installed.
# change to where library/include directories are installed
export HTSLIB_LIBRARY_DIR=/your/path/to/libexec/htslib
export HTSLIB_INCLUDE_DIR=/your/path/to/include/htslib
# NOTE: from same directory as this README
pip install git+https://bitbucket.org/biociphers/majiq_academic.git  # install both majiq and voila
#############################


###dowload a GFF3 annotation file#
#ensembl GFF3 is used
curl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.primary_assembly.annotation.gff3.gz -o gencode.v40.primary_assembly.annotation.gff3.gz
gunzip gencode.v40.primary_assembly.annotation.gff3.gz


####analysis#######

##build a configuration file
cd /your/working/directory
touch config.txt
vim config.txt
[info]
bamdirs=/your/bam/file/directory
genome=hg38
[experiments]
Control=1.Aligned.sortedByCoord.out,2.Aligned.sortedByCoord.out,3.Aligned.sortedByCoord.out,4.Aligned.sortedByCoord.out
condition1=5.Aligned.sortedByCoord.out,6.Aligned.sortedByCoord.out,7.Aligned.sortedByCoord.out
condition2=8.Aligned.sortedByCoord.out,9.Aligned.sortedByCoord.out,10.Aligned.sortedByCoord.out

##Run majiq builder
majiq build /path/to/ggf3 -c /path/to/configfile -j 4 -o /output/path

##psi quantifier
majiq psi -o /output/path -n condition1 /mnt/data/tianyi/CHED/majiq/build/Nb14Rx2Aligned.sortedByCoord.out.majiq 

##delta-PSI quantifier
#quantifies the differential splicing between two different groups (or conditions)
majiq deltapsi -grp1 /path/to/1.Aligned.sortedByCoord.out.majiq /path/to/2.Aligned.sortedByCoord.out.majiq -grp2 /path/to/3.Aligned.sortedByCoord.out.majiq /path/to/4.Aligned.sortedByCoord.out.majiq -j 4 -o /output/path -n Control Condition1


##Visualize results with Interactive VOILA. Voila command:
voila view /path/to/your/build/splicegraph.sql /path/to/deltapsi/output/Control-Condition1.deltapsi.voila

##generate VOILA tsv file
#In order to use MAJIQ/VOILA output as an input for post-quantification analysis you can generate VOILA tsv file. 
#This file is a tab separated value file where each line is a single LSV. As an example, some of the columns are gene id, lsv id, 
#Expected PSI or DeltaPSI, confidence or junction quantification, and many others. Thiis table can be filtered and ordered by different attributes (columns).
voila tsv /path/to/your/build/splicegraph.sql path/to/deltapsi/output/Control-Condition1.deltapsi.voila -f name_of_output_file.tsv

