##install HTSeq to count the number of reads:
conda creare htseq
conda activate htseq
conda install -c bioconda htseq

##Move the .py files from DEXseq R package to your working directory
#to get the original path of your .py files. go to /DEXseq.R
scp /original/path/of/dexseq_count.py /your/working/directory/DEXseq
scp /original/path/of/dexseq_prepare_annotation.py /your/working/directory/DEXseq

##Create DEXseq gff file
python dexseq_prepare_annotation.py -r no /path/to/your/.gtf output/path/and/name/of/your/.gff

##create a list of samples
cd /your/working/directory/DEXseq
touch sample.list
vim sample.list

#build a loop
for f in `cat sample.list`; 
do echo $f \ 
python dexseq_count.py -p yes -r pos -s reverse -f bam /path/to/your/.gff /path/to/your/.bam output/path/sample_name.txt 
done
#tips: do not change the line for the python code, otherwise it doesn't work

#The quote "" lead to DEXseq not recognise the file
#so, have to remove them (https://support.bioconductor.org/p/9143537/#9143593)
#to do so:
for f in `cat sample.list`;
do sed 's/\"//g' ${f}.txt > ${f}_clean.txt
done




