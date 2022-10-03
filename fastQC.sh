#################
#####FastQC######
#################

###Installation
conda activate myenv
#check is java installed
java -version
#install
conda install -c bioconda fastqc

###Run
# https://raw.githubusercontent.com/s-andrews/FastQC/master/INSTALL.txt

# To run FastQC non-interactively you should use the fastqc wrapper script to launch
# the program.  You will probably want to use the zipped install file on every platform
# (even OSX).

# To run non-interactively you simply have to specify a list of files to process
# on the commandline:

fastqc somefile.txt someotherfile.txt

# You can specify as many files to process in a single run as you like.  If you don't
# specify any files to process the program will try to open the interactive application
# which may result in an error if you're running in a non-graphical environment.

# There are a few extra options you can specify when running non-interactively.  Full
# details of these can be found by running 

fastqc --help

# By default, in non-interactive mode FastQC will create an HTML report with embedded
# graphs, but also a zip file containing individual graph files and additional data files
# containing the raw data from which plots were drawn.  The zip file will not be extracted
# by default but you can enable this by adding:

--extract

# To the launch command.

# If you want to save your reports in a folder other than the folder which contained
# your original FastQ files then you can specify an alternative location by setting a
# --outdir value:

--outdir=/some/other/dir/

# If you want to run fastqc on a stream of data to be read from standard input then you
# can do this by specifing 'stdin' as the name of the file to be processed and then 
# streaming uncompressed fastq format data to the program.  For example:

zcat *fastq.gz | fastqc stdin

# If you want the results from a streamed analysis sent to a file with a name other than
# stdin then you can add a colon and put the file name you want, for example:

zcat *fastq.gz | fastqc stdin:my_results

# ..would write results to my_result.html and my_results.zip.

##create a loop
cd /your/working/directory

for f in `cat sample.list`
do echo $f 
zcat ${f}.fq.gz | fastqc stdin:${f} --outdir=/your/output/path
done 





