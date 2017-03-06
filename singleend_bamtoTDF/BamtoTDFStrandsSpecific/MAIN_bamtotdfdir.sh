indir=/scratch/Users/allenma/scripts/singleend/testbam3/
outdir=${indir}/sortedbams/

#this part allows this script to use the other scripts in the directory
current_dir=$(pwd)
scriptdirname=$(dirname $0)
if [ $scriptdirname = '.' ]
then
scriptdirname="$current_dir"
fi
scriptdirname=${scriptdirname}/

#both of these files are organism specific
bedgraphgenomefile=${scriptdirname}human.hg19.genome
#if the bedgraphgenomefile does not exist for this genome you must create it. See program createbedtoolsgenome.py. You will need a fasta file for your organism. 
igvtoolsgenomefile=/opt/igvtools/2.1.24/genomes/hg19.genome
#if the igvtoolsgenomefile does not exist for this genome you must create it. This can be done using igv and a fasta file. Use the tab label genomes and create a genome. You will need a fasta file and hopefully gtf file (or a gff3). 


#title           :MAIN_bamtotdfdir.sh
#description     :This script take in a directory of bam files and make tdf files. All the bam files must be single end. This script uses several other scripts located in the same directory and submits things to tourque. 
#author          :Mary Ann Allen
#date            :20170214
#version         :0.1    
#usage           :First change the in directory above. Second, if you are not using hg19 change the genome files for igv and bedtools. Then "bash MAIN_bamtotdfdir.sh"
#notes           :Edit the in directory at the top of the script to be your directory of bam files.  All the bam files must be single end. This script uses several other scripts located in the same directory and submits jobs to tourque. Does not handle spliced mapped reads correctly.



bedgraphdir=${outdir}genomecoveragebed/
bedgraphfortdfdir=${outdir}genomecoveragebed/fortdf/
tdfdir=${outdir}genomecoveragebed/fortdf/tdfs/


mkdir -p $bedgraphdir
mkdir -p $bedgraphfortdfdir
mkdir -p $tdfdir

module load python_2.7.3
module load samtools_0.1.19 
module load bedtools2_2.22.0
module load igvtools_2.1.24

for pathandfilename in `ls $indir*.bam`; do
rootfilename=`basename $pathandfilename .bam`

qsub -v scriptdirname=$scriptdirname,indir=${indir},rootfilename=$rootfilename,outdir=$outdir,bedgraphdir=$bedgraphdir,bedgraphfortdfdir=$bedgraphfortdfdir,tdfdir=$tdfdir,bedgraphgenomefile=$bedgraphgenomefile,igvtoolsgenomefile=$igvtoolsgenomefile -N pipe${rootfilename} ${scriptdirname}bamtotdffile.sh

done


