indir=/scratch/Users/allenma/scripts/singleend/testbam2/
outdir=${indir}/sortedbams/
current_dir=$(pwd)
scriptdirname=$(dirname $0)

if [ $scriptdirname = '.' ]
then
scriptdirname="$current_dir"
fi
scriptdirname=${scriptdirname}/

#organism specific
bedgraphgenomefile=${scriptdirname}human.hg19.genome
#if the bedgraphgenomefile does not exist for this genome you must create it. See program createbedtoolsgenome.py. You will need a fasta file. 
igvtoolsgenomefile=/opt/igvtools/2.1.24/genomes/hg19.genome
#if the igvtoolsgenomefile does not exist for this genome you must create it. This can be done using igv and a fasta file. You will need a fasta file and a gtf file (or a gff3). 


#title           :bash MAIN_bamtotdfdir.sh
#description     :This script take in a directory of bam files and make tdf files. All the bam files must be single end. This script uses several other scripts located in the same directory and submits things to tourque. 
#author          :Mary Ann Allen
#date            :20170214
#version         :0.1    
#usage           :bash runbamtotdf_singleend.sh
#notes           :Edit the in directory at the top of the script to be your directory of bam files.  All the bam files must be single end. This script uses several other scripts located in the same directory and submits jobs to tourque.



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
echo $rootfilename
echo $scriptdirname
echo $outdir
echo $bedgraphdir
echo $bedgraphfortdfdir
echo $tdfdir
echo $bedgraphgenomefile
echo $igvtoolsgenomefile

qsub -v scriptdirname=$scriptdirname,indir=${indir},rootfilename=$rootfilename,outdir=$outdir,bedgraphdir=$bedgraphdir,bedgraphfortdfdir=$bedgraphfortdfdir,tdfdir=$tdfdir,bedgraphgenomefile=$bedgraphgenomefile,igvtoolsgenomefile=$igvtoolsgenomefile -N pipe${rootfilename} ${scriptdirname}bamtotdffile.sh

done


