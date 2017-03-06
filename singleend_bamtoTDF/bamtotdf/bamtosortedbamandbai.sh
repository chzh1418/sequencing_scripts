#PBS -q long8gb
#PBS -V
#PBS -l nodes=1:ppn=64
#PBS -S /bin/bash

#bamroot=$rootfilename,id=$indir,od=$outdir
cd $PBS_O_WORKDIR
echo Working directory is $PBS_O_WORKDIR

module load samtools_0.1.19

/opt/samtools/0.1.19/samtools sort -m500000000000 ${id}${bamroot}.bam ${od}${bamroot}.sorted
echo finishedsort
/opt/samtools/0.1.19/samtools flagstat ${id}${bamroot}.bam > ${od}${bamroot}.bam.flagstat 2>${od}${bamroot}.bam.flagstat.err
echo finished counting infile
/opt/samtools/0.1.19/samtools index ${od}${bamroot}.sorted.bam
echo finished indexing sorted bam file
/opt/samtools/0.1.19/samtools flagstat ${od}${bamroot}.sorted.bam > ${od}${bamroot}.sorted.bam.flagstat 2>${od}${bamroot}.sorted.bam.flagstat.err
echo finished counting sorted bam file

