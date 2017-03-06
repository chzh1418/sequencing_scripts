#PBS -q long8gb
#PBS -V
#PBS -l nodes=1:ppn=64
#PBS -S /bin/bash



### Switch to the working directory; by default TORQUE launches processes
### from your home directory.  This is a good idea because your -o and -e files 
### will go here
cd $PBS_O_WORKDIR
echo Working directory is $PBS_O_WORKDIR


module load python_2.7.3
module load samtools_0.1.19
module load bedtools2_2.22.0
module load igvtools_2.1.24



bamroot=$rootfilename
id=$indir
od=$outdir
bgd=$bedgraphdir
bgft=$bedgraphfortdfdir
bgg=$bedgraphgenomefile
bd=${scriptdirname}
tdfd=$tdfdir
igvgenome=$igvtoolsgenomefile



/opt/samtools/0.1.19/samtools sort -m500000000000 ${id}${bamroot}.bam ${od}${bamroot}.sorted
echo finished bam sort
/opt/samtools/0.1.19/samtools flagstat ${id}${bamroot}.bam > ${od}${bamroot}.bam.flagstat 2>${od}${bamroot}.bam.flagstat.err
echo finished counting infile
/opt/samtools/0.1.19/samtools index ${od}${bamroot}.sorted.bam
echo finished indexing sorted bam file
/opt/samtools/0.1.19/samtools flagstat ${od}${bamroot}.sorted.bam > ${od}${bamroot}.sorted.bam.flagstat 2>${od}${bamroot}.sorted.bam.flagstat.err
echo finished counting sorted bam file
diff ${outdir}${rootfilename}.sorted.bam.flagstat ${outdir}${rootfilename}.bam.flagstat
echo finished checking that the sorted bam file has the same number of reads as the unsorted
/opt/bedtools/2.22.0/genomeCoverageBed -bg -ibam ${od}${bamroot}.sorted.bam -g $bgg > ${bgft}${bamroot}.BedGraph.temp
/opt/bedtools/2.22.0/sortBed -i ${bgft}${bamroot}.BedGraph.temp >${bgft}${bamroot}.BedGraph
echo finished sorting the Begraph file with both strands
/opt/python-2.7.3/bin/python ${bd}correctBedgraphmillionsmapped.py ${bgft}${bamroot}.BedGraph ${od}${bamroot}.sorted.bam.flagstat
echo finished correcting for millions mapped
/opt/igvtools/2.1.24/igvtools toTDF ${bgft}${bamroot}.BedGraph.mp.BedGraph ${tdfd}${bamroot}.tdf $igvgenome
echo finished making tdf










