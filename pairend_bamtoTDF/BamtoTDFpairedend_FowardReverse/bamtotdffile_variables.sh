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

igvtoolsdir=
samtoolsdir=/opt/samtools/0.1.19/
bedtoolsdir=/opt/bedtools/2.22.0/
pythondir=/opt/python-2.7.3/bin/
bamroot=$rootfilename
id=$indir
od=$outdir
bgd=$bedgraphdir
bgft=$bedgraphfortdfdir
bgg=$bedgraphgenomefile
bd=${scriptdirname}
tdfd=$tdfdir
igvgenome=$igvtoolsgenomefile


#sort the bam file
${samtoolsdir}samtools sort -m500000000000 ${id}${bamroot}.bam ${od}${bamroot}.sorted
echo finished bam sort
#
${samtoolsdir}samtools flagstat ${id}${bamroot}.bam > ${od}${bamroot}.bam.flagstat 2>${od}${bamroot}.bam.flagstat.err
echo finished counting infile
${samtoolsdir}samtools index ${od}${bamroot}.sorted.bam
echo finished indexing sorted bam file
${samtoolsdir}samtools flagstat ${od}${bamroot}.sorted.bam > ${od}${bamroot}.sorted.bam.flagstat 2>${od}${bamroot}.sorted.bam.flagstat.err
#count millions mapped. Flagstat reports to many reads as mapping when you use it on paired end data. 
${samtoolsdir}samtools view -cF 0x100 ${od}${bamroot}.sorted.bam >${od}${bamroot}.sorted.bam.millionsmapped
echo finished counting sorted bam file
#use both flaststats to make sure the samtools conversion didn't lose reads
diff ${outdir}${rootfilename}.sorted.bam.flagstat ${outdir}${rootfilename}.bam.flagstat
echo Does the sorted bam have the same rnumber of reads as the unsored

#http://davetang.org/wiki/tiki-index.php?page=SAMTools#Extracting_only_the_first_read_from_paired_end_BAM_files
${samtoolsdir}samtools view -h -b -f 0x0040 ${od}${bamroot}.sorted.bam > ${bedgraphfortdfdir}${bamroot}.pairfirst.bam
# 0x0040 is hexadecimal for 64 (i.e. 16 * 4), which is binary for 1000000, corresponding to the read in the first read pair.


#need to know the flag for the second strand
# Jess gave me this https://broadinstitute.github.io/picard/explain-flags.html
#128 means second in pair
#128 in hexadecimal is 0x0080
${samtoolsdir}samtools view -h -b -f 0x0080 ${od}${bamroot}.sorted.bam > ${bedgraphfortdfdir}${bamroot}.pairsecond.bam

#That should get me the two parts of pair separate
#Then I need to run genomecoverage and swap the strand info on the second pair


#Then I need to run genomecoverage on each of them
${bedtoolsdir}genomeCoverageBed -bg -split -strand - -ibam ${bedgraphfortdfdir}${bamroot}.pairfirst.bam -g $genome >$outdir/genomecoveragebedpair/${entry}.pairfirst.pos.bed

${bedtoolsdir}genomeCoverageBed -bg  -split -strand + -ibam ${bedgraphfortdfdir}${bamroot}.pairfirst.bam -g $genome | awk -F '\t' -v OFS='\t' '{ $4 = - $4 ; print $0 }'> $outdir/genomecoveragebedpair/${entry}.pairfirst.neg.bed

${bedtoolsdir}genomeCoverageBed -bg -split -strand + -ibam  ${bedgraphfortdfdir}${bamroot}.pairsecond.bam -g $genome > $outdir/genomecoveragebedpair/${entry}.pairsecond.pos.bed

${bedtoolsdir}genomeCoverageBed -bg -split -strand - -ibam ${bedgraphfortdfdir}${bamroot}.pairsecond.bam -g $genome | awk -F '\t' -v OFS='\t' '{ $4 = - $4 ; print $0 }'> $outdir/genomecoveragebedpair/${entry}.pairsecond.neg.bed

#first I need to sort the Bedgraphs
${bedtoolsdir}sortBed -i ${bedgraphfortdfdir}${bamroot}.pairfirst.pos.bed > ${bedgraphfortdfdir}${bamroot}.pairfirst.pos.BedGraph.sort

${bedtoolsdir}sortBed -i ${bedgraphfortdfdir}${bamroot}.pairfirst.neg.bed > ${bedgraphfortdfdir}${bamroot}.pairfirst.neg.BedGraph.sort


${bedtoolsdir}sortBed -i ${bedgraphfortdfdir}${bamroot}.pairsecond.pos.bed > ${bedgraphfortdfdir}${bamroot}.pairsecond.pos.BedGraph.sort

${bedtoolsdir}sortBed -i ${bedgraphfortdfdir}${bamroot}.pairsecond.neg.bed > ${bedgraphfortdfdir}${bamroot}.pairsecond.neg.BedGraph.sort

#Then I need to add two Bedgraphs 

#this should put the values in columns 4 and 5


${bedtoolsdir}unionBedGraphs -i ${bedgraphfortdfdir}${bamroot}.pairfirst.pos.BedGraph.sort ${bedgraphfortdfdir}${bamroot}.pairsecond.pos.BedGraph.sort >${bedgraphfortdfdir}${bamroot}.pos.Bedgraphcol

${bedtoolsdir}unionBedGraphs -i ${bedgraphfortdfdir}${bamroot}.pairfirst.neg.BedGraph.sort ${bedgraphfortdfdir}${bamroot}.pairsecond.neg.BedGraph.sort >${bedgraphfortdfdir}${bamroot}.neg.Bedgraphcol

#then I need to sum cols 4 and 5

awk -F '\t' '{OFS="\t"; print $1,$2,$3,$4+$5;}' ${bedgraphfortdfdir}${bamroot}.pos.Bedgraphcol >${bedgraphfortdfdir}${bamroot}.pos.Bedgraph

awk -F '\t' '{OFS="\t"; print $1,$2,$3,$4+$5;}' ${bedgraphfortdfdir}${bamroot}.neg.Bedgraphcol >${bedgraphfortdfdir}${bamroot}.neg.Bedgraph

#then I need to cat the two Bedgraphs

cat $outdir/genomecoveragebedpair/${entry}.pos.Bedgraph ${bedgraphfortdfdir}${bamroot}.neg.Bedgraph >${bedgraphfortdfdir}${bamroot}.bed

#then I need to sort the final Bedgraph so it can be divided by millions mapped and converted into tdf
${bedtoolsdir}sortBed -i ${bedgraphfortdfdir}${bamroot}.bed >${bedgraphfortdfdir}${bamroot}.Bedgraph



#now correct for millionsmapped
${pythondir}python ${bd}correctBedgraphmillionsmapped.py ${bgft}${bamroot}.BedGraph ${od}${bamroot}.sorted.bam.millionsmapped
echo finished correcting for millions mapped
${igvtoolsdir}igvtools toTDF ${bgft}${bamroot}.BedGraph.mp.BedGraph ${tdfd}${bamroot}.tdf $igvgenome
echo finished making tdf










