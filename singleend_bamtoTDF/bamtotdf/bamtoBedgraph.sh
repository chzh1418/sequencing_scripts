
### Run in desired queue
#PBS -q long8gb

### Use the bash shell
#PBS -S /bin/bash

### Specify the number of nodes and processors for your job
#PBS -l nodes=1:ppn=1

### Switch to the working directory; by default TORQUE launches processes
### from your home directory.  This is a good idea because your -o and -e files 
### will go here
cd $PBS_O_WORKDIR
echo Working directory is $PBS_O_WORKDIR

### Retrieve/use all modules loaded ###
#PBS -V

echo $genome
#bamroot=$rootfilename,od=$outdir,bgd=$bedgraphdir,bgft=$bedgraphfortdfdir,bgg=$bedgraphgenomefile 
module load bedtools2_2.22.0

/opt/bedtools/2.22.0/genomeCoverageBed -bg -strand + -ibam ${od}${bamroot}.sorted.bam -g $bgg > ${bgd}${bamroot}.pos.BedGraph
/opt/bedtools/2.22.0/genomeCoverageBed -bg -strand - -ibam ${od}${bamroot}.sorted.bam -g $bgg | awk -F '\t' -v OFS='\t' '{ $4 = - $4 ; print $0 }'> ${bgd}${bamroot}.neg.BedGraph
cat ${bgd}${bamroot}.pos.BedGraph ${bgd}${bamroot}.neg.BedGraph > ${bgft}${bamroot}.BedGraph.temp
/opt/bedtools/2.22.0/sortBed -i ${bgft}${bamroot}.BedGraph.temp >${bgft}${bamroot}.BedGraph
touch ${bgft}${bamroot}.BedGraph.done







