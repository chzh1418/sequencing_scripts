### Run in desired queue
#PBS -q long2gb

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

module load python_2.7.3
module load samtools_0.1.19
module load bedtools2_2.22.0
module load igvtools_2.1.24


echo starting sortedbam submit
echo $rootfilename
echo $indir
echo $outdir
#qsub -v bamroot=$rootfilename,id=$indir,od=$outdir -N sortedbam${rootfilename} ${scriptdirname}bamtosortedbamandbai.sh
#wait until sortedbam exists
echo complete sortedbam submit
until [ -f ${outdir}${rootfilename}.sorted.bam.flagstat ]
do
     sleep 5
done
echo "File found"
diff ${outdir}${rootfilename}.sorted.bam.flagstat ${outdir}${rootfilename}.bam.flagstat
echo starting bedgraph submit
qsub -v bamroot=$rootfilename,od=$outdir,bgd=$bedgraphdir,bgft=$bedgraphfortdfdir,bgg=$bedgraphgenomefile -N BedGraph${rootfilename} ${scriptdirname}bamtoBedgraph.sh
echo completed bedgraph submit
#wait until this is done
until [ -f ${bedgraphfortdfdir}${rootfilename}.BedGraph.done ]
do     sleep 5
done
echo "File found"

echo starting read depth correction submit
qsub -v bd=${scriptdirname},bamroot=$rootfilename,od=$outdir,bgft=$bedgraphfortdfdir -N readcor${rootfilename} ${scriptdirname}Bedgraph_readcorrection.sh
echo completed read depth correction submit
#wait until this is done
until [ -f ${bedgraphfortdfdir}${rootfilename}.BedGraph.mp.BedGraph.done ]
do     sleep 5
done
echo "File found"
echo starting tdf conversion submit
qsub -v bamroot=$rootfilename,bgft=$bedgraphfortdfdir,tdfd=$tdfdir,igvgenome=$igvtoolsgenomefile -N tdf${rootfilename} ${scriptdirname}BedgraphtoTDF.sh
echo completed tdf conversion submit


