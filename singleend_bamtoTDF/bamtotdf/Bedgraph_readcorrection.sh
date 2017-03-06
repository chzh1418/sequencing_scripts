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

#bamroot=$rootfilename,od=$outdir,bgft=$bedgraphfortdfdir
module load python_2.7.3

/opt/python-2.7.3/bin/python ${bd}correctBedgraphmillionsmapped.py ${bgft}${bamroot}.BedGraph ${od}${bamroot}.sorted.bam.flagstat 
touch ${bgft}${bamroot}.BedGraph.mp.BedGraph.done


