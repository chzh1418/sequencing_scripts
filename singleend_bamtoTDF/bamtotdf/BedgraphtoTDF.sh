### Run in desired queue
#PBS -q long8gb

### Use the bash shell
#PBS -S /bin/bash

### Specify the number of nodes and processors for your job
#PBS -l nodes=1:ppn=32

### Switch to the working directory; by default TORQUE launches processes
### from your home directory.  This is a good idea because your -o and -e files 
### will go here
cd $PBS_O_WORKDIR
echo Working directory is $PBS_O_WORKDIR

### Retrieve/use all modules loaded ###
#PBS -V

module load igvtools_2.1.24

echo 

#bamroot=$rootfilename,bgft=$bedgraphfortdfdir,tdfd=$tdfdir
/opt/igvtools/2.1.24/igvtools toTDF ${bgft}${bamroot}.BedGraph.mp.BedGraph ${tdfd}${bamroot}.tdf $igvgenome
