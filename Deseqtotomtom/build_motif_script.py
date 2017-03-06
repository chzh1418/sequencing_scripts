from __future__ import division
import os
import os.path
import sys
import glob
import pybedtools
from pybedtools import BedTool
import errno

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

def fixbedfile(bedfilelocation, outdir):
	f2 = open(bedfilelocation)
	#find out the number of columns in the bed file
	lenbedfilecols = 0
	for i, linefrombed in enumerate(f2):
		alinefrombed = linefrombed.strip("\n")
		alinefrombed = alinefrombed.split("\t")
		if not linefrombed.startswith("#"):
		        lenbedfile = len(alinefrombed)
		if lenbedfile!=0:
			break
	print "number of columns in bed file", lenbedfile
	f2.close()
	num = 0
	names, lines = [], []
	bedfilelocationsroot = bedfilelocation.split("/")[-1]
	#build a new bed file that has has only uniq names and if didn't have names now has names and removes any lines that are less than 0
	f2 = open(bedfilelocation)
	newfilename = outdir+bedfilelocationsroot+"nameuniq.bed" 
	wf = open(newfilename, "w")
	namenum = {}
        for i, linefrombed in enumerate(f2):
		if linefrombed.startswith("#"):
			wf.write(linefrombed)
		else:	
			alinefrombed = linefrombed.strip("\n")
		        alinefrombed = alinefrombed.split("\t")
			if lenbedfile==3:
				chr, start, stop = alinefrombed
				name = 	chr+":"+start+":"+stop
				if not name in names:
					names.append(name)
					line = [chr, start, stop, name, "0", "."]
					wf.write("\t".join(map(str, line))+"\n")
			elif lenbedfile<=4:
				if lenbedfile==4:
					alinefrombed = alinefrombed+["0", "."]
				elif lenbedfile==5:
					alinefrombed = alinefrombed+["."]
				chr, start, stop, name = alinefrombed[0:4]
				if not name in names:
                                	names.append(name)
					lines.append(linefrombed)
					namenum[name]=1
					wf.write("\t".join(map(str, alinefrombed))+"\n")	
				else:
					nameindex = names.index(name)
					otherline = lines[nameindex]
					if linefrombed==otherline:
						print "repeated line removed", linefrombed
					else:
						num = namenum[name]
						alinefrombed[3] = name+"_"+num
						wf.write("\t".join(map(str, alinefrombed))+"\n")
						names.append(name+"_"+num)
						lines.append(linefrombed)
						num = num +1
						namenum[name] = num
						print "name repeated, but line different"
						print linefrombed
						print otherline
			else:
				print "can't fix this bed file"
	f2.close()
	wf.close()
	return newfilename, lenbedfile






def createDeseqfile(countsfile, conditions, countlanes, outdir, condition1, condition2, cutoff=0.1):
        outfile = countsfile+"."+condition1+condition2
        write_file = outfile+".R"
        print write_file
        wf = open(write_file ,"w")
        R_dump_file = outfile+".Rout"
        graph_file = outfile+".png"
        outfileallinputs = outfile+".res.txt"
        outfilesig = outfile+".resSig.txt"
        outfilesig_orderpval = outfile+".resSig_pvalue.txt"
        wf.write('sink("'+R_dump_file+'")\n')
        wf.write('library( DESeq )\n')
        wf.write('data <- read.delim("'+countsfile+r'", sep="\t", header=FALSE)'+"\n")#need to check that \t comes out like it should. Might write it wrong.
	line = 'countsTable <- subset(data, select=c('
	for i, lane in enumerate(countlanes):
                lanetokeep=int(lane)+1#R starts counting with 1 python with 0
                thecondition = conditions[i]
                if thecondition==condition1:
                        line = line+str(lanetokeep)+","
			print "R will load ", lanetokeep, condition1
                if thecondition==condition2:
                        line = line+str(lanetokeep)+","
			print "R will load ", lanetokeep, condition2
	line = line[0:-1]#removes last comma
	line = line +"))\n"
	wf.write(line)
        wf.write('rownames(countsTable) <- data[[4]]\n')
	line = 'conds <- c('
        for conditionmatch in conditions:
                if conditionmatch==condition1:
                        line = line+'"'+str(conditionmatch)+'",'
			print "R will load ", conditionmatch
                if conditionmatch==condition2:
                        line = line+'"'+str(conditionmatch)+'",'
			print "R will load ", conditionmatch
        line = line[0:-1]#removes last comma
        line = line+")\n"
        wf.write(line)
        wf.write('cds <- newCountDataSet( countsTable, conds )\n')
        wf.write('cds <- estimateSizeFactors( cds )\n')
        wf.write('sizeFactors(cds)\n')
        numcondition1 = conditions.count(condition1)
	numcondition2 = conditions.count(condition2)
	if numcondition1>1 and numcondition2>1:
		wf.write("cds <- estimateDispersions( cds, method='pooled', sharingMode='fit-only' )\n")
	else:
		wf.write("cds <- estimateDispersions( cds, method='blind', sharingMode='fit-only' )\n")
	wf.write('res <- nbinomTest( cds, "'+condition1+'", "'+condition2+'" )\n')
	wf.write('plotDE <- function( res ) plot(res$baseMean,res$log2FoldChange, log="x", pch=20, cex=.1, col = ifelse( res$padj < '+str(cutoff)+', "red", "black" ) )\n')
        wf.write("png('"+graph_file+"')\n")
        wf.write('plotDE( res )\n')
        wf.write('dev.off()\n')
        wf.write('resSig <- res[ res$padj < '+str(cutoff)+', ]\n')
        wf.write('write.table(res, file = "'+outfileallinputs+r'", append = FALSE, sep = "\t")'+"\n")
        wf.write('write.table(resSig, file = "'+outfilesig+r'", append = FALSE, sep = "\t")'+"\n")
        wf.write('write.table(resSig[ order(resSig$pval), ], file = "'+outfilesig_orderpval+r'", append = FALSE, sep = "\t")'+"\n")
        wf.write('sink()\n')
	wf.close()
	return outfile+".R"





def create_sh_script(bamfilescondition1,bamfilescondition2, bedfile, outdir, condition1, condition2, cutoff,queue="long8gb",processers=1):
	print bedfile
	ensure_dir(outdir)
	bedfilelocationsroot = bedfile.split("/")[-1]
	newbedfile = outdir+bedfilelocationsroot+"nameuniq.bed"
	if not os.path.isfile(newbedfile):	
		newbedfile, lenbedfile = fixbedfile(bedfile, outdir)
	wf = open(outdir+"motifpipe.sh", "w")
        wf.write("#!/bin/bash\n")
        wf.write("#PBS -N counting_"+condition1+"_"+condition2+"\n")
        wf.write("#PBS -l nodes=1:ppn="+str(processers)+"\n")
        wf.write("#PBS -q "+queue+"\n")
        wf.write("#PBS -m ae\n")
	wf.write("#PBS -l walltime=96:00:00\n")

	wf.write("#PBS -M Mary.A.Allen@colorado.edu\n")
	wf.write("#PBS -V\n\n\n")
	bamfilecolumns = bamfilescondition1+bamfilescondition2
	bamfilecolumns.sort()
	print bamfilecolumns
	line = "###information used \n"
	wf.write(line)
	oline = [bamfilescondition1,bamfilescondition2, bedfile, newbedfile, outdir, condition1, condition2,"long8gb",1]
	for info in oline:
		print info
		wf.write("###"+ str(info)+"\n")
	line = "### the bam files are in this order in the count table"

	bedfileroot = newbedfile.split("/")[-1]
	countfile = outdir+bedfileroot+".count.bed"
	line = "/opt/bedtools/2.16.2/multiBamCov -bams "
	conditions, countlanes = [],[]
	for i, infile in enumerate(bamfilecolumns):
		column = i+6
		countlanes.append(column)
		if infile in bamfilescondition1:
			conditions.append(condition1)
		if infile in bamfilescondition2:	
			conditions.append(condition2)
		print column, infile	
		baifile = os.path.isfile(infile +".bai") 	
		if baifile!=True:
			print not baifile 
			print "make sure to run"
			print "/opt/samtools/0.1.19/samtools index "+infile
		line = line +infile +" "
	line = line + "-bed "+newbedfile+" >"+countfile+" 2>"+countfile+".err\n\n\n"
	wf.write(line)
	Rfile = createDeseqfile(countfile, conditions, countlanes, outdir, condition1, condition2, cutoff=cutoff)
	print "this is the outdir",outdir
	outfile = countfile+"."+condition1+condition2
	line = "/opt/R/3.1.0/bin/R CMD BATCH --no-save --no-restore "+Rfile+"\n\n\n"
        wf.write(line)
	outfilesig = outfile+".resSig.txt"
	line = "python /projects/dowellLab/groseq/findTFfromGRO/scripts/python_rebuild_the_bed.py "+outfilesig+" "+newbedfile+"\n\n\n"
	wf.write(line)
	outfilessuffixs = [[".bed",".nothit.bed"],[".up.bed",".notup.bed"],[".down.bed",".notdown.bed"]] 
	for suffix,nsuffix in outfilessuffixs:	
		line = "if [ -s "+outfilesig+suffix+" ]\nthen\n"
		wf.write(line)
		line = "/opt/bedtools/2.16.2/fastaFromBed -fi /projects/dowellLab/genomes/human/hg19/hg19ucsc/hg19_all.fa -bed "+outfilesig+suffix+" -fo "+outfilesig+suffix+".fasta\n"
		wf.write(line)
		line = "/opt/bedtools/2.16.2/fastaFromBed -fi /projects/dowellLab/genomes/human/hg19/hg19ucsc/hg19_all.fa -bed "+outfilesig+nsuffix+" -fo "+outfilesig+nsuffix+".fasta\n"
		wf.write(line)
		line = "/opt/meme/4.10.1_4/bin/dreme -p "+outfilesig+suffix+".fasta -n "+outfilesig+nsuffix+".fasta -oc "+outfilesig+suffix+"_dreme_out\n"
		wf.write(line)
		line = "/opt/meme/4.10.1_4/bin/tomtom -no-ssc -oc "+outfilesig+suffix+"_dreme_out/tomtom/ -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 10.0 "+outfilesig+suffix+"_dreme_out/dreme.html /opt/meme/4.10.1_4/db/motif_databases/JASPAR/JASPAR_CORE_2014_vertebrates.meme /opt/meme/4.10.1_4/db/motif_databases/MOUSE/uniprobe_mouse.meme /opt/meme/4.10.1_4/db/motif_databases/MULTI/jolma2013.meme\n"
		wf.write(line)
		line = "fi\n"	
		#/opt/meme/4.10.1_4/bin/tomtom -no-ssc -oc down_dreme_out_control/tomtom/ -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 10.0 down_dreme_out_control/dreme.html /opt/meme/4.10.1_4/db/motif_databases/JASPAR/JASPAR_CORE_2014_vertebrates.meme /opt/meme/4.10.1_4/db/motif_databases/MOUSE/uniprobe_mouse.meme /opt/meme/4.10.1_4/db/motif_databases/MULTI/jolma2013.meme
		wf.write(line)
	wf.close()


