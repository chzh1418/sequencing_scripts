import pandas as pd
import build_motif_script
from build_motif_script import *
import combinebedfiles
from combinebedfiles import *



def useconditionstorundeseqtobam(paper, conditions1, conditions2, outdir, bamdir, deseq_adjpval_cuttoff=0.1):
	print "not done"
	#build a bedfile that is a combination of tfit calls
	#run deseqtobam
	
	create_sh_script(bamfilescondition1,bamfilescondition2, bedfile, outdir, condition1, condition2, deseq_adjpval_cuttoff, queue="long8gb",processers=1)
	


def testshscript1():
	condition1 = "DMSO"
	condition2 = "Nutlin"
	bedfile = p53bedfiles4()
#	bedfile = "/scratch/Users/allenma/scripts/testscript/p53tfitconcatenate.bed"
	outdir = "/scratch/Users/allenma/scripts/testscript8/" 
	bamdir = "/scratch/Shares/dowell/from_jon/Allen2014/bowtie2/sortedbam/"
	bamfilescondition1 = ["SRR1105736", "SRR1105737"]
	bamfilescondition2 = ["SRR1105738", "SRR1105740"]
	bamfilescondition1 = [bamdir+x+".fastqbowtie2.sorted.bam" for x in bamfilescondition1]
	bamfilescondition2 = [bamdir+x+".fastqbowtie2.sorted.bam" for x in bamfilescondition2]
	cutoff = 0.2
	create_sh_script(bamfilescondition1,bamfilescondition2, bedfile, outdir, condition1, condition2, cutoff, queue="long8gb",processers=1)

def p53bedfiles():
	bedfiledir="/scratch/Shares/dowell/md_score_paper/tfit_bed_files/human/081/"
	bedfiles = ["SRR1105736", "SRR1105737", "SRR1105738", "SRR1105739"]
	bedfiles = [bedfiledir+x+"-1_bidir_predictions.bed" for x in bedfiles]	
	prefix ="p53tfit"
	outdir = "/scratch/Users/allenma/scripts/testscript4/" 
	bedfilename = concatenatebedfiles(bedfiles, prefix, outdir, pm=True)
	return bedfilename


def p53bedfiles2():
	bedfiledir="/scratch/Shares/dowell/md_score_paper/tfit_bed_files/human/081/"
        bedfiles = ["SRR1105736", "SRR1105737", "SRR1105738", "SRR1105739"]
        bedfiles = [bedfiledir+x+"-1_bidir_predictions.bed" for x in bedfiles]
        prefix ="p53tfit_atleast2"
        outdir = "/scratch/Users/allenma/scripts/testscript7/"
	bedfilename = concatenatebedfiles_atleast2(bedfiles, prefix, outdir, pm=True)
        return bedfilename


def p53bedfiles3():
	bedfiledir="/scratch/Shares/dowell/md_score_paper/tfit_bed_files/human/081/"
	bedfiles = ["SRR1105736", "SRR1105737", "SRR1105738", "SRR1105739"]
	bedfiles = [bedfiledir+x+"-1_bidir_predictions.bed" for x in bedfiles]	
	prefix ="p53tfit_widen500"
	outdir = "/scratch/Users/allenma/scripts/testscript7/" 
	bedfilename = concatenatebedfiles_widen(bedfiles, prefix, outdir, pm=True)
	return bedfilename


def p53bedfiles4():
	bedfiledir="/scratch/Shares/dowell/md_score_paper/tfit_bed_files/human/081/"
        bedfiles = ["SRR1105736", "SRR1105737", "SRR1105738", "SRR1105739"]
        bedfiles = [bedfiledir+x+"-1_bidir_predictions.bed" for x in bedfiles]
        prefix ="p53tfit_atleast2"
        outdir = "/scratch/Users/allenma/scripts/testscript8/"
	bedfilename = concatenatebedfiles_atleast2_widen(bedfiles, prefix, outdir, pm=True)
        return bedfilename



def Hbedfiles():
	bedfiledir="/scratch/Shares/dowell/md_score_paper/tfit_bed_files/human/081/"
	bedfiles = ["SRR1105736", "SRR1105737", "SRR1105738", "SRR1105739"]
	bedfiles = [bedfiledir+x+"-1_bidir_predictions.bed" for x in bedfiles]	
	prefix ="p53tfit"
	outdir = "/scratch/Users/allenma/scripts/testscript2/" 
	bedfilename = concatenatebedfiles(bedfiles, prefix, outdir)
	return bedfilename

def readconditions(whichorganism="human"):
	con="/scratch/Users/allenma/scripts/Deseqtotomtom/"+"conditions_short_20161103_tentative.txt_20161107-165140.csv"
	f = open(con)	
	title = f.readline
	keywords = []
	for line in f:
		line = line.strip("\n")
		line = line.split(",")
		SRAnumber,organism,tissue,general_celltype,specific_celltype,treatment_code,treated_or_like_treated,repnumber,keyword,exptype,mapped_reads,total_reads,percent_mapped = line
		if organism==whichorganism:
			if keyword not in keywords:
				keywords.append(keyword)
	print keywords

def printinfoonepaper(whichkeyword):
        con="/scratch/Users/allenma/scripts/Deseqtotomtom/"+"conditions_short_20161103_tentative.txt_20161107-165140.csv"
        f = open(con)
        title = f.readline
        keywords = []
	SRRs, tissues, general_celltypes, specific_celltypes, treatment_codes=[],[],[],[],[]
        for line in f:
                line = line.strip("\n")
                line = line.split(",")
                SRAnumber,organism,tissue,general_celltype,specific_celltype,treatment_code,treated_or_like_treated,repnumber,keyword,exptype,mapped_reads,total_reads,percent_mapped = line
                if keyword==whichkeyword:
			print line
			SRRs.append(SRAnumber)
			tissues.append(tissue)
			general_celltypes.append(general_celltype)
			specific_celltypes.append(specific_celltype)
			treatment_codes.append(treatment_code)
	return SRRs, tissues, general_celltypes, specific_celltypes, treatment_codes	

#printinfoonepaper("Heinaniemi2015")
#testshscript1()
#p53bedfiles2()
#readconditions()

def testshscript3():
	condition1 = "DMSO"
	condition2 = "Nutlin"
	bedfile = p53bedfiles_joeysidea()
#	bedfile = "/scratch/Users/allenma/scripts/testscript/p53tfitconcatenate.bed"
	outdir = "/scratch/Users/allenma/scripts/testscript8/" 
	bamdir = "/scratch/Shares/dowell/from_jon/Allen2014/bowtie2/sortedbam/"
	bamfilescondition1 = ["SRR1105736", "SRR1105737"]
	bamfilescondition2 = ["SRR1105738", "SRR1105740"]
	bamfilescondition1 = [bamdir+x+".fastqbowtie2.sorted.bam" for x in bamfilescondition1]
	bamfilescondition2 = [bamdir+x+".fastqbowtie2.sorted.bam" for x in bamfilescondition2]
	cutoff = 0.2
	create_sh_script(bamfilescondition1,bamfilescondition2, bedfile, outdir, condition1, condition2, cutoff, queue="long8gb",processers=1)


