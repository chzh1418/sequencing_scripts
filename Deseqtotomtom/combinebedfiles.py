import build_motif_script
from build_motif_script import *

def concatenatebedfiles(bedfiles, prefix, outdir, pm=True):
	ensure_dir(outdir)
	wf = open(outdir+prefix+"concatenate.info", "w")
	line = "#concatenate and sort the following bedfiles"
	for bedfile in bedfiles:
		line = line +" "+bedfile
	wf.write(line)
	wf.close()
	x = pybedtools.BedTool(bedfiles[0])
	c = x.cat(*bedfiles[1:], postmerge=pm)
	c = c.sort()
	c.saveas(outdir+prefix+"concatenate.bed")
	return outdir+prefix+"concatenate.bed"



def concatenatebedfiles_regionsinatleast2(bedfiles, prefix, outdir, pm=True):
	ensure_dir(outdir)
	wf = open(outdir+prefix+"concatenate.info", "w")
	line = "#concatenate and sort the following bedfiles"
	for bedfile in bedfiles:
		line = line +" "+bedfile+" regions must be in at least 2 of the bedfiles"
	wf.write(line)
	wf.close()
	windowsize = 100
#	x = pybedtools.BedTool(bedfiles[0])
#	c = x.cat(*bedfiles[1:], postmerge=pm)
	a = pybedtools.BedTool(bedfiles[0])
	b = pybedtools.BedTool(bedfiles[1])
	both = a.window(b, w=windowsize)
	both2 = b.window(a, w=windowsize)
	print "in a near b", len(both)
	print "in b near a", len(both2)
	both = both.cat(both2, postmerge=True)
	print "len both", len(both)
	for a in bedfiles:
		for b in bedfiles:
			if a!=b:
				a = pybedtools.BedTool(a)
				b = pybedtools.BedTool(b)
				d = a.window(b, w=windowsize)
				e = b.window(a, w=windowsize)
				both = both.cat(d, postmerge=True)
				print "len both", len(both)
				both = both.cat(e, postmerge=True)
				print "len both", len(both)
	c = both.sort()
	c.saveas(outdir+prefix+"concatenate.bed")
	return outdir+prefix+"concatenate.bed"

def concatenatebedfiles_widen(bedfiles, prefix, outdir, pm=True, widen=500):
	ensure_dir(outdir)
	wf = open(outdir+prefix+"concatenate.info", "w")
	line = "#concatenate and sort the following bedfiles"
	for bedfile in bedfiles:
		line = line +" "+bedfile
	wf.write(line)
	wf.close()
	x = pybedtools.BedTool(bedfiles[0])
	c = x.cat(*bedfiles[1:], postmerge=pm)
	c = c.sort()
	c = c.slop(b=widen, genome="hg19")
	c.saveas(outdir+prefix+"concatenate.bed")
	return outdir+prefix+"concatenate.bed"



def concatenatebedfiles_regionsinatleast2_widen(bedfiles, prefix, outdir, pm=True):
	ensure_dir(outdir)
	wf = open(outdir+prefix+"concatenate.info", "w")
	line = "#concatenate and sort the following bedfiles"
	for bedfile in bedfiles:
		line = line +" "+bedfile+" regions must be in at least 2 of the bedfiles"
	wf.write(line)
	wf.close()
	windowsize = 100
#	x = pybedtools.BedTool(bedfiles[0])
#	c = x.cat(*bedfiles[1:], postmerge=pm)
	a = pybedtools.BedTool(bedfiles[0])
	b = pybedtools.BedTool(bedfiles[1])
	both = a.window(b, w=windowsize)
	both2 = b.window(a, w=windowsize)
	print "in a near b", len(both)
	print "in b near a", len(both2)
	both = both.cat(both2, postmerge=True)
	print "len both", len(both)
	for a in bedfiles:
		for b in bedfiles:
			if a!=b:
				a = pybedtools.BedTool(a)
				b = pybedtools.BedTool(b)
				d = a.window(b, w=windowsize)
				e = b.window(a, w=windowsize)
				both = both.cat(d, postmerge=True)
				print "len both", len(both)
				both = both.cat(e, postmerge=True)
				print "len both", len(both)
	c = both.sort()
	c= c.slop(b=500, genome="hg19")
	c.saveas(outdir+prefix+"concatenate.bed")
	return outdir+prefix+"concatenate.bed"





