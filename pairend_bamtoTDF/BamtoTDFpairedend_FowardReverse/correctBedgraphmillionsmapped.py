from __future__ import division
import sys

def calmp(num_of_reads, total_reads):
        mp = int(num_of_reads)/(int(total_reads)/1000000)
        return mp



def mappedfromflagstat(flagstatfile):
	f = open(flagstatfile)
        lines = f.readlines()
        mapped_reads =lines[2]
        mapped_reads = int(mapped_reads.split(" ")[0])
        print flagstatfile, mapped_reads
	return mapped_reads

def mappedfrommillionsmapped(millionsmappedfile):
	f = open(millionsmappedfile)
	line = f.readline()
	line = line.strip("\n")
	mappedreads = int(line)
	return mappedreads


def main(Bedgraphfile, millionsmappedfile):
	total_reads = mappedfrommillionsmapped(millionsmappedfile)
	f = open(Bedgraphfile)
	bedgraphout = Bedgraphfile+".mp.BedGraph"
	wf = open(bedgraphout, "w")
        print "createing", bedgraphout
        line = f.readline()
	while line:
		line = line.strip("\n")
                line = line.split("\t")
                if len(line)<3:
                	try:
                        	line = line[0].split(" ")
                        except:
                        	print line
                chr, start, stop, num_of_reads = line
                frag = calmp(float(num_of_reads), total_reads)
                newline = "\t".join([chr, start, stop, str(frag)])+"\n"
                wf.write(newline)
                line = f.readline()
	wf.close()
	f.close()




if __name__=="__main__":
	if len(sys.argv)>1:
		Bedgraphfile = sys.argv[1]
		millionsmappedfile = sys.argv[2]
		main(Bedgraphfile, millionsmappedfile)
	else:
		print "python correctBedgraphmillionsmapped.py <Bedgraphfile after cat with pos and negative strands> <millionsmappedfile file>"
 
