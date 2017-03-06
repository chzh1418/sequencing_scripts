# sequencing_scripts

These are some useful scripts for sequncing.

They are all mostly wrappers around other peoples code (except for a couple of easy python scripts). 


Below are two scripts to make visualization files. Bam files are large. So I make tdf files and then use igv to view them.

pairend_bamtoTDF
singleend_bamtoTDF

The following code takes a list of bam files and a list of bed files and does differential expression using Deseq. THe results of the Deseq up and down regions are run though dreme to look for enriched motifs. Then tomtom predicts associated TFs.

Deseqtotomtom
 


