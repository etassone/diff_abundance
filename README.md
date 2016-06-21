# diff_abundance
Pipeline for Differential Abundance analyses

For DeSEQ
run 
Tophat_alignments.sh

samtools_sort.sh

HTSeq.scripts.count (Python program, installed with HTSeq)

##NOTE: The GFF file requirements for HTSeq are very specific. You may need to edit your GFF File to be compatable with HTSeq counts. 
For our research, I have made files available. :) 

Run DeSEQ in R - use the .R program in RStudio to get your results. 

For 2 conditions: DESeq2_twoConditions.R

for 3 conditions: Deseq2_multiple.R
