#!/bin/bash
#Created by Erica Tassone
#May 25, 2015
#Script to sort aligned bam files
#Type PATH to Samtools here if not accessible directly on server
export PATH=$PATH:/path/to/samtools-1.2/

#Type Directory path of where you want your output to go to
DIR='/path/to/output'

#change to directory where your individual alignment folders are located. It's best to have a subfolder for each individual alignment.
cd /path/to/alignments
#Output of tophat is accepted_hits.bam, this file will be in each directory if you aligned with TopHat
for i in */accepted_hits.bam;
do 
	if [ -f"i" ]; then 
		 
		samtools sort -n -m 4G -@ 4 -T tmp -O bam -o $DIR/${i%/*}.sorted.bam $i
		echo Finished sorting ${i%/*}	
fi 
done
 
