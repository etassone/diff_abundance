 #!/bin/bash
 #Erica Tassone
 #May 12, 2016
 
 #Align files to reference for differential abundance using TopHat
 
 #type path to Tophat here so you can access it from the script
 export PATH=$PATH:/path/to/tophat-2.0.14.Linux_x86_64/
  
 #Define some variables for Tophat 
 TRANSCRIPTOME='/path/to/Transcriptome' #Type path to transcriptome here
 OUT='/path/to/outputs' #Path you'd like the outputs to go
 RAW='/path/to/RawData' #Path to the raw fastq files that you want to use for alignments
 
#Loop through raw data and take the prefix for each file - should print only ONE per treatment
for i in `ls $RAW | awk -F "." '{print $1}'|sort|uniq`;
do
	if [ -f"i" ]; then
	#echo $RAW/$i #to test if you're getting the right file names
	mkdir $OUT/$i #Make the output directory
	#Run Tophat. Change settings as needed. 
	tophat --o $OUT/$i -p 6 -G $TRANSCRIPTOME/transcriptome.new.gff $TRANSCRIPTOME/transcriptome.new.fsa $RAW/$i.R1.fastq $RAW/$i.R2.fastq
	echo Finished alignment of $i
fi
done

