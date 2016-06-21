# diff_abundance
Pipeline for Differential Abundance analyses

Scripts to prepare raw fastq sequences for differential abundance analyses using the Tophat/Cufflinks pipeline as well as the DeSeq pipeline 

1: You need to have a transcriptome fasta file and also a GFF file. Make sure you know the PATH to these files.
2: For use on Skynet, you need to export the PATH to the version of the program you want to use. This can be either installed in your working HOME directory somewhere (e.g. /data/user/programs/ OR it can be installed in the data1 directory.

For this pipeline, you will first ALIGN raw FASTQ files of your treatments to the reference using Tophat. This is done using the 
Tophat_alignments.sh script. 

The script needs the PATH to TOPHAT and your TRANSCRIPTOME and RAW DATA
I recommend making a directory for each project that contains the project as the main directory, followed by sub-directories.
e.g. 

mkdir BUGS
cd BUGS
mkdir 1_Tophat_Alignments

then you can use the OUT as /BUGS/1_Tophat_alignments

to organize everything. 

FOR THE RAW FILES:

The code is written to loop through your RAW files if they are labeled in the following format: NAME.R1.fastq NAME.R2.fastq
IF your files have a lowercase R, edit the code from $i.R1.fastq to $i.r1.fastq and it should work all the same.

To run the script you should be able to type:

./Tophat_alignments.sh


Will add more info on file structure soon. 

RawReads - Create a file with RawReads for each experiment. You can change the naming convention but the script expects files to be named using "X.R1.fastq" and "X.R2.fastq". 
