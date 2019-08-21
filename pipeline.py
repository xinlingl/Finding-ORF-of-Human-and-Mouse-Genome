import sys
from os.path import basename
import os
from operator import itemgetter
from Bio import SearchIO
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastpCommandline

bowtie_path = "../bowtie-1.1.2/"
STAR_path = "../STAR/bin/Linux_x86_64/STAR"

# run RiboCode in one step
def riboCode(fasta, fastq, rRNA, gtf, density, counter):
	# skipped the step for trimming adapter sequence for ribo-seq data:
	# Using bowtie :
	print("Removing rRNA contamination using Bowtie 1.1.2 ...... ")
	bowtie1 = bowtie_path+"bowtie-build "+rRNA+" "+"rRNA"
	align_output = raw_input("Please provide the name of alignment output file : ")
	align_output = complete_name(align_output, "align")
	bowtie2 = bowtie_path+"bowtie "+"-p 8 -norc --un un_aligned.fastq rRNA -q "+fastq+" "+align_output
	print(bowtie1)
	print(bowtie2)
	os.system(bowtie1)
	os.system(bowtie2)

	# Using STAR :
	print("Aligning the clean reads to reference genome using STAR ...... ")
	star_index = raw_input("Please provide the name of STAR index output file : ")
	os.system("mkdir "+star_index)
	star1 = STAR_path+" --runThreadN 8 --runMode genomeGenerate --genomeDir " + star_index + " --genomeFastaFiles " + fasta + " --sjdbGTFfile " + gtf
	print(star1)
	unalign = raw_input("Please provide the name for unaligned fastq output file : ")
	unalign = complete_name(unalign, "fastq")
	prefix = raw_input("Please provide the prefix for aligment output files using STAR : ")
	star2 = STAR_path + " --outFilterType BySJout --runThreadN 8 --outFilterMismatchNmax 2 --genomeDir "+star_index+" --readFilesIn "+unalign+" --outFileNamePrefix "+prefix+" --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --alignEndsType EndToEnd"
	print(star2)
	os.system(star1)
	os.system(star2)

	# Using RiboCode:
	print("Identifying translated ORFs using RiboCode ......")
	bam_output = prefix+"Aligned.toTranscriptome.out.bam"
	ribo_out = raw_input("Please provide the name for RiboCode ORF output file : ")
	ribo_cmd = "RiboCode_onestep -g "+gtf+" -f "+fasta+" -r "+ bam_output + " -l no -o "+ribo_out
	print(ribo_cmd)
	os.system(ribo_cmd)
	# Create density plot if user defined:
	#if (density):
	#	print("Plotting the densities of P-sites for predicted ORFs ......")

	#if (counter):
	#	print("Counting the number of RPF reads aligned to ORFs ......")

	print("Finished Successfully")

	return ribo_out

def BLAST(txt_file):
	result=[]
	count=0
	with open(txt_file) as fp:
    	for line in fp:
        	count=count+1
        	line=line.split()
        	if line[1]=="novel":
            	result.append(line)
	result=sorted(result, key=itemgetter(27))
	with open("my_blast.html", "w") as out_handle:
    	for index in range(0,len(result)):
        	result_handle = NCBIWWW.qblast("blastp", "landmark", result[index][28])
        	out_handle.write(result_handle.read())
	return 


def complete_name(name, suffix):
	result = name
	if ("." not in name):
		result = name+"."+suffix
	return result


if __name__ == '__main__':
	if(len(sys.argv) <= 1):
		print("Missing Argument")
	else:	
		# get number of arguements
		# user defined arguements starts from argv[1], i.e., the second arguement
		num_arg = len(sys.argv)
		current = 1
		plot_density = False
		counter = False
		# reading user command line arguments
		# and save all input file names
		while(current < num_arg):
			current_len = len(sys.argv[current])
			if (current_len > 6 and sys.argv[current][current_len-6:current_len] == ".fastq"):
				fastq_input = sys.argv[current]
				current = current + 1
			elif (current_len > 4 and sys.argv[current][current_len-4:current_len] == ".gtf"):
				gtf_input = sys.argv[current]
				current = current + 1
			elif (current_len > 3 and sys.argv[current][current_len-3:current_len] == ".fa"):
				if (sys.argv[current][0:4] == "rRNA"):
					rRNA_fa = sys.argv[current]
					current = current + 1
				else:
					fasta_input = sys.argv[current]
					current = current + 1
			elif (sys.argv[current] == "-d"):
				plot_density = True
			elif (sys.argv[current] == "-c"):
				counter = True
			else:
				print("Invalid Argument:")
				print(sys.argv[current])
				break;




	# print("Input fastq file is:")
	# print(fastq_input)
	# print("Input DNA fasta file is:")
	# print(fasta_input)
	# print("Input rRNA fasta file is:")
	# print(rRNA_fa)
	# print("Input gtf file is:")
	# print(gtf_input)

	output_txt = riboCode(fasta_input,fastq_input, rRNA_fa, gtf_input, plot_density, counter)
	BLAST(output_txt)
