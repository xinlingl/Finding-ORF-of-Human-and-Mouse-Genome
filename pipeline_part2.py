#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 29 23:31:44 2018

"""

import os
import shutil


    
'''
For time-saving, if you have alrealy run the part one of the pipeline, 
you can directly get orfs by typing in terminal:
   python orf_finder.py transcripts_sequence.fa <youroutput.fa>
   
Or if you want to use different gtf files and fasta files, you can just run 
this program and follow the instructions.

Or, if you are not confident in the quality of gene file, you could uncomment 
the data-preprocessing portion, and find the orfs in a de novo manner.
'''



def prepare_transcripts(gtf,fasta):
    print("prepare transcripts......")
    transcript_cmd = "prepare_transcripts -g "+gtf+" -f "+fasta+" -o "+"RiboCode_annot"
    print(transcript_cmd)
    os.system(transcript_cmd)
    print("Transcripts Preparation Finished")

def find_orf(transcript, orfout):
    print("let's start finding orfs......")
    orf_cmd = "python orf_finder.py "+transcript+" " +orfout
    print(orf_cmd)
    os.system(orf_cmd)
    print("YOU HAVE FINISHED ALL THE CONTENT OF THIS PROGRAM!")


bowtie_path = "bowtie-1.1.2/"
STAR_path = "/usr1/home/sgona/STAR/bin/Linux_x86_64/STAR"



#this part is adopted from pipeline.py: riboCode
def Bowtie_STAR(fasta,fastq,rRNA,gtf):
    print("removing rRNA contamination using Bowtie......")
    bowtie1 = bowtie_path+"bowtie-build "+rRNA+" "+"rRNA"
    align_output = raw_input("please provide the name of alignment output file, the suffix is <.align> and must be included: ")
    bowtie2 = bowtie_path+"bowtie "+"-p 8 -norc --un_aligned.fastq rRNA -q "+fastq+" "+align_output
    print(bowtie1)
    print(bowtie2)
    os.system(bowtie1)
    os.system(bowtie2)
    
    print("Aligning the clean reads to reference genome using STAR......")
    star_index = raw_input("please provide the name of STAR index output file:")
    star1 = STAR_path +" --runThreadN 8 --runMode genomeGenerate --genomeDir "+star_index + " --genomeFastaFiles "+fasta+" --sjdbGTFfile "+gtf
    print(star1)
    unalign = raw_input("Please provide the name for unaligned fastq output file, suffix must be included: ")
    prefix = raw_input("Please provide the prefix for aligment output files using STAR : ")
    star2 = STAR_path + " --outFilterType BySJout --runThreadN 8 --outFilterMismatchNmax 2 --genomeDir "+star_index+" --readFilesIn "+unalign+" --outFileNamePrefix "+prefix+" --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outFilterMultimapNmax 1 --outFilterMatchNmin 16 --alignEndsType EndToEnd"
    print(star2)
    os.system(star1)
    os.system(star2)
    
    


if __name__ == '__main__':
    print('''
            ### IMPORTANT NOTES ###
        For saving your time, if you have alrealy run the part one of the pipeline, you can directly get orfs by waiting for the second choice and choose yes.
   
        Or if you want to use different gtf files and fasta files, you can choose yes at the first choice and follow the instructions.

        Or, if you are not confident in the quality of gene file, you can use bowtie and STAR to preprocess the data. However, it demands humorous space and is time-consuming.
        
        If you don't want to excute any part of this program, just choose no at the first choice. 
         ''')
    print("\n")
    work = raw_input("Continue? (y/n)")
    if work == "n" or work =="N":
        exit
    else:
        print("Now would you like to find the orfs from annotated GTF and fasta file?(no STAR or Bowtie needed)")
        go = raw_input("Continue? (y/n)")
        if go == 'y' or go == 'Y':
            gtf_name = raw_input("please provide the full name of GTF file (suffix must be included): ")
            fasta_name = raw_input("please provide the full name of gene fasta file(suffix must be included): ")
            transcripts = "transcripts_sequence.fa"
        
            prepare_transcripts(gtf_name, fasta_name)
            outorfname = raw_input("please provide the name of file you want to write orfs (suffix(fa or fasta) must be included)")
            current_workdir = os.getcwd()
            new_workdir =current_workdir+"/RiboCode_annot/"
            filename = current_workdir+"/RiboCode_annot/transcripts_sequence.fa"
            shutil.move(filename, current_workdir)
            os.chdir(current_workdir)
            find_orf(transcripts,outorfname)
        else:
            print("Or you'd like to find the orf without spending too much time? Let's find it directly!")
            find = raw_input("Continue? (y/n)")
            if find == "y" or find == "Y":
                current_workdir = os.getcwd()
                new_workdir =current_workdir+"/RiboCode_annot/"
                filename = current_workdir+"/RiboCode_annot/transcripts_sequence.fa"
                shutil.move(filename, current_workdir)
                #copy_files(new_workdir, current_workdir)
                os.chdir(current_workdir)
                transcripts = "transcripts_sequence.fa"
                outorfname = raw_input("please provide the name of file you want to write orfs (suffix(fa or fasta) must be included)")
                find_orf(transcripts,outorfname)
            else:
                print("I can't believe you choose the hardest way... but if you want to start from the beginning, it's okay.")
                decision = raw_input("Continue from bowtie? it may take 4 hours or more. Make sure you have backup for your data. (y/n)")
                if decision == "n" or decision == "N":
                    print("glad you save your time. See you then!")
                    exit
                else:
                    fasta = raw_input("please provide the full name of gene fasta file(suffix must be included): ")
                    gtf = raw_input("please provide the full name of GTF file (suffix must be included): ")
                    rRNA = raw_input("please provide the full name of rRNA file (suffix must be included): ")
                    fastq = raw_input("please provide the full name of fastq file (suffix must be included): ")
                    Bowtie_STAR(fasta, fastq, rRNA, gtf)
                    fasta_name = raw_input("please provide the full name of gene fasta file created by STAR(suffix must be included): ")
                    transcripts = "transcripts_sequence.fa"
                    prepare_transcripts(gtf, fasta_name)
                    outorfname = raw_input("please provide the name of file you want to write orfs (suffix(fa or fasta) must be included)")
                    current_workdir = os.getcwd()
                    new_workdir =current_workdir+"/RiboCode_annot/"
                    filename = current_workdir+"/RiboCode_annot/transcripts_sequence.fa"
                    shutil.move(filename, current_workdir)
                    os.chdir(current_workdir)
                    find_orf(transcripts,outorfname)
                    
                    
                    
                    
                    
                    
                    
                
                
    
