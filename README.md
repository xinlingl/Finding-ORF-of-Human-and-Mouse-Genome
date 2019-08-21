# Finding-ORF-of-Human-and-Mouse-Genome
03-713 Group C
Yuting Xiao, Xiaoqi Fang, Xinling Li, Zhenyu Yang
User Manual
Our pipeline is designed to incorporate multiple ORF prediction softwares with different
input datasets. 
The final result returned by our pipeline is chosen and considered with high
quality result. Softwares used in our pipeline including Bowtie, STAR, RiboCode, ORF-Finder
and BLAST. The program can be ran in a single python program named “pipeline.py”.

Part 1: Obtain Result by Running RiboCode
Language: Python 2.7

Installation:
1. For RiboCode
a. Install by pypi:
pip install ribocode
b. Install by conda:
Conda install -c bioconda riboconda
c. Install from source:
git clone https://www.github.com/xzt41/RiboCode
cd RiboCode
python setup.py install
d. Install from local:
pip install RiboCode- *.tar.gz
2. For Bowtie:
a. download bowtie 1.1.2 version
3. For STAR:
a. download software:
i. Method 1:
wget https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz
tar -xzf 2.5.3a.tar.gz
cd STAR-2.5.3a
ii. Method 2:
git clone https://github.com/alexdobin/STAR.git
cd STAR/source
b. Build STAR:
i. In Linux:
make STAR
ii. In Mac OS X:
make STARforMacStatic
iii. All platforms:
make STARforMacStatic CXX=/path/to/gcc
4. For BLAST:
sudo apt-get install samtools
sudo apt-get install seqtk
Required Inputs:
Required input files include:
1. Ribosome profiling data in fastq format
2. Whole genome sequencing data in fasta format
3. rRNA data in fasta format, file name should start with “rRNA”
4. GTF of the second file

Outputs:
Percentage of annotated(overlapping annotated CDSs, have the same stop with annotated
CDS) and novel ORFs(in non-coding genes or non-coding transcripts of coding genes), uORFs,
dORFs, overlap_uORFs, and overlap_dORFs. You may asked to type in the file names while
running the program.
Start and end positions of ORFs in each transcript
ORFs with the same stop codon in different transcript isoforms

Running Pipeline:
All required input files should be gathered into a folder, which is placed in our pipeline main
folder, along with installed Bowtie, STAR, etc. Aftering cd into the data folder, users can run the
pipeline using command line:
python ../pipeline.py input1 input2 input3 input4
Note : When running, user may need to provide file names for intimidated files. It will take
about 10 minutes to run Bowtie, 1 hour and 40 minutes to run STAR, and 15 minutes to
run RiboCode.

How to query on BLAST?
BLAST finds regions of similarity between biological sequences. The program compares
nucleotide or protein sequences to sequence databases and calculates the statistical
significance.
1)Run on your terminal,you will get tmp.txt.:
awk '{if ($2 == "novel") {print}}' RiboCode_ORFs_result.txt > tmp.txt
2) Run on your terminal,you will get sequences.txt.:
awk '{print $29}' tmp.txt > sequences.txt
3)Open a BLAST website:
https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_L
OC=blasthome
Copy the contents in sequences.txt to this web site. You will get the results. It usually
takes days to run this step.


Part 2: Improve the performance of Aho–Corasick algorithm by Ribocode-annotated
transcripts
In this part, user will run a revised aho-corasick algorithm.
User can choose to run from the very beginning(including the bowtie and STAR processing
step, which is really time-consuming), or to create the annotated transcripts using RiboCode
prepare-transcript step and then find orfs, or choose the easiest way, to directly use the
annotation file create by RIboCode in part 1. More details will be revealed when running the
program.
(ps: this step needs some reference files from step 1, so running step 1 in advanced is
recommended.
Installation:
1. Ribocode
We assume you have successfully installed RiboCode in the first step.
2. Biopython
Run this command in terminal to install biopython:
conda install -c anaconda biopython
3. STAR:
You have already installed it in part 1! We will still use STAR and some of its output in
data preprocessing step.
Environment:
Language: python 2.7
Input and output
Input: a fasta file of the genome
Output: a gtf file containing the annotation information
How to run the program
1. Data preparation
You should have 3 py files in this folder:
(1). orf_finder.py
(2). pipeline_part2.py
(3). Prepare_transcripts.py
2. The workflow of the program.
Make sure all the data required are in the same folder.
There are three start point you can choose to run the program. For different ways to run
this program, you will need different files in the folder.
(1). If you want to run pipeline_part2.py from the start, you will need exactly same files in
the part 1. However, this process will takes up to 5~6 hours depend on your dataset and
computer performance.
(2). If you want to start after the STAR step, but before the RiboCode data transcript
preparation step, you need two files in the same folder: the GTF file, and the STAR aligned
FASTA file. (this step may take up to an hour)
(3). If you want to start directly, you need only one file: the annotated transcript created
by RiboCode. Feel free to make it stay in the RiboCode_annot directory, the program will fetch it
for you!
3. How to run the pipeline?
In the terminal, type the following sentence to run the part2 of the pipeline.
python pipeline_part2.py
If you get scanning orf: 0% during the process, don’t panic. The program is running and
all you have to do is to wait.
Part III Final Results:
1. ORFs’ category: It should be a pdf file in ‘RiboCode_ORFs_result_ORFs_category.pdf’.
It usually contains a pie chart which indicates the proportion of annotated ORF, uORF,
dORF, Overlapped, novel PCGs, novel non PCGs.
2. General ORF results: To check you get the right the results, the category of the
columns(first row) should start with ‘ORF_ID’，’ORF_type‘, ‘transcript_id’. The end of the
row should be ‘pval_frame0_vs_frame2’, ‘pval_combined’, ‘AAseq’.
3. Collapsed ORFs result: It is a condensed version of the orf results.
‘RiboCode_ORFs_result_collapsed.txt’ should be observed.
4. ORF.counts. If you run the optional command line, you should observe the counts of the
orfs.
5. BLAST: sample result for query a single protein
sequence(MPVARSWVCRKTYVTPRRPFEKSRLDQELKLIGEYGLRNKREVWRVKFTLA
KIRKAARELLTLDEKDPRRLFEGNALLRRLVRIGVLDEGKMKLDYILGLKIEDFLERRLQT
QVFKLGLAKSIHHARVLIRQRHIRYHLGWAPESSSTCPSDGCPH)
Graphic Summary:
Descriptions:
Alignments:
References:
1. http://bowtie-bio.sourceforge.net/manual.shtml#the-bowtie-build-indexer
2. https://github.com/zkstewart/orf-finder-py
3. https://github.com/xryanglab/RiboCode
