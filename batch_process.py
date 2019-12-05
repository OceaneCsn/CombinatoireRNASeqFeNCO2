import sys
import os
import argparse
import signal

# handle interruptions
def signal_handler(signal, frame):
  sys.exit(0)

signal.signal(signal.SIGINT, signal_handler)

#----------------------------------- PIPELINE RNA-SEQ -----------------------------------#
# traite des données single end afin de les filtrer, les aligner, et les quantifier.
# fait dans un cadre d'analyse différentielle mutli facteur

#----------------------------------- Arguments -----------------------------------#


parser = argparse.ArgumentParser(description='Batch processing to facilitate RNA-Seq analysis')

parser.add_argument('--i', metavar='input_folder', type=str,
                    help='The folder containing input files', default="./filtered_fqs")
parser.add_argument('--o', metavar='output_folder', type=str,
                    help='The folder containing output files', default="./StarMapping")
parser.add_argument('--nThreads', metavar='NThreads', type=int,
                    help='Number of working threads', default=5)
parser.add_argument('--start', metavar='start', type=int,
                    help='Start', default=1)
parser.add_argument('--stop', metavar='stop', type=int,
                    help='Stop', default=24)
parser.add_argument('--mode', type = str,
                    default="mapping",
                    help='Process to execute: please choose between quality, mapping, quantif, or la_totale')

args = parser.parse_args()

#----------------------------------- Process functions -----------------------------------#

print(args.start, args.stop, args.mode, args.nThreads, args.o, args.i)

samples = range(args.start,args.stop+1)

print("Samples to process : ", args.start, "to", args.stop)
def quality():
    # filter paired end reads and produces html quality report
    for sample in samples:
        s = str(sample)
        print("Filtering ", s)
        #fastp can't use more than 16 worker threads 
        os.system("fastp -i "+args.i+"/R"+s+"_1.fq.gz -I "+args.i+"/R"+s+"_2.fq.gz --detect_adapter_for_pe \
             -h ./fastp_reports/R"+s+".html -o "+args.o+"/R"+s+"_1.fq.gz -O \
             "+args.o+"/R"+s+"_2.fq.gz --thread "+str(args.nThreads))


def mapping():
    # maps reads to TAIR10 ref genome using STAR
    for sample in samples:
        s = str(sample)
        print("Mapping ", s)
        #unlimited threads
        os.system("STAR --readFilesIn "+args.i+"/R"+s+"_1.fq.gz "+args.i+"/R"+s+"_2.fq.gz \
            --readFilesCommand gunzip -c --genomeDir /data/ocassan_data/TAIR10/ --outFileNamePrefix \
            "+args.o+"/R"+s+" --outSAMtype BAM Unsorted --runThreadN "+str(args.nThreads))

        print("STAR --readFilesIn "+args.i+"/R"+s+"_1.fq.gz "+args.i+"/R"+s+"_2.fq.gz \
            --readFilesCommand gunzip -c --genomeDir /data/ocassan_data/TAIR10/ --outFileNamePrefix \
            "+args.o+"/R"+s+" --outSAMtype BAM Unsorted --runThreadN "+str(args.nThreads))

def quantif():
    command = "htseq-count -f bam"
    for sample in samples:
        s = str(sample)
        command += args.i+"/R"+s+"Aligned.out.bam "
    command += "/data/ocassan_data/TAIR10/TAIR10_GFF3_genes.gff --idattr=Parent --stranded no > ./quantif_At.txt"
    print(command)

def quantif_parallel():
    for sample in samples:
        s = str(sample)
        print("tmux new-session -d htseq-count -f bam "+args.i+"/R"+s+"Aligned.out.bam \
            /data/ocassan_data/TAIR10/TAIR10_GFF3_genes.gff --idattr=Parent --stranded no > ./Quantif/quantif_R"+s+".txt")
        os.system("tmux new-session -d htseq-count -f bam "+args.i+"/R"+s+"Aligned.out.bam \
            /data/ocassan_data/TAIR10/TAIR10_GFF3_genes.gff --idattr=Parent --stranded no > ./Quantif/quantif_R"+s+".txt")
        
    #os.system(command)

#----------------------------------- Main -----------------------------------#

if args.mode == "quality":
    quality()
if args.mode == "mapping":
    mapping()
if args.mode == "quantif":
    quantif_parallel()
else:
    print("Invalid mode, please choose between quality, mapping, or quantif")