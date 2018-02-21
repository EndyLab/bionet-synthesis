from datetime import datetime
start = datetime.now()
print("Starting run at: ",start)

import numpy as np
import pandas as pd

import json
import os
import glob
import re
import math

import shutil

from Bio import pairwise2
from Bio import Seq
from Bio import SeqIO

BASE_PATH = "/Users/conarymeyer/Desktop/GitHub/bionet-synthesis"
PIPELINE_PATH = BASE_PATH + "/pipeline"
BUILDS_PATH = BASE_PATH + "/builds"
DATA_PATH = BASE_PATH + "/data"

BACKBONE_PATH = BASE_PATH + "/sequencing_files/popen_v1-1_backbone.fasta"
DICTIONARY_PATH = PIPELINE_PATH + "/testing/data_testing/10K_CDS.csv"

import_time = datetime.now()
print("Time to import: ",import_time - start)

build_num = "build"+str(input("Which build: ")).zfill(3)
print("build_num",build_num)

SEQFILE_PATH = "{}/{}/{}_seq_files".format(BUILDS_PATH,build_num,build_num)

forward_primer = "M13-Forward---20-"
reverse_primer = "M13-Reverse"

# Create a dictionary to link the gene name to the corresponding id number
data = pd.read_csv(DICTIONARY_PATH)
dictionary = dict(zip(data['gene_name'], data['idnum']))

count = 0

align_data = []
nan = []
small = []
check = []

backbone_seq = Seq.Seq("ATGCGGTCTTCCGCATCGCCTGGCACGACAGGTTTCCCGACTGGAAAGCGGGCAGTGAGCGCAACGCAATTAATGTGAGTTAGCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATGTTGTGTGGAATTGTGAGCGGATAACAATTTCACACATACTAGAGAAAGAGGAGAAATACTAGATGGCTTCCTCCGAAGATGTTATCAAAGAGTTCATGCGTTTCAAAGTTCGTATGGAAGGTTCCGTTAACGGTCACGAGTTCGAAATCGAAGGTGAAGGTGAAGGTCGTCCGTACGAAGGTACCCAGACCGCTAAACTGAAAGTTACCAAAGGTGGTCCGCTGCCGTTCGCTTGGGACATCCTGTCCCCGCAGTTCCAGTACGGTTCCAAAGCTTACGTTAAACACCCGGCTGACATCCCGGACTACCTGAAACTGTCCTTCCCGGAAGGTTTCAAATGGGAACGTGTTATGAACTTCGAGGACGGTGGTGTTGTTACCGTTACCCAGGACTCCTCCCTGCAAGACGGTGAGTTCATCTACAAAGTTAAACTGCGTGGTACCAACTTCCCGTCCGACGGTCCGGTTATGCAGAAAAAAACCATGGGTTGGGAAGCTTCCACCGAACGTATGTACCCGGAGGACGGTGCTCTGAAAGGTGAAATCAAAATGCGTCTGAAACTGAAAGACGGTGGTCACTACGACGCTGAAGTTAAAACCACCTACATGGCTAAAAAACCGGTTCAGTTACCGGGTGCTTACAAAACCGACATCAAACTGGACATCACCTCCCACAACGAGGACTACACCATCGTTGAACAGTACGAACGTGCTGAAGGTCGTCACTCCACCGGTGCTTAAGCGATGTTGAAGACCATGA")



## Function to read in the sequencing file and trim it based on the phred score
def loadsequencing(file, threshold=0.9):
    seq = SeqIO.read(file, 'abi')

    maxq = np.max(seq.letter_annotations['phred_quality'])
    rolling = pd.Series(seq.letter_annotations['phred_quality']).rolling(window=20).mean()
    start = (rolling > maxq * threshold).idxmax()
    end = (rolling > maxq * threshold)[::-1].idxmax()

    return seq.seq[start:end], seq, start, end, np.mean(seq.letter_annotations['phred_quality'][start:end])

def align_reads(id_num, gene_name, target, forward, reverse, target_seq):
    forward_align = pairwise2.align.globalms(target_seq, forward,1,0,-1,-1, one_alignment_only=True, penalize_end_gaps=False)
    reverse_align = pairwise2.align.globalms(target_seq, reverse,1,0,-1,-1, one_alignment_only=True, penalize_end_gaps=False)
    forward_align = forward_align[0]
    reverse_align = reverse_align[0]
    for_raw = len(forward)
    rev_raw = len(reverse)
    target_length = len(target_seq)

    if for_raw <= target_length:
        for_score = for_raw - forward_align[2]
    else:
        for_score = target_length - forward_align[2]
    if rev_raw <= target_length:
        rev_score = rev_raw - reverse_align[2]
    else:
        rev_score = target_length - reverse_align[2]

    target_length = len(target_seq)

    if for_score == 0:
        print("forward perfect match")
        if rev_score == 0:
            print("reverse perfect match")
            outcome = "Perfect"
            manual = ""
        elif rev_score < 10:
            print("Only mut in rev")
            outcome = "Mutation: {} {}".format(for_score,rev_score)
            manual = "CHECK"
        else:
            print("bad rev read")
            outcome = "Bad Reverse"
            manual = "CHECK"
    elif for_score < 10:
        print("Mut in for read")
        if rev_score == 0:
            print("Only mut in for")
            outcome = "Mutation: {} {}".format(for_score,rev_score)
            manual = "CHECK"
        elif rev_score < 10:
            print("Mut in both reads")
            outcome = "Mutation: {} {}".format(for_score,rev_score)
            manual = "CHECK"
        else:
            print("bad rev read")
            outcome = "Mutation: {} {}".format(for_score,rev_score)
            manual = "CHECK"
    else:
        print("Bad for read")
        if rev_score == 0:
            print("Only bad for read")
            outcome = "Bad Forward"
            manual = "CHECK"
        elif rev_score < 10:
            print("Mut in rev read")
            outcome = "Mutation: {} {}".format(for_score,rev_score)
            manual = "CHECK"
        else:
            print("bad rev read")
            outcome = "Bad clone"
            manual = ""


    return [id_num, gene_name, target, for_raw, forward_align[2], rev_raw, reverse_align[2], target_length, outcome, manual]



old_insert = "ATGCGGTCTTCCGCATCGCCTGGCACGACAGGTTTCCCGACTGGAAAGCGGGCAGTGAGCGCAACGCAATTAATGTGAGTTAGCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATGTTGTGTGGAATTGTGAGCGGATAACAATTTCACACATACTAGAGAAAGAGGAGAAATACTAGATGGCTTCCTCCGAAGATGTTATCAAAGAGTTCATGCGTTTCAAAGTTCGTATGGAAGGTTCCGTTAACGGTCACGAGTTCGAAATCGAAGGTGAAGGTGAAGGTCGTCCGTACGAAGGTACCCAGACCGCTAAACTGAAAGTTACCAAAGGTGGTCCGCTGCCGTTCGCTTGGGACATCCTGTCCCCGCAGTTCCAGTACGGTTCCAAAGCTTACGTTAAACACCCGGCTGACATCCCGGACTACCTGAAACTGTCCTTCCCGGAAGGTTTCAAATGGGAACGTGTTATGAACTTCGAGGACGGTGGTGTTGTTACCGTTACCCAGGACTCCTCCCTGCAAGACGGTGAGTTCATCTACAAAGTTAAACTGCGTGGTACCAACTTCCCGTCCGACGGTCCGGTTATGCAGAAAAAAACCATGGGTTGGGAAGCTTCCACCGAACGTATGTACCCGGAGGACGGTGCTCTGAAAGGTGAAATCAAAATGCGTCTGAAACTGAAAGACGGTGGTCACTACGACGCTGAAGTTAAAACCACCTACATGGCTAAAAAACCGGTTCAGTTACCGGGTGCTTACAAAACCGACATCAAACTGGACATCACCTCCCACAACGAGGACTACACCATCGTTGAACAGTACGAACGTGCTGAAGGTCGTCACTCCACCGGTGCTTAAGCGATGTTGAAGACCATGA"

for forfile in glob.glob("{}/*{}*.ab1".format(SEQFILE_PATH,forward_primer)):
    print(forfile)

    if "BBF10K" in forfile:
        initials, order_number, plate_number, well_number, sample_name, sample_number, well_address = re.match(
            r'.*/([A-Z]+)_([0-9]+)-([0-9])([0-9]+)_([A-Za-z0-9]+)_([0-9]+)_M13-Forward---20-_?([A-H][0-9]{2}).ab1',
            forfile).groups()
        revfile = "{}/{}_{}-{}{}_{}_{}_{}_{}.ab1".format(os.path.dirname(forfile), initials, order_number, (int(plate_number)+1), well_number, sample_name, sample_number, reverse_primer, well_address)
        print(revfile)
        id_num = sample_name + "_" + sample_number
    elif "MMSYN" in forfile:
        initials, order_number, plate_number, well_number, hyphen, sample_name, primer_name, well_address = re.match(
            r'.*/([A-Z]+)_([0-9]+)-([0-9])([0-9]+)_(-?)([A-Za-z0-9_-]+)_([A-Za-z0-9-]+)_([A-H][0-9]{2}).ab1',
            fo`rfile).groups()
        revfile = "{}/{}_{}-{}{}_{}{}_{}_{}.ab1".format(os.path.dirname(forfile), initials, order_number, (int(plate_number)+1), well_number, hyphen, sample_name, reverse_primer, well_address)
        print("Sample_name", sample_name)
        id_num = dictionary[sample_name[:-2]]
        print(id_num)

    new_dir = "{}/{}/{}_seq_files".format(DATA_PATH,id_num,build_num)
#    if os.path.exists(new_dir):
#        print("{} already exists".format(new_dir))
#    else:
#        os.makedirs(new_dir)
#        shutil.copy(forfile, new_dir)
#        shutil.copy(revfile, new_dir)

    with open("{}/{}/{}.json".format(DATA_PATH,id_num,id_num),"r") as json_file:
        data = json.load(json_file)

    gene_seq = data["sequence"]["optimized_sequence"]
    gene_seq = Seq.Seq(gene_seq)
    gene_name = data["gene_name"]

    forward_untrim, _, _, _, forward_qual = loadsequencing(forfile)
    reverse_untrim, revseq, _, _, reverse_qual = loadsequencing(revfile)
    print("Quality", forward_qual, reverse_qual)

    if math.isnan(forward_qual) or math.isnan(reverse_qual):
        forward_untrim, _, _, _, forward_qual = loadsequencing(forfile, threshold=0.8)
        reverse_untrim, revseq, _, _, reverse_qual = loadsequencing(revfile, threshold=0.8)
        nan.append(id_num)
        print("Quality", forward_qual, reverse_qual)
    if len(forward_untrim) < 50 or len(reverse_untrim) < 50:
        forward_untrim, _, _, _, forward_qual = loadsequencing(forfile, threshold=0.8)
        reverse_untrim, revseq, _, _, reverse_qual = loadsequencing(revfile, threshold=0.8)
        small.append(id_num)

    forward = forward_untrim[forward_untrim.find('ATG'):] # Start at start codon
    reverse = reverse_untrim[reverse_untrim.find('TCA'):].reverse_complement() # Stop at stop

    row = align_reads(id_num, gene_name, "gene",forward,reverse,gene_seq)
    align_data.append(row)

    if (row[4] != min(row[3], row[7])) or (row[6] != min(row[5], row[7])):
        print("aligning to backbone")
        row = align_reads(id_num, gene_name, "backbone",forward,reverse,backbone_seq)
        align_data.append(row)

align_data = np.array(align_data)

array = pd.DataFrame({
    "Gene" : align_data[:,0],
    "Gene Name" : align_data[:,1],
    "Target" : align_data[:,2],
    "For Length" : align_data[:,3],
    "For Score" : align_data[:,4],
    "Rev Length" : align_data[:,5],
    "Rev Score" : align_data[:,6],
    "Gene Length" : align_data[:,7],
    "Outcome" : align_data[:,8],
    "Manual" : align_data[:,9]
})

array = array[["Gene","Gene Name","Target","For Length","For Score","Rev Length","Rev Score","Gene Length","Outcome","Manual"]]

array.to_csv("{}/{}/{}_alignment_results.csv".format(BUILDS_PATH,build_num,build_num))

print(array)
print("nan sequences:", nan)
print("small sequences:", small)

new_dir = "{}/{}/{}_fasta_files".format(BUILDS_PATH,build_num,build_num)
print("new_dir", new_dir)
if os.path.exists(new_dir):
    print("{} already exists".format(new_dir))
else:
    os.makedirs(new_dir)
    print("made directory")

for index,row in array.iterrows():
    if row["Outcome"] == "Perfect":
        continue
    for fasta in glob.glob("../data/{}/{}.fasta".format(row["Gene"],row["Gene"])):
        shutil.copy(fasta, new_dir)

print()
stop = datetime.now()
runtime = stop - start
print("Total runtime is: ", runtime)


#
