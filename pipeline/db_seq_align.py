## FOR GIT REPOSITORY -- RUN ALIGNMENTS AND STATE RESULT
## Goes through each directory with .abi files, runs an alignment and then appends the json file with the results

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import json
import os
import glob
import re
import math

import shutil

from Bio import pairwise2
from Bio import AlignIO
from Bio import Align
from Bio import Seq
from Bio import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqIO

## Dependencies:

#State the primers
forward_primer = "M13-Forward---20-"
reverse_primer = "M13-Reverse"

# Paths for backbone and data files
DATA_PATH = "../data/BBF10K_*/*{}*.ab1".format(forward_primer)
DATA_FASTA_PATH = "../data/{}/{}.fasta"
DATA_JSON_PATH = "../data/{}/{}.json"
BACKBONE_PATH = '../sequencing_files/popen_v1-1_backbone.fasta'


#State the vector sequence that the genes are getting built into
backbone_seq = SeqIO.read(BACKBONE_PATH, 'fasta')
backbone_seq = backbone_seq.seq

# Create a dictionary to link the gene name to the corresponding id number
data = pd.read_csv("./testing/data_testing/10K_CDS.csv")
dictionary = dict(zip(data['gene_name'], data['idnum']))

## Function to read in the sequencing file and trim it based on the phred score
def loadsequencing(file):
    seq = SeqIO.read(file, 'abi')

    maxq = np.max(seq.letter_annotations['phred_quality'])
    rolling = pd.Series(seq.letter_annotations['phred_quality']).rolling(window=20).mean()
    start = (rolling > maxq * 0.9).idxmax()
    end = (rolling > maxq * 0.9)[::-1].idxmax()

    return seq.seq[start:end], seq, start, end, np.mean(seq.letter_annotations['phred_quality'][start:end])

#Iterates through every .ab1 file in the specified directories
for seqfile in glob.glob(DATA_PATH):
    print("Running file :{}".format(seqfile))

    i = 0
## This was my way of counting the number of seq files in a folder in case
## there were more than one set of reads to stop it from just overwriting it
#    direct = os.path.dirname(seqfile)
#    for files in glob.glob("{}/*".format(direct)):
#        i = i + 1
#    if i == 4:
#        print("Only four abi files")
#    else:
#        print("Directory {} has {} files".format(direct, i))
#        continue
#    print(i)

    print(seqfile[37])

    goodF = 0
    goodR = 0
    rfpF = 0
    rfpR = 0

    initials = 0

    # Included an if statement to handle the exception cases with (-)MMSYN1
    # The "37" is hard wired based on the path and file name.
    if seqfile[37] == "-":
        initials, order_number, plate_number, well_number, hyphen, sample_name, primer_name, well_address = re.match(
            r'.*/([A-Z]+)_([0-9]+)-([0-9])([0-9]+)_(-)([A-Za-z0-9_-]+)_([A-Za-z0-9-]+)_([A-H][0-9]{2}).ab1',
            seqfile).groups()
        revfile = "{}/{}_{}-{}{}_{}{}_{}_{}.ab1".format(os.path.dirname(seqfile), initials, order_number, (int(plate_number)+1), well_number, hyphen, sample_name, reverse_primer, well_address)
    else:
        initials, order_number, plate_number, well_number, sample_name, primer_name, well_address = re.match(
            r'.*/([A-Z]+)_([0-9]+)-([0-9])([0-9]+)_([A-Za-z0-9_-]+)_([A-Za-z0-9-]+)_([A-H][0-9]{2}).ab1',
            seqfile).groups()
        revfile = "{}/{}_{}-{}{}_{}_{}_{}.ab1".format(os.path.dirname(seqfile), initials, order_number, (int(plate_number)+1), well_number, sample_name, reverse_primer, well_address)

    #if initials == 0:
    #    print("no seq files")
    #    continue

    print(sample_name)

    # Eliminates the fragment number
    split = sample_name.rsplit('_')
    gene = split[0] + "_" + split[1]
    print(gene)

    # stands for ID# and allows you to use the new ID number or the gene name that was given previously
    idnum = dictionary[gene]

    gene_seq = SeqIO.read(DATA_FASTA_PATH.format(idnum,idnum),"fasta")
    gene_seq = gene_seq.seq

    # We assume that for every forward there is a reverse
    forward_untrim, _, _, _, forward_qual = loadsequencing(seqfile)
    reverse_untrim, revseq, _, _, reverse_qual = loadsequencing(revfile)
    print("Quality", forward_qual, reverse_qual)

    if math.isnan(forward_qual) or math.isnan(reverse_qual):
        if math.isnan(forward_qual):
            print("forward is nan")
        if math.isnan(reverse_qual):
            print("reverse is nan")
        continue

    forward = forward_untrim[forward_untrim.find('ATG'):] # Start at start codon
    reverse = reverse_untrim[reverse_untrim.find('TCA'):].reverse_complement() # Stop at stop

    forward_align = pairwise2.align.globalxx(gene_seq, forward, one_alignment_only=True)
    reverse_align = pairwise2.align.globalxx(gene_seq, reverse, one_alignment_only=True)

    forward_align = forward_align[0]
    reverse_align = reverse_align[0]

    if len(forward_align) == 0 or len(reverse_align) == 0:
        if len(forward_align) == 0:
            print("Forward doesn't align")

        if len(reverse_align) == 0:
            print("Reverse doesn't align")

        continue

    ## Checks the forward alignment
    if (forward_align[2] != min(len(forward), len(gene_seq))):
        print("Forward doesn't match gene: scores are {} and {}".format(forward_align[2],len(forward)))
        forward_align_back = pairwise2.align.globalxx(backbone_seq, forward_untrim, one_alignment_only=True)
        forward_align_back = forward_align_back[0]

        ## Checks the forward sequence against the backbone sequence
        if (forward_align_back[2] == min(len(forward_untrim), len(backbone_seq))):
            print("Forward aligns to backbone: scores are {} and {}".format(forward_align_back[2],len(forward_untrim)))
            rfpF = rfpF + 1
        else:
            print("Forward doesn't align to backbone: scores are {} and {}".format(forward_align_back[2],len(forward_untrim)))
    else:
        print("Forward matches gene: scores are {} and {}".format(forward_align[2],len(forward)))
        goodF = goodF + 1

    ## Checks the reverse alignment
    if (reverse_align[2] != min(len(reverse), len(gene_seq))):
        print("Reverse doesn't match gene: scores are {} and {}".format(reverse_align[2],len(reverse)))
        reverse_align_back = pairwise2.align.globalxx(backbone_seq, reverse_untrim, one_alignment_only=True)
        reverse_align_back = reverse_align_back[0]

        ## Checks the reverse sequence against the backbone sequence
        if (reverse_align_back[2] == min(len(reverse_untrim), len(backbone_seq))):
            print("Reverse aligns to backbone: scores are {} and {}".format(reverse_align_back[2],len(reverse_untrim)))
            rfpR = 1
        else:
            print("Reverse doesn't align to backbone: scores are {} and {}".format(reverse_align_back[2],len(reverse_untrim)))
    else:
        print("Reverse matches gene: scores are {} and {}".format(reverse_align[2],len(reverse)))
        goodR = 1

    ## Open the json file to store the saved data
    for file2 in glob.glob("../Data/{}/{}.json".format(idnum,idnum)):
        print(file2)

        with open(file2,"r") as json_file:
            data = json.load(json_file)

        ## Determine if both of the sequences align to the same target and modify the data for the json file
        if goodF + goodR == 2:
            print("good")
            data["status"]["build_complete"] = "Good_Sequence"
        elif rfpF + rfpR == 2:
            print("rfp")
            data["status"]["build_complete"] = "Original_Vector_Sequence"
        else:
            print("unknown")
            data["status"]["build_complete"] = "Unknown_Sequence"

        ## Write to the json file
        with open(file2,'w') as json_file:
            json.dump(data,json_file,indent=2)
