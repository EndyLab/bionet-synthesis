from datetime import datetime
start = datetime.now()

import numpy as np
import pandas as pd

import json
import os
import glob
import re
import math

import shutil
from config import *

from Bio import pairwise2
from Bio import Seq
from Bio import SeqIO
from Bio.Blast import NCBIXML

## ============================================
## ESTABLISH PATHS AND NAMES
## ============================================
PIPELINE_PATH = BASE_PATH + "/pipeline"
BUILDS_PATH = BASE_PATH + "/builds"
DATA_PATH = BASE_PATH + "/data"
BACKBONE_PATH = BASE_PATH + "/sequencing_files/popen_v1-1_backbone.fasta"
DICTIONARY_PATH = PIPELINE_PATH + "/testing/data_testing/10K_CDS.csv"
DATABASE_PATH = BASE_PATH + "/raw_files/BLAST_db/current_BLAST_db.fsa"

forward_primer = "M13-Forward---20-"
reverse_primer = "M13-Reverse"

## ============================================
## TAKE IN BUILD INFORMATION
## ============================================

def choose_build(build_num,path):
    '''Checks to make sure that the build exists and if not asks for a different one'''
    if glob.glob(path):
        return build_num
    else:
        build_num = input("Not a valid build number. Please try again: ")
        build_num = "build{}".format(str(build_num).zfill(3))
        path = "{}/builds/{}/{}_20*.csv".format(BASE_PATH,build_num,build_num)
        return choose_build(build_num,path)

build_num = input("Enter build number: ")
build_num = "build{}".format(str(build_num).zfill(3))
path = "{}/builds/{}/{}_20*.csv".format(BASE_PATH,build_num,build_num)
build_num = choose_build(build_num,path)
SEQFILE_PATH = "{}/{}/{}_seq_files".format(BUILDS_PATH,build_num,build_num)

# Create a dictionary to link the gene name to the corresponding id number
data = pd.read_csv(DICTIONARY_PATH)
dictionary = dict(zip(data['gene_name'], data['idnum']))

## ============================================
## GENERATE BLAST DATABASE
## ============================================
# Create a multientry FASTA file with all seqs from database
db_counter = 1
with open(DATABASE_PATH,"w+") as fsa:
    for file in glob.glob(DATA_PATH + "/*/*.json"):
        with open(file) as json_file:
            data = json.load(json_file)
        gene_id = data["gene_id"]
        sequence = data["sequence"]["optimized_sequence"]
        fsa.write(">{}|BLAST_db|{}\n{}\n".format(gene_id,db_counter,sequence))
        db_counter += 1
        #if db_counter > 10:
        #    break
    # fsa.close()
# Convert the FSA file into a BLAST database
os.system("makeblastdb -in {} -parse_seqids -dbtype nucl".format(DATABASE_PATH))

count = 0
align_data = []
unknown_data = []
nan = []
small = []
check = []

# REVIEW: Change this to take in the fasta file for the backbone
# REVIEW: When using the full backbone it would throw off the alignment so this is just the insert
backbone_seq = Seq.Seq("ATGCGGTCTTCCGCATCGCCTGGCACGACAGGTTTCCCGACTGGAAAGCGGGCAGTGAGCGCAACGCAATTAATGTGAGTTAGCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATGTTGTGTGGAATTGTGAGCGGATAACAATTTCACACATACTAGAGAAAGAGGAGAAATACTAGATGGCTTCCTCCGAAGATGTTATCAAAGAGTTCATGCGTTTCAAAGTTCGTATGGAAGGTTCCGTTAACGGTCACGAGTTCGAAATCGAAGGTGAAGGTGAAGGTCGTCCGTACGAAGGTACCCAGACCGCTAAACTGAAAGTTACCAAAGGTGGTCCGCTGCCGTTCGCTTGGGACATCCTGTCCCCGCAGTTCCAGTACGGTTCCAAAGCTTACGTTAAACACCCGGCTGACATCCCGGACTACCTGAAACTGTCCTTCCCGGAAGGTTTCAAATGGGAACGTGTTATGAACTTCGAGGACGGTGGTGTTGTTACCGTTACCCAGGACTCCTCCCTGCAAGACGGTGAGTTCATCTACAAAGTTAAACTGCGTGGTACCAACTTCCCGTCCGACGGTCCGGTTATGCAGAAAAAAACCATGGGTTGGGAAGCTTCCACCGAACGTATGTACCCGGAGGACGGTGCTCTGAAAGGTGAAATCAAAATGCGTCTGAAACTGAAAGACGGTGGTCACTACGACGCTGAAGTTAAAACCACCTACATGGCTAAAAAACCGGTTCAGTTACCGGGTGCTTACAAAACCGACATCAAACTGGACATCACCTCCCACAACGAGGACTACACCATCGTTGAACAGTACGAACGTGCTGAAGGTCGTCACTCCACCGGTGCTTAAGCGATGTTGAAGACCATGA")

## ============================================
## DEFINE INITIAL FUNCTIONS
## ============================================

## Function to read in the sequencing file and trim it based on the phred score
def loadsequencing(file, threshold=0.9):
    '''Load sequencing reads and trim them based on the specified threshold'''
    seq = SeqIO.read(file, 'abi')

    maxq = np.max(seq.letter_annotations['phred_quality'])
    rolling = pd.Series(seq.letter_annotations['phred_quality']).rolling(window=20).mean()
    start = (rolling > maxq * threshold).idxmax()
    end = (rolling > maxq * threshold)[::-1].idxmax()

    return seq.seq[start:end], seq, start, end, np.mean(seq.letter_annotations['phred_quality'][start:end])

align_data = []
def align_reads(forward, reverse, target_seq):
    '''Generates a forward and reverse alignments and then determines the result'''
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

    if for_raw < 100 and rev_raw < 100:
        outcome = "Bad_reads"
    elif for_score == 0:
        if rev_score == 0:
            outcome = "Perfect"
        elif rev_score < 10:
            outcome = "Mutation: {} {}".format(for_score,rev_score)
        else:
            outcome = "Bad Reverse"
    elif for_score < 10:
        if rev_score == 0:
            outcome = "Mutation: {} {}".format(for_score,rev_score)
        elif rev_score < 10:
            outcome = "Mutation: {} {}".format(for_score,rev_score)
        else:
            outcome = "Mutation: {} {}".format(for_score,rev_score)
    else:
        if rev_score == 0:
            outcome = "Bad Forward"
        elif rev_score < 10:
            outcome = "Mutation: {} {}".format(for_score,rev_score)
        else:
            outcome = "Bad clone"
    return [outcome, for_raw, forward_align[2], rev_raw, reverse_align[2], target_length]

def verify_sequence(id_num,forward,reverse,gene_seq,backbone_seq):
    '''
    Runs alignments against both the target sequence and the backbone and
    returns the result.
    '''
    g_res = align_reads(forward,reverse,gene_seq)
    b_res = align_reads(forward,reverse,backbone_seq)

    if g_res[0] == "Perfect" and b_res[0] == "Bad clone":
        final = "Good_sequence"
    elif "Mutation" in g_res[0] and b_res[0] == "Bad clone":
        final = "Point_mutation"
    elif g_res[0] == "Bad clone" and b_res[0] == "Perfect":
        final = "Original_vector_sequence"
    elif g_res[0] == "Bad clone" and "Mutation" in b_res[0]:
        final = "Original_vector_sequence"
    elif g_res[0] == "Bad clone" and b_res[0] == "Bad clone":
        final = "Unknown_sequence"
    elif g_res[0] == "Bad_reads" and b_res[0] == "Bad_reads":
        final = "Bad_reads"
    else:
        final = "CHECK"
    return [final] + g_res + b_res

def blast_seq(name,build_num,sequence,direction, E_VALUE_THRESH=0.04):
    '''Runs a BLAST search with a read against the database'''
    # Create fasta files to run the BLAST search
    query = "{}/{}/{}.fasta".format(BUILDS_PATH,build_num,direction)
    with open(query,"w+") as seq_file:
        seq_file.write(">{}\n{}".format(name,sequence))

    # Run the BLAST search
    output = "{}/{}/{}_results.xml".format(BUILDS_PATH,build_num,direction)
    os.system("blastn -query {} -db {}  -out {} -evalue 0.001 -outfmt 5".format(query,DATABASE_PATH,output))

    # Take in the blast results
    result_handle = open(output)
    blast_record = NCBIXML.read(result_handle)
    os.remove(output)
    os.remove(query)

    # Return the top BLAST hit
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                hit_id = alignment.title.split("|")
                return hit_id[0]

def align_unknown(name,build_num,forward_sequence,reverse_sequence):
    '''
    Runs a BLAST search on both reads and determines if they share a common hit,
    if they do it preforms an alignment of the reads on the target.
    '''
    print("Using align unknown on {}".format(name))
    for_hit = blast_seq(name,build_num,forward_sequence,'forward')
    rev_hit = blast_seq(name,build_num,reverse_sequence,'reverse')
    if for_hit == rev_hit:
        print("match")
        with open("{}/{}/{}.json".format(DATA_PATH,for_hit,for_hit),"r") as json_file:
            data = json.load(json_file)
        target = data["sequence"]["optimized_sequence"]
        target = Seq.Seq(target)
        new_row = align_reads(forward,reverse,target)
        new_row = [name,for_hit] + new_row
        return new_row
    else:
        return [name,"no hit", "empty", "empty", "empty", "empty", "empty", "empty"]

## ============================================
## TAKE IN ALL OF THE SEQUENCING FILES
## ============================================
for forfile in glob.glob("{}/*{}*.ab1".format(SEQFILE_PATH,forward_primer)):
    # There has been inconsistency in the naming of samples so this is made to account for them
    if "BBF10K" in forfile:
        initials, order_number, plate_number, well_number, sample_name, sample_number, well_address = re.match(
            r'.*/([A-Z]+)_([0-9]+)-([0-9])([0-9]+)_([A-Za-z0-9]+)_([0-9]+)_M13-Forward---20-_?([A-H][0-9]{2}).ab1',
            forfile).groups()
        revfile = "{}/{}_{}-{}{}_{}_{}_{}_{}.ab1".format(os.path.dirname(forfile), initials, order_number, (int(plate_number)+1), well_number, sample_name, sample_number, reverse_primer, well_address)
        id_num = sample_name + "_" + sample_number
        unknown = False
    elif "MMSYN" in forfile:
        initials, order_number, plate_number, well_number, hyphen, sample_name, primer_name, well_address = re.match(
            r'.*/([A-Z]+)_([0-9]+)-([0-9])([0-9]+)_(-?)([A-Za-z0-9_-]+)_([A-Za-z0-9-]+)_([A-H][0-9]{2}).ab1',
            forfile).groups()
        revfile = "{}/{}_{}-{}{}_{}{}_{}_{}.ab1".format(os.path.dirname(forfile), initials, order_number, (int(plate_number)+1), well_number, hyphen, sample_name, reverse_primer, well_address)
        id_num = dictionary[sample_name[:-2]]
        unknown = False
    elif "Unk" in forfile:
        initials, order_number, plate_number, well_number, sample_name, sample_number, well_address = re.match(
            r'.*/([A-Z]+)_([0-9]+)-([0-9])([0-9]+)_([A-Za-z0-9]+)_([0-9]+)_M13-Forward---20-_?([A-H][0-9]{2}).ab1',
            forfile).groups()
        revfile = "{}/{}_{}-{}{}_{}_{}_{}_{}.ab1".format(os.path.dirname(forfile), initials, order_number, (int(plate_number)+1), well_number, sample_name, sample_number, reverse_primer, well_address)
        id_num = sample_name + "_" + sample_number
        unknown = True

    # Record the sequencing file names so that they can be used to realign at a later time
    for_seq_file = forfile.split("/")[-1]
    print("Forward: ",for_seq_file)
    rev_seq_file = revfile.split("/")[-1]
    print("Reverse: ",rev_seq_file)

    # Generates a new directory for each gene with their reads for a specific build
    #    new_dir = "{}/{}/{}_seq_files".format(DATA_PATH,id_num,build_num)
    #    if os.path.exists(new_dir):
    #        print("{} already exists".format(new_dir))
    #    else:
    #        os.makedirs(new_dir)
    #        shutil.copy(forfile, new_dir)
    #        shutil.copy(revfile, new_dir)

    # Trim the reads and check their quality
    forward_untrim, _, _, _, forward_qual = loadsequencing(forfile)
    reverse_untrim, revseq, _, _, reverse_qual = loadsequencing(revfile)
    print(id_num, "Quality", forward_qual, reverse_qual)

    # Check that the reads have sufficient length and if not decrease the threshold and reload
    if math.isnan(forward_qual) or math.isnan(reverse_qual):
        forward_untrim, _, _, _, forward_qual = loadsequencing(forfile, threshold=0.8)
        reverse_untrim, revseq, _, _, reverse_qual = loadsequencing(revfile, threshold=0.8)
        nan.append(id_num)
        print("Quality", forward_qual, reverse_qual)
    if len(forward_untrim) < 50 or len(reverse_untrim) < 50:
        forward_untrim, _, _, _, forward_qual = loadsequencing(forfile, threshold=0.8)
        reverse_untrim, revseq, _, _, reverse_qual = loadsequencing(revfile, threshold=0.8)
        small.append(id_num)

    # REVIEW: Not all sequences are CDS's and so we need to grab bases from the beginning and
    # of the target sequence or something to trim it.

    # Trim the reads to start at the start or stop codons
    forward = forward_untrim[forward_untrim.find('ATG'):] # Start at start codon
    reverse = reverse_untrim[reverse_untrim.find('TCA'):].reverse_complement() # Stop at stop

    if unknown:
        print("_____________UNKNOWN__________________")
        row = align_unknown(id_num,build_num,forward,reverse) + [well_address,for_seq_file,rev_seq_file]
        unknown_data.append(row)

    else:
        with open("{}/{}/{}.json".format(DATA_PATH,id_num,id_num),"r") as json_file:
            data = json.load(json_file)

        gene_seq = data["sequence"]["optimized_sequence"]
        gene_seq = Seq.Seq(gene_seq)
        gene_name = data["gene_name"]

        row = verify_sequence(id_num,forward,reverse,gene_seq,backbone_seq)
        if row[0] == "Unknown_sequence":
            row = align_unknown(id_num,build_num,forward,reverse) + [well_address,for_seq_file,rev_seq_file]
            unknown_data.append(row)

        else:
            row = [id_num] + [gene_name] + row + [well_address,for_seq_file,rev_seq_file]
            align_data.append(row)

# Build a dataframe containing all of the results from the unknown sequences
unknown_data = np.array(unknown_data)
unknown_df = pd.DataFrame({
    "Gene ID" : unknown_data[:,0],
    "Hit ID" : unknown_data[:,1],
    "Hit Result" : unknown_data[:,2],
    "Hit For Length" : unknown_data[:,3],
    "Hit For Score" : unknown_data[:,4],
    "Hit Rev Length" : unknown_data[:,5],
    "Hit Rev Score" : unknown_data[:,6],
    "Hit Length" : unknown_data[:,7],
    "Well" : unknown_data[:,8],
    "For Read" : unknown_data[:,9],
    "Rev Read" : unknown_data[:,10]
})
unknown_df = unknown_df[["Gene ID","Hit ID","Hit Result","Well","Hit Result","Hit For Length","Hit For Score","Hit Rev Length","Hit Rev Score","Hit Length","For Read","Rev Read"]]
unknown_df.to_csv("{}/{}/{}_alignment_results-unknown.csv".format(BUILDS_PATH,build_num,build_num))

# Split apart the results of the unknown sequences into the intended targets and the actual hits
intended_unknown = pd.DataFrame({
    "Gene ID" : unknown_data[:,0],
    "Well" : unknown_data[:,8],
    "Outcome" : "Failed",
    "For Read" : unknown_data[:,9],
    "Rev Read" : unknown_data[:,10],
    "Unintended" : True
})
actual_unknown = pd.DataFrame({
    "Gene ID" : unknown_data[:,1],
    "Well" : unknown_data[:,8],
    "Outcome" : unknown_data[:,2],
    "For Read" : unknown_data[:,9],
    "Rev Read" : unknown_data[:,10],
    "Unintended" : True
})
complete_unknown = pd.concat([intended_unknown,actual_unknown])

# Build a dataframe containing all of the results from the normal sequences
align_data = np.array(align_data)
array = pd.DataFrame({
    "Gene ID" : align_data[:,0],
    "Gene Name" : align_data[:,1],
    "Outcome" : align_data[:,2],
    "Gene Result" : align_data[:,3],
    "Gene For Length" : align_data[:,4],
    "Gene For Score" : align_data[:,5],
    "Gene Rev Length" : align_data[:,6],
    "Gene Rev Score" : align_data[:,7],
    "Gene Length" : align_data[:,8],
    "Template Result" : align_data[:,9],
    "Template For Length" : align_data[:,10],
    "Template For Score" : align_data[:,11],
    "Template Rev Length" : align_data[:,12],
    "Template Rev Score" : align_data[:,13],
    "Template Length" : align_data[:,14],
    "Well" : align_data[:,15],
    "For Read" : align_data[:,16],
    "Rev Read" : align_data[:,17]
})
array = array[["Gene ID","Gene Name","Outcome","Well","Gene Result","For Read","Rev Read","Gene For Length","Gene For Score","Gene Rev Length","Gene Rev Score","Gene Length","Template Result","Template For Length","Template For Score" ,"Template Rev Length","Template Rev Score" ,"Template Length"]]
array.to_csv("{}/{}/{}_alignment_results-array.csv".format(BUILDS_PATH,build_num,build_num))
outcomes = pd.DataFrame({
    "Well" : array["Well"],
    "Gene ID" : array["Gene ID"],
    "Outcome" : array["Outcome"],
    "Unintended" : False,
    "For Read" : align_data[:,16],
    "Rev Read" : align_data[:,17]
})

# Combine all of the dataframes to create a comprehensive dataframe to update the database
complete = pd.concat([outcomes,complete_unknown])
other = []
## ============================================
## UPDATE THE DATABASE
## ============================================
for index, row in complete.iterrows():
    if "BBF10K_" not in row["Gene ID"]:
        print("Skipped: ", row["Gene ID"])
        continue
    with open("{}/{}/{}.json".format(DATA_PATH,row["Gene ID"],row["Gene ID"]),"r") as json_file:
        data = json.load(json_file)
    for build in reversed(data["status"]["build_attempts"]):

        # Converts 'A2' -> 'A02'
        if len(build["build_well"]) == 2:
            build["build_well"] = str(build["build_well"][0]+"0"+build["build_well"][1])

        # Updates the build attempt with the outcomes and seq files
        if build["build_number"] == build_num:
            if build["build_well"] == row["Well"]:
                build["build_outcome"] = row["Outcome"]
                build["forward_read"] = row["For Read"]
                build["reverse_read"] = row["Rev Read"]

            # Creates a new dict for the unintended build
            elif row["Unintended"]:
                data["status"]["build_attempts"].append({
                    "build_well" : row["Well"],
                    "build_number" : build_num,
                    "build_outcome" : row["Outcome"],
                    "forward_read" : row["For Read"],
                    "reverse_read" : row["Rev Read"]
                })
            else:
                other.append(row["Gene ID"])
    if row["Outcome"] == "Good_sequence" or row["Outcome"] == "Perfect":
        data["status"]["build_complete"] = True
    else:
        data["status"]["build_complete"] = False
    with open("{}/{}/{}.json".format(DATA_PATH,row["Gene ID"],row["Gene ID"]),"w") as json_file:
        json.dump(data,json_file,indent=2)

complete.to_csv("{}/{}/{}_alignment_results-complete.csv".format(BUILDS_PATH,build_num,build_num))


print(array)
print("nan sequences:", nan)
print("small sequences:", small)
print("other",other)

input("continue?")

# REVIEW: Currently not pulling all of the FASTA files because soon there will only be Genbank files

# new_dir = "{}/{}/{}_fasta_files".format(BUILDS_PATH,build_num,build_num)
# if os.path.exists(new_dir):
#     print("{} already exists".format(new_dir))
# else:
#     os.makedirs(new_dir)
#     print("made directory")
#
# for index,row in array.iterrows():
#     for fasta in glob.glob("{}/{}/{}.fasta".format(DATA_PATH,row["Gene ID"],row["Gene ID"])):
#         shutil.copy(fasta, new_dir)

# print()
# print(array['Outcome'].value_counts())
# print()


stop = datetime.now()
runtime = stop - start
print("Total runtime is: ", runtime)
