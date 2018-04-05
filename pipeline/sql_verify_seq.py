import argparse
import sys

import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

import os
import re
import math
import glob
import json
import numpy as np
import pandas as pd
from datetime import datetime

from Bio import pairwise2
from Bio import Seq
from Bio import SeqIO
from Bio.Blast import NCBIXML

from config import *
from db_config import *

# TODO: Turn into a function

def verify_seq():
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
    # Present all available plates to resuspend
    print("Which build would you like to analyze?")
    build_names = []
    for index,build in enumerate(session.query(Build).filter(Build.status == 'sequencing')):
        print("{}. {}".format(index,build.build_name))
        build_names.append(build.build_name)

    # Asks the user for a number corresponding to the plate they want to resuspend
    number = int(input("Enter build: "))
    target = session.query(Build).filter(Build.build_name == build_names[number]).first()
    seq_plate = [plate for plate in target.plates if plate.plate_type == 'seq_plate'][0]
    print("Will analyze build: ",target.build_name)

    path = "{}/builds/{}/{}_20*.csv".format(BASE_PATH,target.build_name,target.build_name)
    SEQFILE_PATH = "{}/{}/{}_seq_files".format(BUILDS_PATH,target.build_name,target.build_name)

    # Create a dictionary to link the gene name to the corresponding id number
    data = pd.read_csv(DICTIONARY_PATH)
    dictionary = dict(zip(data['idnum'],data['gene_name']))

    ## ============================================
    ## GENERATE BLAST DATABASE
    ## ============================================
    # REVIEW: Change this to take in the fasta file for the backbone
    # REVIEW: When using the full backbone it would throw off the alignment so this is just the insert
    backbone = "ATGCGGTCTTCCGCATCGCCTGGCACGACAGGTTTCCCGACTGGAAAGCGGGCAGTGAGCGCAACGCAATTAATGTGAGTTAGCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATGTTGTGTGGAATTGTGAGCGGATAACAATTTCACACATACTAGAGAAAGAGGAGAAATACTAGATGGCTTCCTCCGAAGATGTTATCAAAGAGTTCATGCGTTTCAAAGTTCGTATGGAAGGTTCCGTTAACGGTCACGAGTTCGAAATCGAAGGTGAAGGTGAAGGTCGTCCGTACGAAGGTACCCAGACCGCTAAACTGAAAGTTACCAAAGGTGGTCCGCTGCCGTTCGCTTGGGACATCCTGTCCCCGCAGTTCCAGTACGGTTCCAAAGCTTACGTTAAACACCCGGCTGACATCCCGGACTACCTGAAACTGTCCTTCCCGGAAGGTTTCAAATGGGAACGTGTTATGAACTTCGAGGACGGTGGTGTTGTTACCGTTACCCAGGACTCCTCCCTGCAAGACGGTGAGTTCATCTACAAAGTTAAACTGCGTGGTACCAACTTCCCGTCCGACGGTCCGGTTATGCAGAAAAAAACCATGGGTTGGGAAGCTTCCACCGAACGTATGTACCCGGAGGACGGTGCTCTGAAAGGTGAAATCAAAATGCGTCTGAAACTGAAAGACGGTGGTCACTACGACGCTGAAGTTAAAACCACCTACATGGCTAAAAAACCGGTTCAGTTACCGGGTGCTTACAAAACCGACATCAAACTGGACATCACCTCCCACAACGAGGACTACACCATCGTTGAACAGTACGAACGTGCTGAAGGTCGTCACTCCACCGGTGCTTAAGCGATGTTGAAGACCATGA"
    backbone_seq = Seq.Seq(backbone)
    # Create a multientry FASTA file with all seqs from database
    db_counter = 1
    with open(DATABASE_PATH,"w+") as fsa:
        for part in session.query(Part):
            fsa.write(">{}|BLAST_db|{}\n{}\n".format(part.part_id,db_counter,part.seq))
            db_counter += 1
        fsa.write(">{}|BLAST_db|{}\n{}\n".format('backbone',db_counter,backbone))



    # Convert the FSA file into a BLAST database
    os.system("makeblastdb -in {} -parse_seqids -dbtype nucl".format(DATABASE_PATH))

    count = 0
    align_data = []
    unknown_data = []

    nan = []
    small = []
    check = []
    no_seq = []
    no_hit = []
    strange = []
    outcomes = []

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
                outcome = "Bad Reverse"
        else:
            if rev_score == 0:
                outcome = "Bad Forward"
            elif rev_score < 10:
                outcome = "Bad Forward".format(for_score,rev_score)
            else:
                outcome = "Bad clone"
        return outcome

    def verify_sequence(id_num,forward,reverse,gene_seq,backbone_seq,strange):
        '''
        Runs alignments against both the target sequence and the backbone and
        returns the result.
        '''
        g_res = align_reads(forward,reverse,gene_seq)
        b_res = align_reads(forward,reverse,backbone_seq)

        if g_res == "Perfect":
            if b_res == "Bad clone":
                final = "sequence_confirmed"
            elif b_res == "Bad_reads":
                final = "sequence_failure"
            elif "Bad" in b_res:
                final = "sequence_failure"
            else:
                print("==== Not Sure ====",g_res,b_res)
                strange += [[id_num,g_res,b_res]]
                final = "sequence_confirmed"
        elif "Mutation" in g_res:
            if b_res == "Bad clone":
                final = "cloning_mutation"
            elif "Mutation" in b_res:
                final = "cloning_failure"
            elif "Bad" in b_res:
                final = "cloning_mutation"
            else:
                print("==== Not Sure ====",g_res,b_res)
                strange += [[id_num,g_res,b_res]]
                final = "cloning_mutation"
        elif g_res == "Bad clone":
            if b_res == "Perfect":
                final = "cloning_failure"
            elif "Mutation" in b_res:
                final = "cloning_failure"
            elif b_res == "Bad clone":
                final = "Unknown_sequence"
            elif "Bad" in b_res:
                final = "cloning_mutation"
            else:
                print("==== Not Sure ====",g_res,b_res)
                strange += [[id_num,g_res,b_res]]
                final = "Unknown_sequence"
        elif g_res == "Bad_reads" and b_res == "Bad_reads":
            final = "sequence_failure"
        elif "Bad" in g_res:
            final = "sequence_failure"
        else:
            final = "{} - {}".format(g_res,b_res)
        return final,strange

    def blast_seq(name,sequence,direction, E_VALUE_THRESH=0.04):
        '''Runs a BLAST search with a read against the database'''
        # Create fasta files to run the BLAST search
        query = BASE_PATH + "/raw_files/BLAST_db/{}.fasta".format(direction)
        with open(query,"w+") as seq_file:
            seq_file.write(">{}\n{}".format(name,sequence))

        # Run the BLAST search
        output = BASE_PATH + "/raw_files/BLAST_db/results.xml"
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

    def align_unknown(name,forward_sequence,reverse_sequence,well,plate):
        '''
        Runs a BLAST search on both reads and determines if they share a common hit,
        if they do it preforms an alignment of the reads on the target.
        '''
        print("Using align unknown on {}".format(name))
        for_hit = blast_seq(name,forward_sequence,'forward')
        rev_hit = blast_seq(name,reverse_sequence,'reverse')
        if for_hit == rev_hit:
            target_part = session.query(Part).filter(Part.part_id == for_hit).first()
            try:
                target = Seq.Seq(target_part.seq)
                print("Matched: ",target_part.part_id)
                g_res = align_reads(forward,reverse,target)
                print('Unknown align: ',g_res)
                if g_res == "Perfect":
                    final = "sequence_confirmed"
                    new_well = Well('seq_plate',target_part,well.address,seq_outcome=final)
                    new_wells.append(new_well)
                    target_part.wells.append(new_well)
                elif "Mutation" in g_res:
                    final = "cloning_mutation"
                    new_well = Well('seq_plate',target_part,well.address,seq_outcome=final)
                    new_wells.append(new_well)
                    target_part.wells.append(new_well)
                elif g_res == "Bad clone":
                    final = "cloning_failure"
                elif "Bad" in g_res:
                    final = "sequence_failure"
                else:
                    final = "{} - {}".format(g_res,b_res)
                print('finished try')
                return 'cloning_error'
            except:
                print("Couldn't find: ",for_hit)
                no_hit.append(name)
                no_hit.append(for_hit)
                return 'cloning_failure'
        else:
            return 'cloning_failure'
    counter = 0
    # targets = ['build007']
    # for target in session.query(Build).filter(Build.status == 'sequencing').filter(Build.build_name.in_(targets)).order_by(Build.id):
    for target in session.query(Build).filter(Build.status == 'sequencing').order_by(Build.id):

        print('================= {} ================='.format(target.build_name))
        SEQFILE_PATH = "{}/{}/{}_seq_files".format(BUILDS_PATH,target.build_name,target.build_name)
        seq_plate = [plate for plate in target.plates if plate.plate_type == 'seq_plate'][0]
        new_wells = []
        for well in seq_plate.wells:
            print(counter)
            if well.seq_outcome != '':
                print(well.seq_outcome)
                continue
            id_num = well.parts.part_id
            print(id_num)
            if glob.glob("{}/*{}*{}*.ab1".format(SEQFILE_PATH,well.parts.part_id,forward_primer)):
                forfile = glob.glob("{}/*{}*{}*.ab1".format(SEQFILE_PATH,well.parts.part_id,forward_primer))[0]
                revfile = glob.glob("{}/*{}*{}*.ab1".format(SEQFILE_PATH,well.parts.part_id,reverse_primer))[0]
            elif glob.glob("{}/*{}*{}*.ab1".format(SEQFILE_PATH,dictionary[well.parts.part_id],forward_primer)):
                forfile = glob.glob("{}/*{}*{}*.ab1".format(SEQFILE_PATH,dictionary[well.parts.part_id],forward_primer))[0]
                revfile = glob.glob("{}/*{}*{}*.ab1".format(SEQFILE_PATH,dictionary[well.parts.part_id],reverse_primer))[0]
            else:
                print("Can't find seq_files",dictionary[well.parts.part_id])
                no_seq.append(id_num)
                continue

            for_seq_file = forfile.split("/")[-1]
            well.for_read = for_seq_file
            print("Forward: ",for_seq_file)
            rev_seq_file = revfile.split("/")[-1]
            well.rev_read = rev_seq_file
            print("Reverse: ",rev_seq_file)

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

            gene_seq = Seq.Seq(well.parts.seq)

            outcome, strange = verify_sequence(id_num,forward,reverse,gene_seq,backbone_seq,strange)

            if outcome == 'Unknown_sequence':
                outcome = align_unknown(id_num,forward,reverse,well,seq_plate)
            elif " - " in outcome:
                check.append(id_num)
                check.append(outcome)

            outcomes.append(outcome)
            well.seq_outcome = outcome
            well.parts.status = outcome

            print("Status: ",well.parts.status)
            print()
            counter += 1
        seq_plate.wells += new_wells

    for part in session.query(Part):
        part.eval_status()

    print("nan:")
    print(nan)
    print("small:")
    print(small)
    print("to be checked:")
    print(check)
    print("no seq files:")
    print(no_seq)
    print("no hit:")
    print(no_hit)
    print("strange:")
    print(strange)

    outcomes = pd.Series(outcomes)
    print(outcomes.value_counts())

    commit = int(input("Commit changes (1-yes, 2-no): "))
    if commit == 1:
        session.commit()
    return

if __name__ == "__main__":
    verify_seq()
















#
