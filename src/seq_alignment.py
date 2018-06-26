import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

import sys
import os
import re
import math
import glob
import json
import numpy as np
import pandas as pd
from datetime import datetime
from matplotlib import pyplot as plt

from Bio import pairwise2
from Bio import Seq
from Bio import SeqIO
from Bio import SeqFeature
from Bio.Blast import NCBIXML
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist
import shutil

from config import *
from db_config import *
import ot_functions as ot

session,engine = connect_db()

## Query the database for the parts to be sequenced and also pulls in their fragment sequences
builds = []
print('Builds awaiting a manual check:')
for i,build in enumerate(session.query(Build).filter(Build.status == 'sequencing').order_by(Build.build_name)):
    print('{}: {}'.format(i,build.build_name))
    builds.append(build)

ans = ot.request_info('Which build would you like to check: ',type='int',select_from=[x for x in range(len(builds))])
tar_build = builds[ans]

target = tar_build.build_name

length = shutil.get_terminal_size().columns


fasta_path = "{}/builds/{}/{}_fasta_files".format(BASE_PATH,target,target)

query_wells = "SELECT parts.part_id,parts.cloning_enzyme,parts.seq,parts.part_name,fragments.fragment_name,fragments.seq as frag_seq,wells.address,wells.plate_type,wells.vector,builds.build_name,wells.misplaced FROM parts \
        INNER JOIN wells ON parts.id = wells.part_id\
        INNER JOIN plates ON wells.plate_id = plates.id\
        INNER JOIN builds ON plates.build_id = builds.id\
        INNER JOIN part_frag ON parts.id = part_frag.part_id\
        INNER JOIN fragments ON part_frag.fragment_id = fragments.id\
        WHERE builds.build_name = '{}'\
            AND plates.plate_type = 'seq_plate'\
            AND wells.misplaced != 'True'".format(target)

df = pd.read_sql_query(query_wells,con=engine)


enzyme_dict = {
    'BbsI' : {
        'seq' : 'GAAGAC',
        'offset' : 2
    },
    'BtgZI' : {
        'seq' : 'GCGATG',
        'offset' : 10
    },
}

def build_vector(row):
    '''
    Imports and parses a genbank file describing a vector. It assumes that there is
    a forward and reverse sequencing primers annotated on the correct strands with the
    annotation type of 'primer_bind'
    '''
    vector = row['vector']
    enzyme = row['cloning_enzyme']
    seq = row['full_seq']
    record = SeqIO.read(glob.glob("{}/raw_files/vectors/*{}*".format(BASE_PATH,vector))[0], "genbank")
    for f in record.features:
        if f.type == 'primer_bind' and f.strand == -1:
            reverse = reverse_complement(str(f.extract(record.seq)))
        elif f.type == 'primer_bind' and f.strand == 1:
            forward = f.extract(record.seq)
    vector_seq = record.seq

    for_ind = str(vector_seq).find(str(forward))
    rev_ind = str(vector_seq).find(reverse) + len(reverse)

    seq_region = str(vector_seq)[for_ind:rev_ind]

    beginning = seq_region[:seq_region.find(reverse_complement(enzyme_dict[enzyme]['seq'])) - enzyme_dict[enzyme]['offset']]
    ending = seq_region[seq_region.find(enzyme_dict[enzyme]['seq']) + enzyme_dict[enzyme]['offset'] + len(enzyme)+6:]
    return beginning+seq+ending,str(forward),str(reverse)

def loadsequencing(file, threshold=0.9):
    '''Load sequencing reads and trim them based on the specified threshold'''
    seq = SeqIO.read(file, 'abi')

    maxq = np.max(seq.letter_annotations['phred_quality'])
    rolling = pd.Series(seq.letter_annotations['phred_quality']).rolling(window=20).mean()
    start = (rolling > maxq * threshold).idxmax()
    end = (rolling > maxq * threshold)[::-1].idxmax()

    return seq.seq[start:end], seq, start, end, np.mean(seq.letter_annotations['phred_quality'][start:end])

def reverse_complement(seq):
    return seq.translate(str.maketrans("ATGC","TACG"))[::-1]

def digest(seq,enzyme):
    '''
    Cuts the sequence at the specified restriction site. It makes several assupmtions:
    The cutsites are pointing inward toward the center of the sequence and it
    eliminates the overhang on the 5' end but retains the 'sticky end' on the 3'
    '''
    return seq[seq.find(enzyme_dict[enzyme]['seq'])+len(enzyme_dict[enzyme]['seq'])+enzyme_dict[enzyme]['offset']+4:seq.find(reverse_complement(enzyme_dict[enzyme]['seq']))-enzyme_dict[enzyme]['offset']]

def find_reads(row):
    forfile = glob.glob("{}/builds/{}/{}_seq_files/*{}*{}*.ab1".format(BASE_PATH,target,target,row['part_id'],'For'))[0]
    revfile = glob.glob("{}/builds/{}/{}_seq_files/*{}*{}*.ab1".format(BASE_PATH,target,target,row['part_id'],'Rev'))[0]
    return forfile,revfile

def import_reads(row):
    forward_untrim, _, _, _, forward_qual = loadsequencing(row['forfile'])
    reverse_untrim, revseq, _, _, reverse_qual = loadsequencing(row['revfile'])

    if len(forward_untrim) < 50 or len(reverse_untrim) < 50:
        forward_untrim, _, _, _, forward_qual = loadsequencing(row['forfile'], threshold=0.8)
        reverse_untrim, revseq, _, _, reverse_qual = loadsequencing(row['revfile'], threshold=0.8)
        small = True
    else:
        small = False
    return str(forward_untrim),str(reverse_untrim),forward_qual,reverse_qual,small

def trim_reads(row):
    forward = row['forward'][:row['forward'].find(reverse_complement(row['rev_primer']))] # Start at start codon
    reverse = reverse_complement(row['reverse'][:row['reverse'].find(reverse_complement(row['for_primer']))]) # Stop at stop
    return forward,reverse

def create_alignment_str(form):
    ori_align,res,tar_align,_,_ = form.split("\n")
    ori_seg = [ori_align[n:n+length] for n in range(0,len(ori_align),length)]
    res_seg = [res[n:n+length] for n in range(0,len(res),length)]
    tar_seg = [tar_align[n:n+length] for n in range(0,len(tar_align),length)]
    align = ''
    full_align = zip(ori_seg,res_seg,tar_seg)
    for ori,res,tar in full_align:
        align += ori+'\n'+res+'\n'+tar+'\n'+'\u001b[0m'
    return align

def align_reads(row):
    forward_align = pairwise2.align.globalms(row['full_seq'], row['forward'],1,0,-1,-1, one_alignment_only=True, penalize_end_gaps=False)[0]
    reverse_align = pairwise2.align.globalms(row['full_seq'], row['reverse'],1,0,-1,-1, one_alignment_only=True, penalize_end_gaps=False)[0]

    for_align_str = create_alignment_str(format_alignment(*forward_align))
    rev_align_str = create_alignment_str(format_alignment(*reverse_align))

    for_raw = len(row['forward'])
    rev_raw = len(row['reverse'])
    target_length = len(row['full_seq'][row['full_seq'].find(row['for_primer'])+len(row['for_primer']):row['full_seq'].find(reverse_complement(row['rev_primer']))])

    if for_raw <= target_length:
        for_score = for_raw - forward_align[2]
    else:
        for_score = target_length - forward_align[2]
    if rev_raw <= target_length:
        rev_score = rev_raw - reverse_align[2]
    else:
        rev_score = target_length - reverse_align[2]

    return for_raw,forward_align[2],for_score,rev_raw,reverse_align[2],rev_score,target_length,for_align_str,rev_align_str

def analyze_align(row):
    if row['for_raw'] < 100 and row['rev_raw'] < 100:
        outcome = "bad_reads"
    elif row['for_raw'] < 100 and row['rev_score'] > 10 or row['rev_raw'] < 100 and row['for_score'] > 10:
        outcome = "cloning_failure"
    elif row['for_score'] == 0:
        if row['rev_score'] == 0:
            outcome = "perfect"
        elif row['rev_score'] < 10:
            outcome = "mutation_{}-{}".format(row['for_score'],row['rev_score'])
        else:
            outcome = "bad_reverse"
    elif row['for_score'] < 10:
        if row['rev_score'] == 0:
            outcome = "mutation_{}-{}".format(row['for_score'],row['rev_score'])
        elif row['rev_score'] < 10:
            outcome = "mutation_{}-{}".format(row['for_score'],row['rev_score'])
        else:
            outcome = "bad_reverse"
    else:
        if row['rev_score'] == 0:
            outcome = "bad_forward"
        elif row['rev_score'] < 10:
            outcome = "bad_forward".format(row['for_score'],row['rev_score'])
        else:
            outcome = "bad_clone"
    return outcome

def export_fasta(row):
    with open('{}/{}.fasta'.format(fasta_path,row['part_id']),'w') as fasta:
        fasta.write('>{}\n{}'.format(row['part_id'],row['full_seq']))

def check_alignment(row):
    if row['outcome'] == 'perfect':
        return row['outcome']
    ot.print_center('================== \u001b[1mForward\u001b[0m ==================\n')
    forward = row['for_align_str'].replace('.','\u001b[41m.\u001b[0m')
    forward = re.sub(r'(?<=[A-Z])-(?=[A-Z])', '\u001b[41m-\u001b[0m',forward)
    print(forward)
    ot.print_center('================== \u001b[1mReverse\u001b[0m ==================\n')
    reverse = row['rev_align_str'].replace('.','\u001b[41m.\u001b[0m')
    reverse = re.sub(r'(?<=[A-Z])-(?=[A-Z])', '\u001b[41m-\u001b[0m',reverse)
    print(reverse)
    ot.print_center('{}: \u001b[31;1m\u001b[1m{}\u001b[0m'.format(row['part_id'],row['outcome']))
    print('Options: k:keep, c:check, p:perfect, m:mutation, b:bad_reads,f:cloning_failure')
    ans = ot.request_info('Please state the outcome: ',type='string',select_from=['c','k','p','m','b','f'])
    if ans == 'k':
        return row['outcome']
    elif ans == 'c':
        return 'CHECK'
    elif ans == 'p':
        return 'perfect'
    elif ans == 'm':
        return 'mutation'
    elif ans == 'b':
        return 'bad_reads'
    elif ans == 'f':
        return 'cloning_failure'

print('Starting analysis:',datetime.now())

df['frag_seq'] = df.apply(lambda row: digest(row['frag_seq'],row['cloning_enzyme']),axis=1)

grouped = df.groupby('part_id')
full_seq = grouped.agg({'frag_seq': lambda x: "%s" % ''.join(x)})

df = df[['part_id','seq','address','build_name','vector','cloning_enzyme']]
df = df.set_index('part_id')
df = df.drop_duplicates()
df = pd.concat([df,full_seq],axis=1)
df = df.reset_index()
df = df.rename(columns={'index': 'part_id','frag_seq': 'full_seq'})

df['full_seq'],df['for_primer'],df['rev_primer'] = zip(*df.apply(build_vector,axis=1))
df['correct'] = df.apply(lambda row: True if row['seq'] in row['full_seq'] else False,axis=1)
df['forfile'],df['revfile'] = zip(*df.apply(find_reads,axis=1))
df['forward'],df['reverse'],df['forward_qual'],df['reverse_qual'],df['small'] = zip(*df.apply(import_reads,axis=1))
df['forward'],df['reverse'] = zip(*df.apply(trim_reads,axis=1))
df['for_raw'],df['for_align'],df['for_score'],df['rev_raw'],df['rev_align'],df['rev_score'],df['seq_length'],df['for_align_str'],df['rev_align_str'] = zip(*df.apply(align_reads,axis=1))
df['outcome'] = df.apply(analyze_align,axis=1)

## ===============================
## Finding unknowns
## ===============================

BLAST_db_path = '{}/raw_files/BLAST_db/current_BLAST_db.fsa'.format(BASE_PATH)

def create_BLAST_database(file_path='{}/raw_files/BLAST_db/current_BLAST_db.fsa'.format(BASE_PATH)):
    '''Generates a multi-input fasta file to store all of the sequences in the
    in the BLAST database'''
    ot.print_center('...Generating BLAST database...')
    db_counter = 1
    with open(file_path,"w+") as fsa:
        for part in session.query(Part).order_by(Part.part_id):
            fsa.write(">{}|BLAST_db|{}\n{}\n".format(part.part_id,db_counter,part.seq))
            db_counter += 1
    # Convert the FSA file into a BLAST database
    os.system("makeblastdb -in {} -parse_seqids -dbtype nucl".format(file_path))

create_BLAST_database()

def blast_seq(name,sequence, E_VALUE_THRESH=0.04):
    '''Runs a BLAST search with a read against the database and returns the ID of
    the part that is identified'''
    # Create fasta files to run the BLAST search
    query = BASE_PATH + "/raw_files/BLAST_db/query.fasta"
    with open(query,"w+") as seq_file:
        seq_file.write(">{}\n{}".format(name,sequence))

    # Run the BLAST search
    output = BASE_PATH + "/raw_files/BLAST_db/results.xml"
    os.system("blastn -query {} -db {}  -out {} -evalue 0.001 -outfmt 5".format(query,BLAST_db_path,output))

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

def find_unknown(row):
    '''Uses both reads to search the database for any relevant matches'''
    if row['outcome'] != 'bad_clone':
        return None
    for_hit = blast_seq(row['part_id'],row['forward'])
    rev_hit = blast_seq(row['part_id'],row['reverse'])
    if for_hit == rev_hit:
        return for_hit
    else:
        return 'unknown'

df['unknown'] = df.apply(find_unknown,axis=1)
df_unknown = df[df.outcome == 'bad_clone']
df_unknown = df_unknown[['unknown','address','vector','cloning_enzyme','build_name','for_primer','rev_primer','forfile','revfile']]
df_unknown = df_unknown.rename(columns={'unknown':'part_id'})

def gen_full_unknown_seq(row):
    '''Digests and assembles the complete sequence from all fragments'''
    full_seq = ''
    for frag in session.query(Fragment).join(Part,Fragment.parts).filter(Part.part_id == row['part_id']):
        frag_seq = digest(frag.seq,row['cloning_enzyme'])
        full_seq += frag_seq
    return full_seq

df_unknown['full_seq'] = df_unknown.apply(gen_full_unknown_seq,axis=1)
df_unknown['full_seq'],df_unknown['for_primer'],df_unknown['rev_primer'] = zip(*df_unknown.apply(build_vector,axis=1))
df_unknown['forward'],df_unknown['reverse'],df_unknown['forward_qual'],df_unknown['reverse_qual'],df_unknown['small'] = zip(*df_unknown.apply(import_reads,axis=1))
df_unknown['forward'],df_unknown['reverse'] = zip(*df_unknown.apply(trim_reads,axis=1))
df_unknown['for_raw'],df_unknown['for_align'],df_unknown['for_score'],df_unknown['rev_raw'],df_unknown['rev_align'],df_unknown['rev_score'],df_unknown['seq_length'],df_unknown['for_align_str'],df_unknown['rev_align_str'] = zip(*df_unknown.apply(align_reads,axis=1))
df_unknown['outcome'] = df_unknown.apply(analyze_align,axis=1)



fasta_path = "{}/builds/{}/{}_fasta_files".format(BASE_PATH,target,target)
ot.make_directory(fasta_path)

df.apply(export_fasta,axis=1)

print('Completed analysis: ',datetime.now())
print()
print(df.outcome.value_counts())
print()
ot.print_center('Press Enter to review questionable sequences')
input()

df['outcome'] = df.apply(check_alignment,axis=1)
df_unknown['outcome'] = df_unknown.apply(check_alignment,axis=1)

print('New outcome totals:')
print(df.outcome.value_counts())
print(df_unknown.outcome.value_counts())

ot.print_center('...Updating the database...')

def update_database(row):
    '''Adds all of the alignment information into the database'''
    tar_well = session.query(Well).join(Plate,Well.plates).join(Build,Plate.builds)\
        .join(Part,Well.parts).filter(Build.build_name == row['build_name'])\
        .filter(Part.part_id == row['part_id']).filter(Plate.plate_type == 'seq_plate').one()
    tar_well.for_read = row['forfile']
    tar_well.rev_read = row['revfile']
    tar_well.seq_outcome = row['outcome']
    tar_well.misplaced = row['unknown']
    return

df.apply(update_database,axis=1)

def add_unknowns(row):
    '''Generates a new well for the accidental parts and adds in the necessary information'''
    if row['outcome'] == 'bad_clone':
        return
    tar_part = session.query(Part).filter(Part.part_id == row['part_id']).one()
    tar_plate = session.query(Plate).join(Build,Plate.builds).filter(Build.build_name == row['build_name']).filter(Plate.plate_type == 'seq_plate').one()
    new_well = Well('seq_plate',tar_part,row['address'],for_read=row['forfile'],rev_read=row['revfile'],seq_outcome=row['outcome'])
    tar_plate.wells.append(new_well)
    return

df_unknown.apply(add_unknowns,axis=1)

full = pd.concat([df,df_unknown],ignore_index=True)
tar_build = session.query(Build).filter(Build.build_name == target).one()

if len(full[full.outcome == 'CHECK']) > 0:
    ot.print_center('...Generating csv to manually check...')
    full = full[full.outcome == 'CHECK']
    full['manual'] = ''
    full = full[['part_id','outcome','manual','address','forfile','revfile','forward_qual','reverse_qual','for_raw','for_align','for_score','rev_raw','rev_align','rev_score','unknown']]
    full.to_csv('{}/builds/{}/{}_manual_check.csv'.format(BASE_PATH,target,target),index=False)
    tar_build.status = 'manual_check_pending'
    ot.print_center('===== MANUAL CHECKS REQUIRED =====')
    print('Enter any of the following in the "manual" column: k:keep, p:perfect, m:mutation, f:cloning_failure')
    print('When finished run the `manual_align.py` script to add results into the database')
else:
    ot.print_center('No sequences to review, analysis is complete')
    tar_build.status = 'complete'

ot.print_center('...Committing results to the database...')
session.commit()


#
