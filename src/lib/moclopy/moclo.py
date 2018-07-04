import numpy
import time
import pandas as pd
from synbiolib import codon
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqUtils import MeltingTemp as mt
from Bio import SeqIO
from io import StringIO
from Bio.Alphabet import IUPAC
import unittest
import collections
from numpy.random import choice,randint
import logging

## =======
## Logging
## =======

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logging.basicConfig(filename='moclopy.log',level=logging.DEBUG)


## ==========
## Parameters
## ==========

### Fragger parameters
synthesis_max = 1500
enzymes = {
        "AarI" : {
            "seq" : "CACCTGC",
            "jump" : 4,
            "overhang" : 4
            },
        "BbsI" : {
            "seq" : "GAAGAC",
            "jump" : 2,
            "overhang" : 4
            },
        "BfuAI" : { 
            "seq" : "ACCTGC",
            "jump" : 4,
            "overhang" : 4
            },
        "BsaI" : {
            "seq" : "GGTCTC",
            "jump" : 1,
            "overhang" : 4
            },
        "BsmBI" : {
            "seq" : "CGTCTC",
            "jump" : 1,
            "overhang" : 4
            },
        "BtgZI" : { 
            "seq" : "GCGATG",
            "jump" : 10,
            "overhang" : 4
            },
        "SapI" : { 
            "seq" : "GCTCTTC",
            "jump" : 1,
            "overhang" : 3
            }
        }

FG_part_types = {
        "cds" : {
            "prefix" : "GGTCTCNA",
            "suffix" : "AGAGCTTNGAGACC"
            },
        "eukaryotic_promoter" : {
            "prefix" : "GGTCTCNGGAG",
            "suffix" : "AATGNGAGACC"
            },
        "prokaryotic_promoter" : {
            "prefix" : "GGTCTCNGGAG",
            "suffix" : "TACTNGAGACC"
            },
        "rbs" : { 
            "prefix" : "GGTCTCNTACT",
            "suffix" : "AATGNGAGACC"
            },
        "terminator" : {
            "prefix" : "GGTCTCNGCTT",
            "suffix" : "CGCTNGAGACC"
            },
        "operon" : { 
            "prefix" : "GGTCTCNGGAG",
            "suffix" : "CGCTNGAGACC"
            },
        "cds_aari" : {
            "prefix" : "CACCTGCNNNNGGAG",
            "suffix" : "CGCTNNNNGCAGGTG"
            },
        "cds_operon" : {
            "prefix" : "GGTCTCNA",
            "suffix" : "AGAGCTTNGAGACC"
            },
        "rp_selection" : {
            "prefix": "CACCTGCNNNNGGAG",
            "suffix": "CGCTNNNNGCAGGTG"
            }
        }

### Fixer parameters
default_table = codon.load_codon_table(taxonomy_id="custom_1", custom=True)

FG_end_codons = {
    'M': 'ATG',
    'W': 'TGG',
    'F': 'TTT',
    'L': 'CTG',
    'I': 'ATT',
    'V': 'GTG',
    'S': 'TCC',
    'P': 'CCA',
    'T': 'ACC',
    'A': 'GCC',
    'Y': 'TAC',
    'H': 'CAT',
    'Q': 'CAG',
    'N': 'AAC',
    'K': 'AAG',
    'D': 'GAT',
    'E': 'GAG',
    'C': 'TGC',
    'R': 'CGC',
    'G': 'GGC',
    '*': 'TGA'
}

cut_sites= {
        "BbsI": "GAAGAC",
        "BtgZI": "GCGATG",
        "BsaI": "GGTCTC",
        "BsmBI": "CGTCTC",
        "AarI": "CACCTGC",
        "BfuAI": "ACCTGC",
        "SapI": "GAAGAGC"
        }

## /\/\/\/\/\/\
## ============
## Fix sequence
## ============
## /\/\/\/\/\/\

def reverse_complement(seq):
    '''Creates reverse complement of a DNA string.'''
    seq = seq.upper()
    return seq.translate(str.maketrans("ATGCN","TACGN"))[::-1]

def transl_checker(DNA, protein):
    '''Checks translation of a protein based on the standard table.'''
    if Seq(DNA, IUPAC.unambiguous_dna).translate() == protein:
        return True
    else:
        return False

def translate(seq):
    '''Translates a string to the standard table'''
    return str(Seq(seq, generic_dna).translate(table=11))

def random_dna_sequence(length):
    '''Generates a random DNA sequence of a given length.'''
    return ''.join(np.random.choice(('A', 'C', 'T', 'G')) for _ in range(length))

def enzyme_finder(seq, cut_sites=cut_sites):
    '''Finds if a restriction enzyme is in a sequence and returns the enzyme name. If no enzyme is found, returns False'''
    seq = str(seq).upper()
    for enzyme,sequence in cut_sites.items():
        if sequence in seq:
            location = seq.find(cut_sites[enzyme])
            logger.debug("Found {} at location {}".format(enzyme,location))
            return location, len(cut_sites[enzyme])
        elif reverse_complement(sequence) in seq:
            location = seq.find(reverse_complement(cut_sites[enzyme]))
            logger.debug("Found {} at location {}".format(enzyme,location))
            return location, len(cut_sites[enzyme])
    return False

def homopolymer_finder(seq, max_homopolymer=7):
    '''Finds largest homopolymer in a sequence. If homopolymer is over threshold, returns position and length of homopolymer. Else, returns False'''
    seq = seq.upper()
    prev = ''
    count = 0
    max_count = 0
    loc = None
    for i, char in enumerate(seq):
        if char == prev:
            count += 1
        else:
            count = 0
        if count > max_count:
            max_count = count
            loc = i - count
        prev = char
    if max_count > max_homopolymer:
        logger.debug("Found {}bp homopolymer at location {}".format(max_count+1,loc))
        return loc, max_count + 1
    else:
        return False

def repeat_finder(seq):
    '''Finds largest repeat in a sequence and returns the position of the repeat.'''
    string = seq.upper()
    l = list(string)
    d = collections.deque(string[1:])
    match = []
    longest_match = []
    while d:
        for i, item in enumerate(d):
            if l[i]==item:
                match.append(item)
            else:
                if len(longest_match) < len(match):
                    longest_match = match
                match = []
        d.popleft()
    repeat_sequence = ''.join(longest_match)
    if len(repeat_sequence) > 20:
        location = string.find(repeat_sequence)
        length = len(repeat_sequence)
        logger.debug("Found long {}bp repeat at location {}".format(length,location))
        return location,length
    if len(repeat_sequence) < 20 and int(mt.Tm_Wallace(Seq(repeat_sequence))) > 55:
        location = string.find(repeat_sequence)
        length = len(repeat_sequence)
        logger.debug("Found long {}bp repeat at location {}".format(length,location))
        return location,length
    else:
        return False

def FG_MoClo_ends(seq):
    '''Changes end of any protein sequence to be FG MoClo compatible.'''
    return str(seq[:-6] + ''.join(list(map(lambda x: FG_end_codons[str(Seq(x, IUPAC.unambiguous_dna).translate())], [seq[-6:-3],seq[-3:]]))))

def gc_range_finder(seq):
    '''Finds highest and lowest 50bp GC window in a sequence. Returns sequence ot change position'''
    seq = seq.upper()
    highest_gc = (0, False)
    lowest_gc = (1, False)
    average = gc_content(seq)
    for start_base in range(len(seq)-50):
        gc_seq = seq[start_base:start_base+50]
        gc_percent = gc_content(gc_seq)
        if gc_percent > highest_gc[0]:
            highest_gc = (gc_percent, gc_seq)
        if gc_percent < lowest_gc[0]:
            lowest_gc = (gc_percent, gc_seq)
    if highest_gc[0] - lowest_gc[0] > .51: # Finds if a change needs to be made
        if abs(highest_gc[0] - average) > abs(lowest_gc[0] - average): # Finds whether it should be lowest or highest gc region
            return seq.find(highest_gc[1]), len(highest_gc[1]), "AT"
        else:
            return seq.find(lowest_gc[1]), len(lowest_gc[1]), "GC"
    else:
        return False

def gc_content(seq):
    return (seq.count("G") + seq.count("C")) / len(seq)
def gc_in_range(seq, min_gc=0.3, max_gc=0.65 ):
    return (gc_content(seq) > min_gc) & (gc_content(seq) <= max_gc)

## ==========
## DNA checks
## ==========

def cds_check(gene_id,seq):
    def is_seq(seq):
        return seq.replace(r'\s +', '').strip() != ""

    def is_triplet(seq):
        return len(seq) % 3 == 0

    def translate(seq):
        return str(Seq(seq, generic_dna).translate(table=11))

    def dna_matches_aa(seq, aa):
        prot = translate(seq)
        if prot[-1] == "*":
            prot = prot[:-1]
        return prot == aa

    def has_start(seq):
        return seq[:3] == "ATG"

    def has_stop(seq):
        return translate(seq)[-1] == "*"

    def has_no_internal_stops(seq):
        return not "*" in translate(seq)[:-1]

    basic_tests = [is_seq, is_triplet, has_start, has_stop, has_no_internal_stops]
    for test in basic_tests:
        if not test(seq):
            raise ValueError("{} failed on {}".format(gene_id, test.__name__))
    return "Clear"

def input_checker(gene_id,seq):
    seq = seq.upper()
    if "A"*11 in seq:
        raise ValueError('{} failed on long "A" homopolymer'.format(gene_id))
    elif "T"*11 in seq:
        raise ValueError('{} failed on long "T" homopolymer'.format(gene_id))
    elif "G"*11 in seq:
        raise ValueError('{} failed on long "G" homopolymer'.format(gene_id))
    elif "C"*11 in seq:
        raise ValueError('{} failed on long "C" homopolymer'.format(gene_id))
    elif "GGTCTC" in seq or "GAGACC" in seq:
        print(seq)
        raise ValueError('{} failed on has_bsaI'.format(gene_id))
    elif "GAAGAC" in seq or "GTCTTC" in seq:
        raise ValueError('{} failed on has_bbsI'.format(gene_id))
    elif len(seq) < 10:
        raise ValueError('{} failed on small sequence (<10bp)'.format(gene_id))
    elif len(seq) > 20000:
        raise ValueError('{} failed on large sequence (>20kbp)'.format(gene_id))
    else:
        return "Clear"

## ===============
## Recode sequence
## ===============

def generate_codon_dict(seq):
    '''Creates a dictionary to store the codons that have been attempted thus far'''
    coding_dna = Seq(seq, generic_dna)
    translation = coding_dna.translate()
    return dict(enumerate(zip(translation,[[seq[n:n+3]] for n in range(0, len(seq), 3)])))


def replace_codons(table,seq,change_list,pos_list,makelist=True,halfway=True,GC_bias="None"):
    seq = seq.upper()
    if len(pos_list) == 0:
        logger.debug("len(pos_list) = 0. Reseting pos_list")
        return seq, generate_codon_dict(seq)
    
    # If taking standard input from x_finder functions, create initial position list
    if makelist:
        position = pos_list[0] - (pos_list[0] % 3)
        pos_list = list(map(lambda x: int(x/3), list(range(position, position + (pos_list[1] // 3) * 3, 3))))

    # Choose position to change from middle of list, rather than edge
    if halfway:
        position = pos_list[int(len(pos_list)/2)]
        if position == 0: # Handles single positions
            position = pos_list[0]
    else:
        position = pos_list[0]
    
    # Setup pandas table
    table = table.reset_index()
    section = table[table.AA == change_list[position][0]]
    
    # Allow GC bias in changes
    if GC_bias == "GC":
        bad_rows = section[section.Triplet.str.endswith('A')].Triplet
        bad_rows = bad_rows.append(section[section.Triplet.str.endswith('T')].Triplet)
        section = section[-section.Triplet.isin(bad_rows)]
    elif GC_bias == "AT":
        bad_rows = section[section.Triplet.str.endswith('G')].Triplet
        bad_rows = bad_rows.append(section[section.Triplet.str.endswith('C')].Triplet)
        section = section[-section.Triplet.isin(bad_rows)]
    
    # Checks for available codons. If none, skip to next codon
    available = section[~section.Triplet.isin(change_list[position][1])]
    available = available[available['Fraction'] != 0.000000]
    if len(available) == 0:
        pos_list.remove(position)
        return replace_codons(table,seq,change_list,pos_list,makelist=False,GC_bias=GC_bias)
    
    # Randomly choose new codon based on frequencies
    options = available.Fraction.tolist()
    norm = [float(i)/sum(options) for i in options]
    codons = available.Triplet.tolist()
    new_codon = codons[int(choice(len(options),1,p=norm))]
    change_list[position][1].append(new_codon)
    logger.debug("{} -> {} at position {}\n".format(seq[position*3:(position*3)+3],new_codon,position*3))
    seq = ''.join((seq[:position*3],new_codon,seq[(position*3)+3:]))
    return seq, change_list

## ============================
## Recode list and fix sequence
## ============================

def recode_list(table,gene_id,seq,change_list={},fix_attempts=0):
    fix_attempts = fix_attempts + 1
    max_attempts = 500
    seq = FG_MoClo_ends(seq)
    if fix_attempts > max_attempts:
        logger.debug("{} could not be fixed. Please manually check this sequence.".format(gene_id))
        logger.debug("Terminating program")
        raise ValueError('Max attempts reached')
    # Check cutsites
    cutsites = enzyme_finder(seq)
    if cutsites:
        logger.debug("Fixing cutsites in {}".format(gene_id))
        seq, change_list = replace_codons(table,seq,change_list,cutsites)
        return recode_list(table, gene_id, seq, change_list, fix_attempts=fix_attempts)
    repeats = repeat_finder(seq)
    if repeats:
        logger.debug("Fixing repeat regions in {}".format(gene_id))
        seq, change_list = replace_codons(table,seq,change_list,repeats)
        return recode_list(table, gene_id, seq, change_list, fix_attempts=fix_attempts)
    homopolymers = homopolymer_finder(seq)
    if homopolymers:
        logger.debug("Fixing homopolymer regions in {}".format(gene_id))
        seq, change_list = replace_codons(table,seq,change_list,homopolymers)
        return recode_list(table, gene_id, seq, change_list, fix_attempts=fix_attempts)
    gc_range = gc_range_finder(seq)
    if gc_range:
        logger.debug("Fixing gc_range regions in {}".format(gene_id))
        gc_range_position = gc_range[0], gc_range[1]
        gc_bias = gc_range[2]
        seq, change_list = replace_codons(table,seq,change_list,gc_range_position,GC_bias=gc_bias)
        return recode_list(table, gene_id, seq, change_list, fix_attempts=fix_attempts)
    return seq 

def fix_sequence(gene_id,seq,table=default_table):
    print("Fixing {}".format(gene_id))
    try:
        fixed_gene = recode_list(table,gene_id,seq,generate_codon_dict(seq))
    except Exception as e:
        print("{} FAILED ON {}".format(gene_id,e))
        fixed_gene = "FAILED"
    return fixed_gene

## /\/\/\/\/\/\/\/\/
## =================
## Fragment sequence
## =================
## /\/\/\/\/\/\/\/\/
def error_finder(seq,overhangs=[]):
    error_list=[]
    seq = seq.upper()
    homopolymer = homopolymer_finder(seq)
    if homopolymer:
        print("homopolymer")
        middle = int(homopolymer[1]/2) + homopolymer[0]
        error_list.append(middle)
    repeat = repeat_finder(seq)
    if repeat:
        print("repeat")
        end_repeat = repeat[0]+repeat[1]
        error_list.append(end_repeat)
    gc_range = gc_range_finder(seq)
    #if gc_range:
    #    print("gcrange")
    #    middle = int(gc_range[1]/2) + gc_range[0]
    #    error_list.append(middle)
    error_list = sorted(error_list)
    if error_list == []:
        return False
    else:
        return error_list[0]

def part_type_preparer(part_type, seq, suffix="", prefix=""):
    if part_type == "vector":
        seq = seq + seq[-4:]
        return seq
    if part_type == "custom":
        seq = prefix + seq + suffix
        return seq
    else:
        part = FG_part_types[part_type]
        if len(seq) < 300:
            N_replace = "G"
        else:
            N_replace = "A"
        seq = part["prefix"].replace("N", N_replace) + seq + part["suffix"].replace("N", N_replace)
        seq = "GGAG" + seq + "CGCT"
        return seq

def position_fragmenter(seq, pos_list, overhangs=[]):
    seq = seq.upper()
    frags = []
    for position in pos_list:
        overhang = seq[position-2:position+2]
        if overhang in overhangs:
            return position_fragmenter(seq, [position+1], overhangs)
        frags += [seq[0:position+2], seq[position-2:]]
    return frags, overhang

def fragger(seq,overhangs=["GGAG","CGCT"],synthesis_max=1500):
    error_checked = []
    frags = []
    error = error_finder(seq)
    if error == False:
        error_checked.append(seq)
    while error != False:
        frag, overhang = position_fragmenter(seq,[error],overhangs)
        error_checked += frag
        overhangs.append(overhang)
        seq = frag[1]
        error = error_finder(seq)
    for seq in error_checked:
        if len(seq) > synthesis_max:
            num_frags = len(seq) // synthesis_max + 1
            frag_len = len(seq) // num_frags
            pos_list = []
            for fragment in range(num_frags)[1:]:
                pos_list.append(fragment*frag_len)
            frag, overhang = position_fragmenter(seq,pos_list,overhangs)
            frags += frag
            overhangs.append(overhang)
        else:
            frags.append(seq)
    return frags

def fragment_sequence(gene_id,seq,part_type,cloning_enzyme_prefix="GAAGACTT",cloning_enzyme_suffix="GCGTCTTC",synthesis_max=1500):
    seq = seq.upper()
    seq = part_type_preparer(part_type,seq)
    fragments = fragger(seq,synthesis_max=synthesis_max)
    fragments = list(map(lambda x: cloning_enzyme_prefix + x + cloning_enzyme_suffix, fragments))
    return fragments

## /\/\/\/\/\/\/\
## ==============
## Sequence_input
## ==============
## /\/\/\/\/\/\/\

def sequence_input(gene_id,seq,part_type):
    success = True
    seq = seq.upper()
    part_type = part_type.lower() 
    # Check cds for proper format. If it's alright, fix it 
    if part_type == 'cds':
        try:
            cds_checker(gene_id,seq)
        except Exception as e:
            print("{} failed on {}".format(gene_id,e))
            return "FAILED"
        seq = fix_sequence(gene_id,seq)

    # Check all seqs for bad elements
    try:
        input_checker(gene_id,seq)
    except Exception as e:
        print('{} failed on {}'.format(gene_id,e))
        return "FAILED"

    # Fragment sequence
    frags = fragment_sequence(gene_id,seq,part_type)
    print("{} successfully checked and ready to submit".format(gene_id))
    return frags



