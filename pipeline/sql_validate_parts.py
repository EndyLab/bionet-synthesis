import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from difflib import ndiff
import re

import io
import sys
import math

from config import *
import db_config
session,engine = db_config.connect_db()

def validate_part(part):
    print(part.part_id)
    frag_num = [[frag,frag.fragment_name[-1],frag.cloning_method] for frag in part.fragments]
    frag_num = sorted(frag_num, key=lambda tup: tup[1])
    # print(frag_num)
    if len(frag_num) == 1:
        print('single')
        parts = [gene for gene in frag_num[0][0].parts]
        if len(parts) > 1:
            print('linked')
            if frag_num[0][0].cloning_method == 'entry_2':
                if frag_num[0][0].seq[8:12] != 'AATG':
                    input("5' failed")
                if frag_num[0][0].seq[-20:-16] != 'TGAA':
                    input("3' failed")

            elif frag_num[0][0].cloning_method == 'entry_1':
                if frag_num[0][0].seq[8:12] != 'GGAG':
                    input("5' failed")
                if frag_num[0][0].seq[-20:-16] != 'CGCT':
                    input("3' failed")
            seq = frag_num[0][0].seq
            if seq[seq.find('GTCTTC')-6:seq.find('GTCTTC')-2] != 'TGAA':
                print("Found:",seq[seq.find('GTCTTC')-6:seq.find('GTCTTC')-2])
            else:
                print('bbsI good')
            if seq[seq.find('GCGATG')+16:seq.find('GCGATG')+20] != 'AATG':
                print("Found:",seq[seq.find('GCGATG')+16:seq.find('GCGATG')+20])
            else:
                print('btgZI good')


        elif frag_num[0][0].cloning_method == 'entry_2':
            if frag_num[0][0].seq[8:12] != 'AATG':
                input("5' failed")
            if frag_num[0][0].seq[-12:-8] != 'TGAA':
                input("3' failed")

        elif frag_num[0][0].cloning_method == 'entry_1':
            if frag_num[0][0].seq[8:12] != 'GGAG':
                input("5' failed")
            if frag_num[0][0].seq[-12:-8] != 'CGCT':
                input("3' failed")
        else:
            input("neither cloning method")
    else:
        print('multi')
        if frag_num[0][0].cloning_method == 'entry_2':
            if frag_num[0][0].seq[8:12] != 'AATG':
                input("5' failed")
            if frag_num[-1][0].seq[-12:-8] != 'TGAA':
                input("3' failed")

        elif frag_num[0][0].cloning_method == 'entry_1':
            if frag_num[0][0].seq[8:12] != 'GGAG':
                input("5' failed")
            if frag_num[-1][0].seq[-12:-8] != 'CGCT':
                input("3' failed")
        pairs = [[i,j] for i,j in zip(frag_num,frag_num[1:])]
        for pair in pairs:
            if pair[0][0].seq[-12:-8] != pair[1][0].seq[8:12]:
                input("Junction failed")


if __name__ == "__main__":
    counter = 0
    for part in session.query(Part).order_by(Part.part_id):
        if counter == 750:
            break
        counter += 1
        validate_part(part)












#
