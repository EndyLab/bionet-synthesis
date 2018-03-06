def fragment_gene(sequence,entry_type):
    import numpy as np
    import pandas as pd

    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna
    from difflib import ndiff
    import re

    import json
    import io
    import sys
    import math

    ## ==========================================================
    ## Configurations
    ## ==========================================================
    synthsis_config = {
        "max_length" : 1500,
        "min_length" : 300
    }
    pcr_config = {
        "max_length" : 5000
    }
    standard_flanks = {
        "prefix" : "GGAG",
        "suffix" : "CGCT"
    }
    enzyme_configuration = {
        "AarI" : {
            "prefix" : "N",
            "suffix" : "N"
            },
        "BbsI" : {
             "prefix" : "N",
             "suffix" : "N"
             },
        "BfuAI" : {
            "prefix" : "N",
            "suffix" : "N"
            },
        "BsaI" : {
            "prefix" : "N",
            "suffix" : "N"
            },
        "BsmBI" : {
            "prefix" : "N",
            "suffix" : "N"
            },
        "BtgZI" : {
            "prefix" : "N",
            "suffix" : "N"
            },
        "SapI" : {
            "prefix" : "N",
            "suffix" : "N"
            },
    }
    # Establish part types
    part_type = dict(
        eukaryotic_promoter = {
            # FreeGenes MoClo Eukaryotic Promoter definition
            "cloning_enzyme" : "BbsI",
            "small_cloning_enzyme" : "BtgZI",
            "retrieval_enzyme" : "BsaI",
            "prefix" : "GGAG",
            "suffix" : "AATG"
        },
        prokaryotic_promoter = {
            # FreeGenes MoClo Prokaryotic Promoter definition
            "cloning_enzyme" : "BbsI",
            "small_cloning_enzyme" : "BtgZI",
            "retrieval_enzyme" : "BsaI",
            "prefix" : "GGAG",
            "suffix" : "TACT"
        },
        rbs = {
            # FreeGenes MoClo RBS definition
            "cloning_enzyme" : "BbsI",
            "small_cloning_enzyme" : "BtgZI",
            "retrieval_enzyme" : "BsaI",
            "prefix" : "TACT",
            "suffix" : "AATG"
        },
        terminator = {
            # FreeGenes MoClo Prokaryotic Promoter definition
            "cloning_enzyme" : "BbsI",
            "small_cloning_enzyme" : "BtgZI",
            "retrieval_enzyme" : "BsaI",
            "prefix" : "GCTT",
            "suffix" : "CGCT"
        }
    )
    # Config fragments
    BsaI_prefix = "GGAGGTCTCA".lower()
    BsaI_suffix = "CGAGACCGCT".lower()

    ## ==========================================================
    ## Add Prefix and Suffix
    ## ==========================================================
    if entry_type not in part_type.keys():
        print("not a valid type")
        return
    part = part_type[entry_type]
    full_sequence = BsaI_prefix + part["prefix"] + sequence + part["suffix"] + BsaI_suffix
    print(full_sequence)
    print("length",len(full_sequence))


## Test the function
import pandas as pd

data = pd.read_csv("./gene.csv")
for index,row in data.iterrows():
    fragment_gene(row["Sequence"],row["Part_Type"])












#
