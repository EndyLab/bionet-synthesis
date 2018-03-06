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
            "prefix" : "CACCTGCCCTA",
            "suffix" : "ACTCGCAGGTG"
            },
        "BbsI" : {
             "prefix" : "GAAGACTA",
             "suffix" : "ACGTCTTC"
             },
        "BfuAI" : {
            "prefix" : "ACCTGCCCTA",
            "suffix" : "ACTCGCAGGT"
            },
        "BsaI" : {
            "prefix" : "GGTCTCA",
            "suffix" : "CGAGACC"
            },
        "BsmBI" : {
            "prefix" : "CGTCTCA",
            "suffix" : "CGAGACG"
            },
        "BtgZI" : {
            "prefix" : "GCGATGCATGCTTGCA",
            "suffix" : "CCTCTGAATACATCGC"
            },
        "SapI" : {
            "prefix" : "GCTCTTCA",
            "suffix" : "GCTCTTCA"
            },
    }
    # Establish part types
    part_type = dict(
        cds = {
            # FreeGenes MoClo CDS definition
            "cloning enzyme": "BbsI",
            "small_cloning_enzyme" : "BtgZI",
            "retrieval_enzyme": "BsaI",
            "prefix": "A",
            "suffix": "AGAGCTT"
        },
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
    ## ==========================================================
    ## Add Prefix and Suffix
    ## ==========================================================
    if entry_type not in part_type.keys():
        print("not a valid type")
        return
    part = part_type[entry_type]
    full_sequence = standard_flanks["prefix"] + enzyme_configuration["BsaI"]["prefix"] + part["prefix"] + sequence + part["suffix"] + enzyme_configuration["BsaI"]["suffix"] + standard_flanks["suffix"] # Have a programmatic way to condense prefix / suffix rather than just defining them
    print(full_sequence)
    print("length",len(full_sequence))


## Test the function
import pandas as pd

data = pd.read_csv("./gene.csv")
for index,row in data.iterrows():
    fragment_gene(row["Sequence"],row["Part_Type"])












#
