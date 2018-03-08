def fragment_gene(sequence,entry_type):
    ## ==========================================================
    ## Configurations
    ## ==========================================================
    synthesis_configuration = {
        "max_length" : 1500,
        "min_length" : 300
    }
    pcr_configuration = {
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
        "BsaI-methylated-prefix" : {
            "prefix" : "CCGGTCTCA",
            "suffix" : "CGAGACC"
            },
        "BsaI-methylated-suffix" : {
            "prefix" : "GGTCTCA",
            "suffix" : "CGAGACCGG"
            },

    }
    # Establish part types
    part_type = dict(
        cds = {
            # FreeGenes MoClo CDS definition
            "cloning_enzyme": "BbsI",
            "small_cloning_enzyme" : "BtgZI",
            "retrieval_enzyme": "BsaI-methylated-suffix",
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
    ## Add Retrieval prefix and suffix
    ## ==========================================================
    if entry_type not in part_type.keys():
        print("not a valid type")
        return
    part = part_type[entry_type]
    retrieval_enzyme = part["retrieval_enzyme"]
    full_sequence = enzyme_configuration[retrieval_enzyme]["prefix"] + part["prefix"] + sequence + part["suffix"] + enzyme_configuration[retrieval_enzyme]["suffix"] # Have a programmatic way to condense prefix / suffix rather than just defining them

    ## ==========================================================
    ## Add Cloning prefix and suffix
    ## ==========================================================
    seq = standard_flanks["prefix"] + full_sequence + standard_flanks["suffix"]
    num_frags = len(seq) // synthesis_configuration["max_length"] + 1
    frag_len = len(seq) // num_frags

    frags = []
    if len(seq) < 300:
        cloning_enzyme = part["small_cloning_enzyme"]
        frag = enzyme_configuration[cloning_enzyme]["prefix"] + seq  + enzyme_configuration[cloning_enzyme]["suffix"]
        frags.append(frag)
    else:
        for i in range(num_frags):
            cloning_enzyme = part["cloning_enzyme"]
            frag = seq[max(0, i * frag_len - 2):min((i+1) * frag_len + 2,len(seq))]
            frag = enzyme_configuration[cloning_enzyme]["prefix"] + frag + enzyme_configuration[cloning_enzyme]["suffix"]
            frags.append(frag)
    return frags,cloning_enzyme




# ## Test the function
# import pandas as pd
#
# data = pd.read_csv("./gene.csv")
# for index,row in data.iterrows():
#     print (fragment_gene(row["Sequence"],row["Part_Type"]))












#
