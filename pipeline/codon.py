import pandas as pd
import numpy as np
import os
import logging

CODON_USAGE_DB = os.path.dirname(__file__) + "/codon_table_data/codon_usage.spsum"
COMMON_SPECIES = {
    'E.coli': "83333",
    'S.cerevisiae':  "4932",
    'H.sapiens': "9606",
    'B.subtilis': "1432"
}

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

codons = ['CGA', 'CGC', 'CGG', 'CGT', 'AGA', 'AGG', 'CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG', 'TCA', 'TCC', 'TCG', 'TCT', 'AGC', 'AGT', 'ACA', 'ACC', 'ACG', 'ACT', 'CCA', 'CCC', 'CCG', 'CCT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'AAA', 'AAG', 'AAC', 'AAT', 'CAA', 'CAG', 'CAC', 'CAT', 'GAA', 'GAG', 'GAC', 'GAT', 'TAC', 'TAT', 'TGC', 'TGT', 'TTC', 'TTT', 'ATA', 'ATC', 'ATT', 'ATG', 'TGG', 'TAA', 'TAG', 'TGA']
standard_genetic_code = ['R', 'R', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L', 'L', 'L', 'S', 'S', 'S', 'S', 'S', 'S', 'T', 'T', 'T', 'T', 'P', 'P', 'P', 'P', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', 'K', 'K', 'N', 'N', 'Q', 'Q', 'H', 'H', 'E', 'E', 'D', 'D', 'Y', 'Y', 'C', 'C', 'F', 'F', 'I', 'I', 'I', 'M', 'W', '*', '*', '*']
tables = dict()

def load_codon_table(species=None, taxonomy_id=None):
    """Load a codon table based on the organism's species ID."""
    print("in codon_table")
    if taxonomy_id in COMMON_SPECIES:
        print("found in COMMON_SPECIES")
        taxonomy_id = COMMON_SPECIES[taxonomy_id]
        print("New id: ",taxonomy_id)
    print("in codon_table")

    # if species in COMMON_SPECIES:
    #     taxonomy_id = COMMON_SPECIES[species]

    taxonomy_id = str(taxonomy_id)

    if taxonomy_id in tables:
        logger.debug("Returning codon table from cache.")
        return tables[taxonomy_id]

    logger.debug("Loading codon table from {}.".format(CODON_USAGE_DB))
    with open(CODON_USAGE_DB) as f:
        for header in f:
            codon_counts = f.readline()

            taxid, species, cds_number = header.strip().split(":")[:3]

            if taxonomy_id and taxonomy_id != taxid:
                continue

            logger.debug("Loaded {} {}".format(taxid, species))

            print("before first table")
            table = list(zip(codons, standard_genetic_code, [int(x) for x in codon_counts.split()]))
            print("after first table")
            table = pd.DataFrame(table, columns=['Triplet', 'AA', 'Number'])
            table.set_index(['AA', 'Triplet'], inplace=True)
            table.sort_index(inplace=True)

            table['Fraction'] = table.groupby('AA').transform(lambda x: x / x.sum())

            tables[taxid] = table
            tables[species] = table
            break

    return table

def codon_table_10plus(table):
    """Return a codon table only representing codons with > 10% occurrence frequency."""
    print("in 10plus")
    table = table.ix[table.Fraction >= 0.1]
    table = table.groupby(level=0).transform(lambda x: x / x.sum())

    return table

def reverse_complement(table, seq):
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(str.maketrans("ATGC","TACG"))[::-1]

def recode_sequence(table, seq, rep):
    """Recode a sequence to replace certain sequences, using a given codon table."""
    pos = seq.find(rep)

    if pos < 0:
        return seq

    pos -= pos % 3

    for i in range(pos, pos + (len(rep) // 3 + 1) * 3, 3):
        codon = seq[i:i+3]
        choices = codon_table_10plus(table).ix[table.ix[codon]]
        choices = choices[choices.index != codon]

        if choices.shape[0] > 0:
            newcodon = choices.iloc[(choices.Fraction.cumsum() / choices.Fraction.cumsum().max() < np.random.rand()).sum()].name # Stochastically allocate codon
#             newcodon = choices[choices.index != codon].Fraction.idxmax() # Deterministically allocate codons by using the most frequence one
            break

    logger.warn("{} -> {}".format(codon, newcodon))

    return seq[:i] + newcodon + seq[i+3:]

def remove_cutsites(cut_sites, seq):
    """Remove cutsites from a sequence."""
    changes = 0

    for enzyme, cut in cut_sites + [(e, reverse_complement(c)) for e, c in cut_sites]:
        while cut in seq:
            logger.warn("{} cuts ({})".format(enzyme, cut))
            changes += 1
            seq = recode_sequence(seq, cut)

    return seq, changes

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    import IPython; IPython.embed()
