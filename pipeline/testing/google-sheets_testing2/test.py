import codon

table = codon.load_codon_table(taxonomy_id=4932)
print(codon.recode_sequence(table, "ATGGGCTAA", "GGC"))


