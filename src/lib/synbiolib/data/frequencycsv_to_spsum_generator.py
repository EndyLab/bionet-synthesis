import pandas as pd

# Grab table from user
table = input("Desired table : ")

# Sets up spsum format from SPSUM_LABEL ftp://ftp.kazusa.or.jp/pub/codon/current/SPSUM_LABEL
# Stop codons are removed, added at end
spsum_format = "CGA CGC CGG CGT AGA AGG CTA CTC CTG CTT TTA TTG TCA TCC TCG TCT AGC AGT ACA ACC ACG ACT CCA CCC CCG CCT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT AAA AAG AAC AAT CAA CAG CAC CAT GAA GAG GAC GAT TAC TAT TGC TGT TTC TTT ATA ATC ATT ATG TGG".split()

# Read in table as pandas dataframe
csv_data = pd.read_csv("frequency_example.csv")
for index, row in csv_data.iterrows():
    spsum_format = [ x if not x == row["codon"] else str(int(float(row[table])*100)) for x in spsum_format ]

# Add in stop codons
spsum_format.extend(['100', '0', '0'])
properspsum = ' '.join(spsum_format)

# Add taxid (custom_x)
taxid = input("Taxid of table? : ")


# Write table
with open("custom_table.spsum", "a") as custom_table:
    custom_table.write(taxid + ":" + table + ":" + str(sum(list(map(int, spsum_format)))) + "\n" + properspsum)

