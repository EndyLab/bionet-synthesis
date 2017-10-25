# Pipeline tools for Twist synthesis

## Quickstart
1. Get yo' gene coding sequences into a CSV with the format Gene,Sequence
2. Optimize: `cat genes.csv | python optimize-genes.py > optimized.csv`
   This removes invalid restriction sites, converts stops to TGA for MoClo, and sets the penultimate codon to our fixed list
3. Fragment: `cat optimized.csv | python fragment-genes.py > fragments.csv`
   This breaks DNA into <1.8kb fragments for synthesis, adds assembly cut sequences, and pairs small genes together. Output is in Twist's order format.


## Optimize
Optimize does the following:
* Runs basic tests for proper coding sequences:
  * Sequence isn't empty
  * Sequence is a triplet
  * Sequence begins with start codon
  * Sequence ends with stop codon
  * Sequence doesn't contain internal stops
* Forces stop codons to TGA
  * We need this for the TGAA MoClo site
* Optimizes the sequence to meet several requirements:
  * Remove RE sites for BfuAI, AarI, BtgZI, BbsI, BsmBI, SapI, BsaI
  * Remove homopolymer runs of >= 6 base pairs
  * Checks (but does not attempt to fix) out of spec GC
* Forces the final (non-stop) codon to our fixed list for SapI cloning compatibility
  * We then check if this change introduced a new RE site: we don't attempt to fix if this is so, but throw an error

## Fragment
* Add an 'A' nucleotide before and after each sequence, to match CDSes for MoClo (aATG / TGAa)
* Break big genes into fragments
  * We break genes into approximately equal pieces based on their size
  * We add BbsI recognition sequences to the beginning and end of each sequence
    Note that this assumes that a gene broken into three fragments won't have identical overhangs on, for example, the end of fragments 1 and 2 (so they could be swapped on assembly)
* Pair small (<300 bp) genes so they go above Twist's minimum synthesis limit
  * We pair the biggest small gene with the smallest (eg, 299 with 50, then 250 with 70, if fragments of those sizes existed).
  * The genes are connected with a linker for BbsI / BtgZI gene selection on assembly
* We output the fragments in the format Twist wants for their order spreadsheet
  * Genes are named [Gene Name]_[Fragment Number]
  * Small genes are named [Gene 1 Name]_link_[Gene 2 Name]_1

## Validation
The validation script, `validate-fragments.py`, takes Twist-format input from
`fragment-genes.py` and does a virtual digestion and assembly. It outputs
expected gene sequences for cross-checking with the input. It's pretty rough and
ready, but should give a good indication of whether the fragmentation worked
properly.
