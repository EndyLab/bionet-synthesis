# Bionet synthesis KG bash scripts

find_recursive : Recursively gives another script each file in the directory. Sets global-sequence to 0.
- $1 = data directory
- $2 = type of file to find
- $3 = script to run
- $4 = second parameter for above script
- $5 = third parameter for above script

json_read.sh : Reads a section of a json file
- $1 file name
- $2 section to read

json_write.sh : Writes a section of a json file
- $1 file name
- $2 section to edit
- $3 edits to section

sequence_retriever.sh : Retrieves sequence from FASTA file
- $1 file name

sequence_length.sh : Prints length of FASTA file
- $1 file name
- $2 If = "write", write number to a file called sequence-total.txt

sequence_length-global.sh : Prints final number of base pairs from directory
- $1 data directory


