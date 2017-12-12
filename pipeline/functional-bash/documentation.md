# Bionet synthesis KG bash scripts

find_recursive : Recursively gives another script each file in the directory. Sets global-sequence to 0.
- $1 = data directory
- $2 = type of file to find
- $3 = script to run
- $4 = second parameter for above script
- $5 = third parameter for above script
- $6 = fourth parameter for above script
- $7 = fifth parameter for above script

json_read.sh : Reads a section of a json file
- $1 file name
- $2 section to read
- $3 If = "write", write number to a file called json_read.txt

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

parameter_number.sh : Prints the total number of occurances of a pattern 
- $1 data directory
- $2 section to read
- $3 pattern to read

sequence_match.sh : Prints number of base pairs synthesized if a field is matched
- $1 data directory 
- $2 section to read
- $3 what to match

global.sh : Generalized commands to run on command line
- global.sh edit findif : search for a string in a json section, if it exists, write something to a different section 

