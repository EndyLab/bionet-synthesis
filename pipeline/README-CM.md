# Pipeline tools for gene manipulation, tracking and cloning using the OT-one

## Quickstart
0. Compile all of the submission information in the 10K_CDS.csv file
1. Populate: `cat 10K_CDS.csv & template.json | python populate.py > dir BBF10K_`
  Creates new directories and json files for each new gene to be synthesized and updates the template.json file with all of the available data present in the .csv file provided.
2. Optimize Genes: `python optimize-genes.py > Updated json files`
   This removes invalid restriction sites, converts stops to TGA for MoClo, and sets the penultimate codon to our fixed list
3. Fragment Genes: `python fragment_genes.py > Updated json files`
   This breaks DNA into <1.8kb fragments for synthesis, adds assembly cut sequences, and pairs small genes together. Output is in Twist's order format.
4. Add Fragments: `cat 10K_CDS.csv & fragments1-5.csv | python add_frag.py > Updated json files`
  ONLY FOR FIRST ORDER. Takes in the .csv file detailing all of the fragments for each gene and then creates keys for each in their corresponding .json file.
5. Fragment Location Assignment: `cat 10K_CDS.csv | python frag_loc_assign.py > Updated json files`
  Takes in a plate map from the plate_maps directory and assigns the locations of each fragment to their corresponding json file.
6. Validate Fragments: `python validate_fragments.py`
  Preforms an in silico cloning reaction to determine if the fragments will form the desired sequence.
7. Abandoned: `python abandoned.py`
  Determines which genes are now abandoned based on the last completed order.
8. Database Status: `python db_status.py > sankey.txt`
  Queries the database and counts occurrences of specific useful criteria and outputs a text file to generate a Sankey diagram.
9. Submission Status: `python sub_status.py`
  Same as the db_status script but only searches for genes from a specific order.
10. Plate Resuspension: `python plate_resuspension.py`
  Calculates the amount of water to required to resuspend all of the wells of a specific plate to the desired concentration.
11. Build Golden Gate Reactions: `python build.py > build_.csv `
  Queries the database for buildable genes, specifies the required reagents and setup, executes the OT-One protocol, and generates a csv detailing the new plate map.
12. Transformation Plating: `cat build_.csv | python plating.py > agar_maps.csv`
  Takes in the build_map and executes a defined dilution series with the OT-One to plate all of the transformations
13. Sequence Data Sorter: `cat 10K_CDS.csv & Sequencing files | python seq_sort.py`
  Transfers all of the sequencing files to their corresponding directories.
14. Database Sequence Alignment: `python seq_align_db.py`
  Generates a sequence alignment of all of the sequence file pairs and the corresponding FASTA file in each directory and scores the alignment to determine if it is a successful clone or not.
15. Sagacious Submission: `python sagacious_input.py`
  Pulls genes from the database and generates a new csv file to be sent off to Sagacious for review.
16. Modify JSON: `python modify_json.py`
  Short script for modification of the json file for each gene.


## Populate
Populate does the following:
* Reads in a .csv file detailing all of the currently known information about each gene.
* Reads in the template.json files to initialize file with the correct starting conditions
* Assess if the gene has already been entered by checking the "idnum" column and will only create new directories and files for rows that have "None" in the idnum column
* Assigns a the next unused id number to each new entry
* Generates a new directory titled with the id number
* Creates a new FASTA file using the sequence data provided
* Also creates a json file based on the template provided and appeding all of the new Information
* Updates the .csv file with all of the new id numbers

## Optimize Genes
NEED TO MAKE IT INTERFACE WITH THE DATABASE
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

## Fragment Genes
NEED TO MAKE IT INTERFACE WITH THE DATABASE
* Add an 'A' nucleotide before and after each sequence, to match CDSes for MoClo (aATG / TGAa)
* Break big genes into fragments
  * We break genes into approximately equal pieces based on their size
  * We add BbsI recognition sequences to the beginning and end of each sequence
    Note that this assumes that a gene broken into three fragments won't have identical overhangs on, for example, the end of fragments 1 and 2 (so they could be swapped on assembly)
* Pair small (<300 bp) genes so they go above Twist's minimum synthesis limit
  * We pair the biggest small gene with the smallest (eg, 299 with 50, then 250 with 70, if fragments of those sizes existed).
  * The genes are connected with a linker for BbsI / BtgZI gene selection on assembly
* We output the fragments in the format Twist wants for their order spreadsheet
  * Genes are named `[Gene Name]_[Fragment Number]`
  * Small genes are named `[Gene 1 Name]_link_[Gene 2 Name]_1`

## Add Fragments
add_frag does the following:
* ONLY FOR THE FIRST 5 SUBMISSIONS
* Reads in the .csv file containing all of the current data used to populate the database and builds a dictionary with the gene name and ID#
* Pulls in the .csv file detailing all of the fragments for each gene
* Splits up the linked genes
* Iterates over the fragments, searches for the corresponding directory and json file, and appends the fragment name and sequence

## Fragment Location Assignment
frag_loc_assign does the following:
* Reads in the .csv file containing all of the current data used to populate the database and builds a dictionary with the gene name and ID#
* Pulls in the plate maps that are available and requests which one to use
* Iterates through the rows of the dataframe and searches for multi-fragment parts
* Rebuilds the dataframe with the multiple fragments entered as separate items
* Uses the dictionary to convert fragments titled with its gene name to ID# instead
* Searches for the corresponding directory and .json file and then appends the fragment location and its sequence to the file
* Determines if the gene is buildable by counting the number of fragments with locations and comparing it to the number of fragments that there should be
* Reports all of the genes that were not found in the dictionary or the directories

## Validate Fragments
NEED TO MAKE IT INTERFACE WITH THE DATABASE
The validation script, `validate-fragments.py`, takes Twist-format input from
`fragment-genes.py` and does a virtual digestion and assembly. It outputs
expected gene sequences for cross-checking with the input. It's pretty rough and
ready, but should give a good indication of whether the fragmentation worked
properly.

## Abandoned
abandoned does the following:
* Opens up every file in the database and checks if the gene was ordered prior to the last completed order and if it is marked buildable.
* If the it is not buildable and was ordered earlier than the last completed order it is marked as true under the key abandoned.

## Database Status
db_status does the following:
* Queries every gene's json file in the database
* Counts the number of:
  * Genes in each submission
  * Ordered genes
  * Buildable genes
  * Build attempts
  * Sequence verified clones
  * Original vector clones
  * Unknown sequence clones
  * Contributors
* Generates a txt file that incorporates the data in a form that can be copied and pasted in the http://sankeymatic.com/build/ page to generate a Sankey diagram.
  * Change the settings to Width:1000px, Vertical_space:30px, and Node_width:0px

## Submission Status
sub_status does the following:
* Queries every gene's json file in the database, but only counts genes from the specified submission.
* Counts the number of:
  * Genes in each submission
  * Ordered genes
  * Buildable genes
  * Build attempts
  * Sequence verified clones
  * Original vector clones
  * Unknown sequence clones
  * Contributors

## Plate Resuspension
CHANGE HOW IT IS PRINTING THE DECK MAP
plate_resuspension does the following:
* Takes in the argument `-r` to tell it to send the protocol to the OT-one
* Queries the plate_maps directory and asks which submission to look into and then displays all of the plates present in that csv file
* Calculates the amount of water to resuspend the DNA in each well to the desired concentration
* Prints out the required deck map
* Tells the OT-One what volume of water to add into each well
  * Currently it is not set to mix because that dramatically increases the run time and we can just seal it and then vortex to resuspend
* Reports the runtime at the end of the protocol

## Build Golden Gate Reactions
Build will do the following:
* Takes in the argument `-r` to tell it to run the protocol on the OT-one
* Takes in the argument `-m` to tell it to manual adjust the parameters for the selection of genes
* Queries the database for genes marked as buildable
  * Will either pull from user specified plates or it will go through and pull from the first set of plates that it come across
  * Search criteria is bounded by the number of plates that you want to pull from (max = 3), number of reactions (best to stick to 96), and the max number of fragments you want to assemble in a single reaction
* Generates two dataframes, one specifying the amount of master mix to add to each well and the other to specify where each fragment needs to be transferred to
* Displays a layout table detailing where everything needs to be placed on the deck
* Displays the required ingredients to make the required amount of master mix
* Commands the OT-One to transfer the master mix to a PCR tube strip using the p200 and then uses the p10 multichannel to aliquot the liquid into the wells
  * An extra reaction is added to each well for each extra fragment added. This allows the a single master mix to be used for all reactions
* Transfers each of the fragments to their specified well
  * PULL OFF THE ALUMINUM FOIL FROM THE SOURCE PLATES AND REPLACE WITH THE AEROSEAL
* Looks for previous build maps to determine which build number this run is
* Asks the user if it was successful and should be recorded. If successful:
  * Generates a build map linking the wells with the constructs
  * Updates the database with the information from the run  
* Reports the runtime at the end of the protocol


## Transformation Plating
CHANGE THE DILUTION CALCULATION TO ALLOW FOR LARGER DILUTIONS AND HAVE IT GENERATE A PLATE MAP FOR THE AGAR PLATES
Plating will do the following:
* Takes in the argument `-r` to tell it to run the protocol on the OT-one
* Pulls the build maps available in the build directory and asks the user which build they would like to plate
* Displays the required layout for the OT-One and tells the user how many agar plates will be required and what to name them
* Repeatedly plates and dilutes the cells until it has completed the serial dilution
* Reports the runtime at the end of the protocol


## Sequence Data Sorter
seq_sort does the following:
* Reads in the .csv file containing all of the current data used to populate the database and builds a dictionary with the gene name and ID#
* Iterates through all of the sequencing files that are present in the sequencing_files directory and copies all of the absorbance files to the corresponding data directory

## Database Sequence Alignment
db_seq_align does the following:
* This program assumes that the sequencing files have already been sorted into their respective directories and that
* Takes in numerous initial strings:
  * Assumes that the names of the primers are, M13-Forward---20- and M13-Reverse based on the names given by Elim
  * Takes in paths for the the locations of the sequencing files, their corresponding json files and the path to the sequence of the original vector
* Generates a dictionary relating the gene names to their corresponding gene ID#s
* Searches for all sequencing files within the data directories and parses the name to link it to a specific gene and determine which primer was used to generate it

## Sagacious Submission
sagacious_input does the following:
* Goes through the database and selects genes that have been successfully sequence verified and have not been previously submitted
* Adds all of those genes to a dataframe and then returns a new csv file
  * Files are titled: "sagacious_order_order_number_date"
* Updates the json files of the submitted genes with the order number that it was submitted in

## Modify JSON
modify_json does the following:
* Opens every gene directory and JSON file currently in the database
* Modifies the JSON file however the user specifies and rewrites the file with the edits

* space
