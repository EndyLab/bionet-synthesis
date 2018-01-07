# Pipeline tools for gene tracking and cloning using the OT-one

## Quickstart
1. Populate: `cat 10K_CDS.csv & template.json | python populate.py > dir BBF10K_`
  Creates new directories and json files for each new gene to be synthesized and updates the template.json file with all of the available data present in the .csv file provided
2. 




## Populate
Populate does the following:
* Reads in a .csv file detailing all of the currently known information about each gene. Information should include:
  * gene_name
  * sequence
  * author
  * author_email
  * author_affiliation
  * author_project
  * cloning_method
  * part_type
  * build_type
  * safety
  * collection
  * other_tags
  * ordered
  * will_build
  * date_ordered
  * order_number
  * idnum
    * New submissions should enter "None" in this space
* Reads in the template.json files to initialize file with the correct starting conditions
* Assess if the gene has already been entered by checking the "idnum" column and will only create new directories and files for rows that have "None" in the idnum column
* Assigns a the next unused id number to each new entry
* Generates a new directory titled with the id number
* Creates a new FASTA file using the sequence data provided
* Also creates a json file based on the template provided and appeding all of the new Information
* Updates the .csv file with all of the new id numbers












* something
