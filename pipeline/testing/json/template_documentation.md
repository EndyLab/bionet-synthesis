# Template JSON Documentation

## Version
* Describes the current verison of the template. Utilizes Semver 2.0 standard for versioning

## Gene ID
* Unique identifier for each gene that is submitted
* Begins with "BBF10K\_" and then increments in number with padding to 6 significant figures (i.e. BBF10K_000004)

## Author
* Takes in strings for:
  * **the author name** (first and last)
  * **email**
  * **affiliation**
  * **ORCID** - Optional

## Info
### Documentation
* **Gene name** - the user specified name of the part
* **Who** - Answer to "Who would be interested or who would they be useful to?"
* **What** - Answer to "What exactly is this gene? What does it do?"
* **Where** - Answer to "Where does this gene come from? How was it discovered?"
* **Why** - Answer to "Why are you synthesizing this gene?"
* **Future Applications** - Answer to "What are some future applications of this gene?"
* **Database Links** - Contains a list with available links as strings

### Gene Metadata
#### Cloning
* **Part Type** - Strings listed below are the exact accepted entry types
  * 'cds'
  * 'eukaryotic_promoter'
  * 'prokaryotic_promoter'
  * 'rbs'
  * 'terminator'
  * 'vector'
* **Cloning Enzyme** - Enzyme used to cut the synthesized fragments and clone into the standard pOpen vector
* **Retrieval Enzyme** - Enzyme used to cut the part out of the standard vector and yield MoClo compatible sticky ends
* **Optimize** - Boolean value stating whether or not to optimize. True means to optimize.
* **Target Organism** - Specifies which organism's codon table to look up with codon optimizing
  * **Organism name** - Human readable name specifying the organism the part was optimized for
  * **Taxid** - String containing the number corresponding to the desired organism

#### Collection ID
> ID number used to group batch submissions together

* Every submission is given a unique collection ID such that we can pull out a complete submission
* *ADD IN EXAMPLE OF A COLECTION ID*

#### IP
> Specifically refers to the results found in the Sagacious IP search

* Takes entries for the submission number and the outcome of that search
* *WILL NEED TO UPDATE WHEN WE GET RESULTS BACK FROM SAGACIOUS*

### Safety
* We assume that Twist will catch any harmful sequences and so this is just a spot for a user to voice any concerns about the sequence or its use

### Order Number
* Specifies which order this part was sent out in.
* Takes the form "submission000" and increments sequentially with three significant figures (i.e. submission_008)

## Sequence
> Stores all relevant sequences associated with the part

* **Original sequence** - The sequence that was originally submitted
* **Optimized sequence** - The sequence after the optimization
* **Fragment sequences** - A dictionary with the fragment names as keys with their sequences as their corresponding keys

## Status
### **Current status**
> Stores a string that corresponds to a unique state depicting the status of the part at a give time

1. **'submitted'**- Starting state after being entered into the database
2. **'optimized'** - The user has cleared the changes made during the optimization
3. **'ordered'** - The order has been sent to Twist
4. **'synthesis_abandoned'** - *ENDPOINT* - Twist wasn't able to synthesize one of the fragments
5. **'received'** - All of the fragments need to assemble the part have been received
6. **'building'** - The part is currently in production
7. **'sequencing'** - The cloning run has finished and the part is out for sequencing
8. **'sequence_confirmed'** - *ENDPOINT* - The part has been cloned and successfully sequence verified
9. **'cloning_mutation'** - The sequence revealed a small flaw in the sequence. Another colony will be sequenced
10. **'bad_reads'** - The read quality was too low to deduce a sequence
11. **'cloning_failed'** - The sequence yielded the original vector sequence, a partial sequence, or some other unintended sequence
12. **'cloning_abandoned'** - *ENDPOINT* - The part has failed to clone after two rounds of cloning and so it has been abandoned

### Build Attempts
> A list of dictionaries, each entry representing a different build attempt

* **Building** - Boolean indicating if it is currently being built
* **Build well** - The well containing the part in the 96-well plate
* **Build date** - Specifies the day that the build was run
* **Build number** - Takes the form "build000" and increments sequentially with three significant figures (i.e. build010)
* **Build outcome** - The results from the sequence verification, will be used to update current_status. This result will provide more detail to failures than current_status. Outcomes include:
  * 'squence_confirmed' - Same as current status
  * 'cloning_mutation' - Same as current status
  * 'bad_reads' - Same as current status
  * 'original_vector' - The sequence matches the destination vector used in the cloning reaction
* **Forward read** - The name of the forward sequencing file used in the alignment
* **Reverse read** - The name of the reverse sequencing file used in the alignment

## Location
* **Fragments** - A nested dictionary with the fragment names as the keys and their location as the corresponding values. The locations are stored as strings containing both the plate and well that the fragment is located in (i.e. '"plate"\_"well"')
* **Cloned** - Contains a list of dictionaries that correspond to the instances of verified successful clones. Has not actually been used in any of the code thus far
  * **Vector** - The vector that the part was cloned into
  * **Organism** - The organism that is propagating this plasmid
  * **Location** - Specifies a list of plates and wells that this strain is in

## Dates
> All dates are in the following format: "Year.Month.Day"

* **Submitted** - Date when the entry was submitted through the google form
* **Ordered** - Date when the part was ordered
* **Received** - Date when all of the fragments have been received
* **Completed** - Date when the part was sequence verified to be correct
