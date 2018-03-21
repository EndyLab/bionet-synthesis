Template JSON Documentation.

See `template.json` for example.

Leaf node data types are all `string` unless otherwise specified.

# version

Describes the current verison of the template. Utilizes Semver 2.0 standard for versioning

# gene_id

Unique identifier for each gene that is submitted.

Begins with "BBF10K\_" and then increments in number with padding to 6 significant figures (i.e. BBF10K_000004).

# author
* `author_name` (first and last)
* `email`
* `affiliation`
* `orcid` - Optional

# info

## info.documentation
* `gene_name` - the user specified name of the part
* `who` - Answer to "Who would be interested or who would they be useful to?"
* `what` - Answer to "What exactly is this gene? What does it do?"
* `where` - Answer to "Where does this gene come from? How was it discovered?"
* `why` - Answer to "Why are you synthesizing this gene?"
* `future_app` - Answer to "What are some future applications of this gene?"
* `database_links` - Contains a list with available links as strings

## info.gene_metadata

### info.gene_metadata.cloning
* `part_type` - Strings listed below are the exact accepted entry types
* `cloning_enzyme` - Enzyme used to cut the synthesized fragments and clone into the standard pOpen vector
* `retrieval_enzyme` - Enzyme used to cut the part out of the standard vector and yield MoClo compatible sticky ends
* `optimize` - Boolean value stating whether or not to optimize. True means to optimize.

The following are valid values for `part_type`:

* `'cds'`
* `'eukaryotic_promoter'`
* `'prokaryotic_promoter'`
* `'rbs'`
* `'terminator'`
* `'vector'`

#### info.gene_metadata.cloning.target_organism 

Specifies which organism's codon table to look up with codon optimizing

* `organism_name` - Human readable name specifying the organism the part was optimized for
* `taxid` - String containing the number corresponding to the desired organism

### info.gene_metadata.collection_id

ID number used to group batch submissions together

Every submission is given a unique collection ID such that we can pull out a complete submission

TODO: Add example of collection_id

### info.gene_metadata.IP

Specifically refers to the results found in the Sagacious IP search.

Takes entries for the submission number and the outcome of that search.

TODO: Update this when we get results back from Sagacious.

### info.gene_metadata.safety

We assume that Twist will catch any harmful sequences and so this is just a spot for a user to voice any concerns about the sequence or its use.

## info.order_number

Specifies which order this part was sent out in.

Takes the form "submission000" and increments sequentially with three significant figures (i.e. submission_008).

# sequence

Stores all relevant sequences associated with the part.

* `original_sequence` - The sequence that was originally submitted
* `optimized_sequence` - The sequence after the optimization
* `fragment_sequences` - An array with the fragment names as keys with their sequences as their corresponding keys

# status
## current_status

Stores a string that corresponds to a unique state depicting the status of the part at a give time.

1. `'submitted'`- Starting state after being entered into the database
2. `'optimized'` - The user has cleared the changes made during the optimization
3. `'ordered'` - The order has been sent to Twist
4. `'synthesis_abandoned'` - *ENDPOINT* - Twist wasn't able to synthesize one of the fragments
5. `'received'` - All of the fragments need to assemble the part have been received
6. `'building'` - The part is currently in production
7. `'sequencing'` - The cloning run has finished and the part is out for sequencing
8. `'sequence_confirmed'` - *ENDPOINT* - The part has been cloned and successfully sequence verified
9. `'cloning_mutation'` - The sequence revealed a small flaw in the sequence. Another colony will be sequenced
10. `'bad_reads'` - The read quality was too low to deduce a sequence
11. `'cloning_failed'` - The sequence yielded the original vector sequence, a partial sequence, or some other unintended sequence
12. `'cloning_abandoned'` - *ENDPOINT* - The part has failed to clone after two rounds of cloning and so it has been abandoned

## build_attempts

An array, each entry representing a different build attempt.

Each entry in the array has the following entries:

* `building` - Boolean indicating if it is currently being built
* `build_well` - The well containing the part in the 96-well plate
* `build_date` - Specifies the day that the build was run
* `build_number` - Takes the form "build000" and increments sequentially with three significant figures (i.e. build010)
* `build_outcome` - The results from the sequence verification, will be used to update current_status. This result will provide more detail to failures than current_status. Valid values specified below.
* `Forward read` - The name of the forward sequencing file used in the alignment
* `Reverse read` - The name of the reverse sequencing file used in the alignment

Valid values for `build_outcome` are:

* `'sequence_confirmed'`` - Same as `current_status`
* `'cloning_mutation'` - Same as `current_status`
* `'bad_reads'` - Same as `current_status`
* `'original_vector'` - The sequence matches the destination vector used in the cloning reaction

# location

* `fragments` - An object with the fragment names as the keys and their location as the corresponding values. The locations are stored as strings containing both the plate and well that the fragment is located in (i.e. '"plate"\_"well"').

## location.cloned

Contains an array where each entry corresponds to the instances of verified successful clones. Has not actually been used in any of the code thus far.

Each array entry has the following properties:

* `vector` - The vector that the part was cloned into
* `organism` - The organism that is propagating this plasmid
* `location` - An array of plates and wells that this strain is in

TODO specify properties of the entries in `location`.

# dates

All dates are in the format `YYYY.MM.DD`

* `submitted` - Date when the entry was submitted through the google form
* `ordered` - Date when the part was ordered
* `received` - Date when all of the fragments have been received
* `completed` - Date when the part was sequence verified to be correct
