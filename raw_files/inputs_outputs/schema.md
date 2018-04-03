# FreeGenes process

## FreeGene retrieve.py (manages documentation and parts)
- Get orders and documentation
- Generate json documentation for Bionet 
- Generate list of Parts

## FreeGene fragmentation.py (manages fragment generation / linking)
- From list of Parts, generate Fragments compatible with Twist input

## FreeGene OpenFoundry.py (manages OpenFoundry upload)
- Upload list of Parts and Fragments -> OpenFoundry

## FreeGene Twist.py (manages Twist orders)
- Upload list of Fragments to Twist
- Track json/txt database of Twist orders
- When Twist plates arrive, upload Fragment-plate -> OpenFoundry



# OpenFoundry 
- Asks to begin builds whenever there are enough Fragments from Fragment-plates to make x number of Parts
- Data passed -> OpenFoundry (before sequencing files and the such, those are from suppliers)

## Part
- Part ID - primary key
- Part name
- Desired sequence 
- Status
- Mastermix information 

## Many - to - Many part/frag

## Fragment
- Fragment ID - primary key
- Fragment name
- Sequence?

## Fragment-plate 
- Fragment ID - primary key
- Fragment name

## Mastermix (limited quantity, manually written)
- Mastermix ID - primary key
- Vector information 
- Enzyme concentrations
- Fragment / Part concentrations










