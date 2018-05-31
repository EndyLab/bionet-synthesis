# Tutorial for using the OpenFoundry

## Overall Workflow
0. Generate a list of all of the sequences to be synthesized
- [ ] Provide a template csv for adding new sequences
1. Optimize all of the sequences present using the optimization script
2. When ready to order, fragment genes and send for sequencing
- [ ] Provide information about why we are fragmenting and what the options are
3. When the plates come in run the location assignment script to parse the newly added plate maps and add locations to all of the fragments. This will also determine which genes are now able to be built.
4. When the necessary fragments arrive, run the build script
5. Once the golden gate assemblies are made, place in the thermocycler and run the designated program
6. Run the transformation protocol
7. Run the plating protocol
8. Place plates in the 37C incubator overnight
9. Retrieve the plates and photograph them within the image box
10. Fill up a 96 deep well plate with desired media and run the colony picking script
11. Place the deep well block in a plate shaker at 37C overnight
12. Either run a 96 well miniprep or send the whole block out for sample preparation and sequencing
13. If using Sanger sequencing, use the reads to run sequence alignments against all of the corresponding target sequences and determine the outcomes of the cloning reactions

## General information for OpenTrons code
- All scripts interacting with the OpenTrons share a similar process and ideology
  - The intent is for the robot to keep track of where objects should be and have it tell the user what materials it needs and where to put them.
    - This liberates the user from having to track which plates and wells are important and eliminates the need to translate that information to the robot.
  - When any of the OpenTrons scripts are run alone they will simulate the run by generating a virtual robot to send commands too.
    - It is recommended to quickly simulate the run prior to full execution to make sure that the robot will to the indended task.
  - Once the simulated run is complete pass the `-r` argument after the name of the script to tell it to run the protocol.
  - All of the scripts will print out a deck layout that instructs the user where all of the needed materials need to be placed.


## Optimization

### Purpose

### Specifics

### Code Breakdown



## Fragmentation

### Purpose

### Specifics

### Code Breakdown



## Fragment Location Assignment

### Purpose
Assigns the well locations of all of the received DNA fragments to the corresponding fragments in the database such that when building golden gate reactions, the OpenTrons knows which wells to pull from to make a certain gene.

### Requirements
- All new plate maps must be added to the plate_maps directory

### Code Breakdown
1. Queries the database to find all of the existing plates that have been entered thus far
2. Iterates through all of the plate_maps in the plate_maps directory
3. Reads the csv file in as a pandas data frame
4. Determines which version of plate map it is
  - Currently checks for two different versions that Twist Biosciences uses.
  - If the version doesn't match the plate maps being used change the script such that the keys correspond to the column names in the plate map.
  - You will specifically need the plate ID number, well location, name of the fragments, the sequence of the fragment, and the synthesis yield
5. Iterates through each of the unique plates in the plate map
6. Generates a new Plate object in the database
7. Goes through all of the fragments in each plate and searches the database for the fragment object that matches the sequence of the fragment in the plate map
8. When found the fragment is added to the well listed in the plate map to the plate object
  - If not found, the fragment is added to a list that will be present at the end so that you can trouble shoot why the fragments are not in the database
9. Goes through all of the parts and checks all of their corresponding fragments to determine which parts have all of their fragments and can therefore be constructed


## Plate Resuspension

### Purpose
Resuspend all of the fragments in a synthesis plate such that all of the fragments are at the same concentration

### Requirements
- Plate to resuspend
- Tubes of water
- 1x aluminum seal

### User Intervention
1. Pull the seal off of the plate to be resuspended and place it on the deck where specified by the script
2. Fill up as many tubes of water as is specified
3. After resuspension, reseal with a new aluminum seal and vortex lightly to ensure the DNA is resuspended

### Code Breakdown
0. If run directly, it will return a list of all of the plates that are currently not suspended
  - Must enter `-r` after the script name for the script to connect to the robot, otherwise it simulates the protocol
1. Determines the deck layout and initializes the OpenTrons
2. Calculates the volume needed to resuspend each fragment
  - Currently calculates the volume based on the desired concentration of 40 femtomoles/µL
3. Sums the volume needed for each fragment to determine how much water is needed to carry out the protocol
4. Aliquots the calculated volume of water to each well
  - Tracks the volume in each tube of water so it knows when to switch tubes
  - Changes pipettes depending on the volume to be pipetted


## Golden Gate Assembly

### Purpose
Determines which genes to build, calculates the needed master mix and instructs the robot to add all of the fragments needed to assemble each gene

### Requirements
- The necessary enzymes and buffers
  - Run a simulation of the build to see the golden gate reaction mix needed
- Genes that have all of the necessary fragments
- 2x 48 well PCR blocks
- 12x 8-lid strips
- 1x 8-tube PCR strip
- 1-2x 10µL tip racks
- 1x 200µL tip rack
- Tubes of water

### User Intervention
1. Mark the two 48 well PCR blocks with the build number, note the orientation of the plates so you know where the 'A1' position is, and add them to the slot labelled `DEST_PLATE`.
2. Add the PCR tube strip to the first row in the rack labelled `PCR-strip-tall` in the deck layout.
3. Take the synthesis plates required for the build out of the -20C freezer to thaw.
4. Once thawed, vortex lightly to mix the resuspensions
5. Replace the aluminum seals with permeable seals
  - We cover the plates with permeable seals to dramatically decrease the evaporation rate by covering the wells but also allows the Opentrons to puncture through with ease
  - The OpenTrons can stab through the aluminum foil directly but the aluminum forms a tight seal with the tip as it penetrates the foil making it far more difficult for the pipette to aspirate the media, which is why we use the permeable seal
3. Assemble the master mix based on the volumes provided by the script, mix it well and add it to the 'A1' position (the bottom left spot) of the centrifuge tube rack
  - The volumes of each reagent for an 8µL reaction of master mix is detailed below:

| Reagent | Amount per reaction     |
| ------------- | :-------------: |
|   Cutsmart Buffer    |    1 µL    |
|   ATP    |    1 µL    |
|   Vector (100ng/µL)    |    0.25 µL    |
|   T4 DNA Ligase    |    0.5 µL    |
|   Restriction Enzyme    |    0.25 µL    |
|   Water      |   5 µL    |

4. Add a microfuge tube of water to the 'B1' position in case any of the fragments must be diluted.
5. When a synthesis plate is pipetted from for the first time, it will ask the user to check the calibration and allows the user to move the pipette up or down.
6. Once the run is completed, cover the PCR tube plates with the lids, flick mix them, and quickly spin them down to collect the liquid at the bottom.
7. Place the plates in a thermocycler and run the following protocol:

| Step | Temperature (C) | Time (min)  |
| ---- | ------------- | ------------- |
| 1 | 37 C | 5:00 |
| 2 | 16 C | 5:00 |
| 3 | GOTO 1 | 29 times |
| 4 | 37 C | 10:00 |  
| 5 | 4 C | Indefinitely |

8. Replace the permeable seal from the synthesis plates with an aluminum seal and store in a -20C freezer

### Code Breakdown
1. Asks the user how many reactions they would like to build as well as which enzyme to use
2. Queries the database for all of the parts that are buildable
  - A single build is currently limited to 96 samples
  - Currently preferences parts whose fragments are in the earliest synthesis plates added to the database
3. Determines the build number of the current build
4. Generates a new build object in the database and added a new build plate and placed all of the parts in their corresponding wells
  - Also generates a csv file detailing which gene was in each well
5. Iterates through all of the genes to be built and pulls all of the synthesis plates required to carry out the build
6. Groups the plates by three as that is the number of available slots on the OpenTrons given the current plate map
7. Checks each plate to see if it has already been resuspended. If the plate has not, it will enter into the plate resuspension protocol.
8. Once all of the plates are resuspended, it will generate and display the deck layout
9. Calculates the number of reactions needed to run the build
  - Currently uses an extra reaction volume for each fragment that is needed for a certain gene to maintain the proportion of DNA to master mix
10. Calculates and presents the volume of each reagent needed to generate the master mix
11. Once ready it will initialize the robot, all of the containers and the pipettes
12. Aliquots the master mix into PCR tube strip with the p200
13. Adds the 8µL of master mix for each reaction in the PCR tube plates using the p10 multichannel
14. If there is an unfilled row, then it uses the p10 single channel to add master mix to each of remaining wells
15. Adds an extra reaction volume of master mix for each extra fragment in each well
16. Reorganizes the build plan such that it will pipette all of the needed wells from a single synthesis plate before moving to the next
  - Minimizes the need to swap out plates, in the event that more plates are needed than can fit on the deck at one time
17. Iterates through each row of the build plan and transfers the fragments to their target wells
  - If the recorded volume is lower than the specified amount then it will add an extra 5 µL of water to the well prior to pipetting
  - If the plate has not been pipetted from yet then it will enter into a manual calibration step where you can drive the motor up or down to make sure that it is set to the right height


## Transformation

### Purpose
Transform chemically competent *E. coli* cells with the golden gate reactions

### Requirements
- PCR Tube Plate(s) with competent cells in as many tubes as needed
  - We use 7.5 µL of competent cells in each tube
  - We make competent Top10 *E. coli* cells using the Zymo Mix and Go kit detailed in the material list
- DNA to transform
- Container of ice

### User Intervention
1. Thaw the plate of competent cells on ice for 5-10 minutes
2. Once all of the tubes have thawed, label it such that the 'A1' well is clearly defined and add it to the spot labelled `Transformation`
3. Add the Golden Gate reactions to the slot labelled `Build_plate`
4. After the DNA is transferred, move the transformation plate to the container of ice
5. After waiting 30 minutes, heat shock the cells in a 42C water bath for 30 seconds
6. Proceed directly to the plating protocol

### Code Breakdown
1. Asks the user which build they would like to transform
2. Initializes the robot, all of the containers and the pipettes
3. Uses the p10 multichannel to transfer the DNA from the build plate to the corresponding row in transformation plate
  - If there is an incomplete row of DNA to transform, it will use the p10 single channel to add the remaining constructs


## Serial Dilution & Plating

### Purpose
Dilute and plate the cells to yield pickable transformants for each golden gate assembly

### Requirements
- 4x agar plates with the appropriate antibiotic resistance
  - We use rectangular petri plates with 40 mL of LB agar to maintain a consistent agar thickness
  - Dry the plates out in a 37C incubator prior to plating to make sure that the droplets can soak into the media
- Transformed cells
- 1x 8-tube PCR strip
- 2x microfuge tubes with 1.2 mL of LB

### User Intervention
1. Place the 8-tube PCR strip in the bottom row of the rack in the `PCR-strip-tall` slot
2. Place the transformation tube in the slot labelled `Transformation`
3. Place the two tubes of LB in the 'A1' and 'B1' wells in the rack in the slot labelled `Tube_rack`
4. When prompted, enter the build that you would like to plate
5. Write the labels presented by the script on the agar plates and place them on the deck where specified
6. When transferring cells to a new plate, the script will enter into a calibration step to ensure that the droplets containing the cells make contact with the agar
  - If you change the calibration at all it will enter into the calibration step again to double check that the new calibration is good. Once no change occurs during a calibration, it will no longer ask to recalibrate and will apply the new z-coordinate to all of the proceeding rows.
7. After the plating is complete, allow the plates to sit facing up and covered on the bench for at least 15 minutes to allow the droplets to fully soak into the agar.
8. Move the plates to 37C incubator and place them upside down

### Code Breakdown
1. Queries the database to find which builds have not yet been plated and presents the results to the user.
2. Asks the user which build they would like to plate.
3. Pulls in the information about the build
  - Currently, if the build consists of more than 48 constructs it will ask the user which half they would like to plate as the deck layout is unable to hold 4 agar plates at once
4. Determines the number of plates needed and assigns them names
  - The names are assigned as '[build_number]\_p[plate_number]' where the plate_number is assigned base on which fraction of the wells are plated on each (i.e. wells 'A4'-'H7' would be on plate 'build000_p2')
5. The agar plates are then assigned to the open positions on the deck.
6. Initialize the robot, containers and pipettes
7. Uses the p200 single channel to transfer LB to the PCR tube strip
8. Iterates through a dilution and plating cycle consisting of the following steps:
  1. Dilutes the cells by using the p10 multichannel to aliquot LB into the target row
  2. Mix the LB with the cells
  3. Transfer the cells to the agar plate
    - Dispenses the cells prior to reaching the agar and then stabs the tips slightly into the agar such that the droplets that pull up on the side of the tips make contact with the agar
      - This avoids the need for perfect calibration of the z-axis
    - When initially plating on a plate it will enter into a z-axis calibration step to ensure that all of the droplets make contact with the plate.
  4. Pipettes off extra cells from each tube to allow for sufficient dilution factors
  5. Repeats the process for the number of dilutions desired for each transformation
    - We have found that 4 serial dilutions leads to an effective dilution series with the volume of competent cells and Golden Gate reaction we add. Might need to be tuned depending on cloning efficiency
9. Once a plate is completed, it will aliquot out more LB into the PCR tube strip

## Colony Picking

### Purpose
Inoculating cultures from colonies to propagate the clones in order to sequence verify the constructs

### Requirements
- 1x deep well plate
  - We use a 2.2 mL deep well plate but a lower volume plate should work too
- Permeable seal
- Agar plates with matured colonies
- LB with the desired antibiotic
- Calibrated USB webcam
  - Calibrate the camera using the chessboard calibration system
    - [ ] Clean up the calibration scripts so that it is very easy to do the initial calibration
- Image box
- Transilluminator

### User Intervention
1. Retrieve the agar plates from the incubator
2. Run the `capture_build.py` script
3. Specify the build that you would like to plate
4. Place the plates inside the imaging box with the transilluminator on, as specified.
  - The plate should be placed face up, with the lid off, oriented with the 'A1' or equivalent well on the bottom right
  - Press `q` while the camera is streaming to capture the current frame
5. Once all of the images are taken, aliquot LB+antibiotic into a deep well plate for as wells as there are transformants.
  - We use a p1250 multichannel to quickly fill a 96 deep well plate with 1250 µL per well as it is far faster to do that then to have the OpenTrons pipet it with the current pipettes that are mounted
  - We pour about 60 mLs in a sterile disposable trough to aliquot and then refill it.
  - We use a final volume of 1.25 mL in each well so that we can use 200 µL to make 2x glycerol stocks and have a mL left to prep/sequence
6. Run the `colony_picking.py` script
7. Specify the build that you would like to pick from
- [ ] Clean up the colony picking code to display the deck layout
- [ ] Modify the colony picking script such that it allows the user to do all of the calibrations at the start
8. Place the deep well plate in slot 'C2'
9. Place the transformation plate in slot 'B2'
10. Follow the prompts to calibrate the system. This includes:
  - Confirming the location of the agar edge
  - Calibrating the xy coordinates for the four outer most colonies
  - Calibrating the z coordinate for a single colony
11. Once the protocol is finished, seal the deep well plate with a permeable seal
12. Place the deep well plate in the 37C plate shaking incubator and grow overnight

### Code Breakdown

#### Image Capture
1. In the first script `capture_build.py`, it will display and ask which build the user would like to photograph.
2. Once chosen, it will create a new directory to contain the newly captured images
3. It then analyzes the build, determines which plates need to be photographed and prompts the user to place a specific plate into the image box
4. The script will then stream the webcam within the image box
5. It will continue to stream until the user presses `q` which will save the image to the new directory using the name of the plate
5. The script the repeats this cycle for all of the plates in the build

#### Image Analysis
1. The second script, `colony_picking.py`, begins by asking the user which build they would like to pick colonies from
2. Once the build is specified, it will open the images generated by  `capture_build.py`
  - All further steps occur on each image individually
3. The plate name is extracted from the file name and used to identify which wells the plate contains
4. Undistorts the image to account for the radial and tangential distortions intrinsic to the camera
  - Utilizes a previously determined set of distortion coefficients intrinsic to the camera
  - See the camera calibration page for further details on determining these coefficients
5. Rotates the image so that it reflects the orientation of the plate when placed on the OpenTrons deck
6. Finds the edges of the plate by setting a low threshold to pull out the entire plate from the background and then approximate the shape with a polygon
7. It then uses the corners of the approximated shape outlining the plate to apply a geometric transform to justify and crop the image to only include the plate
8. From the modified image it guesses where the edge of the agar are within the plate and asks the user to correct them if needed
  - The edge of the agar appears as a subtle feature in the image and is variable across images so currently it assumes that the corners will always be in the same place and so they are hard coded in the script
  - As the locations of these points are used to map the pixels to the OT coordinates, it is important that they are accurate
9. Once found, it will apply another transform that crops the image to just the agar slab
10. It then segments the image based on clusters of colonies on the plate to generate a grid
  - The script is currently set to have an 12 rows and 8 columns like a standard 96 well plate
  - It does so by calculating the sum of pixels in a given row/column of pixels and finds all of the low points which should represent the gaps between the clusters of colonies
11. Using the intersection points in the grid, it then groups the individual sections with other sections corresponding to the same construct
  - Currently we use 4 serial dilutions
12. The groups are then reordered such that they reflect the order that they are plated
  - The cells are plates in a left to right, bottom to top manner, i.e. 'A1' would be the bottom left section and 'H3' would be the top right section
13. Assigns each section of the grid the corresponding well in the target deep well plate, based on which plate is currently being analyzed and where the section is located on the plate
14. Looks at each section and chooses the best candidate colony available
  - Uses an adaptive thresholding system to isolate the colonies from the background
  - Analyzes several different attributes of each colony and checks them against empirically determined criteria to choose the best colony
  - If no colonies are found with the stringent criteria, it will broaden the acceptable ranges of each attribute to find a reasonable colony
  - If that also fails, it will save that there were no pickable colonies and ask the user to deal with it
  - Passing the parameter `-c` when running the script will allow the user to go through and choose each individual colony
15. Calculates the center of mass of all of the desired colonies

#### Colony Picking
1. Takes the centers of all of the colonies and converts them to xy coordinates that the OpenTrons can use
  - Maps the coordinates by generating a proportion of the measured dimensions of the petri plate and the size of the image in pixels
  - The reference point on the robot is the bottom left corner and so the coordinates are modified to reflect the distance in both axis from that point
2. Then calibrates the x and y axis to ensure that the robot will accurately pick all of the colonies
  - The size of a colony can be less than 1mm making the tolerance low which makes it important to be very accurate
  - It will grab a tip and then hover over the bottom left colony
  - Once over where it thinks the colony is, the user is able to drive the motor around to account for any offset.
  - This process repeats for all colonies closest to each corner to store the x and y offsets across the plate
  - It then applies a bilinear interpolation to modify all of the coordinates to account for the offset
    - It will also present a heat map displaying the offsets in both dimensions to visualize the corrections.
3. The robot will then return to the first colony and allow the user to calibrate the Z axis
4. Once the pipette tip has sufficient contact with the colony it will transfer the tip to the desired well in the deep well plate, mix the tip around and then drop it in the well to ensure that the colony inoculates the well
5. The colony picking process then continues for all of the remaining colonies


## Sequence Alignment
- [ ] Clean up and generalize this code a lot. It still uses the 10K_CDS csv to account for the very early sequences
- [ ] convert the backbone sequence to a fasta file in a different directory

### Purpose
Assign an outcome to each part that was attempted in a build based on the sequencing results.

### Requirements
- Sanger sequencing reads
- Backbone sequence (Truncated to the region that should be replaced by the insert)

### User Intervention
1. Transfer the sequencing files to a directory labelled '[build number]\_seq\_files' within the corresponding build directory
2. Check to make sure that the beckbone sequence, and names of the sequencing primers are correct
3. Specify which build you would like to analyze
4. Allow the alignments to run
  - Can take a few minutes if aligning hundreds of sequences

### Code Breakdown
- [ ] Add the code breakdown

















* space
