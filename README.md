# OpenFoundry

The OpenFoundry is a collection of automation and inventory management tools that enable high-throughput cloning of synthesized DNA. The system utilizes OpenTrons robotics which are inexpensive and open source liquid handling robots. The pipeline automates every step in the cloning process and tracks the progress of each sample. In our hands, we have been able to maintain a rate of cloning 500 genes per robot per week.

## Getting Started
Clone the OpenFoundry repository
`git clone https://github.com/EndyLab/bionet-synthesis.git`

Install all of the requirements
`cd bionet-synthesis`
`pip install -r requirements.txt`

## Workflow Overview
1. Use the template.csv file as a guide to generate a csv in the `./src/` directory with the sequences that you would like to synthesize and clone.
2. Optimize all of the sequences present using `optimize.py`.
3. When ready to order, run `fragment.py` to fragment the genes and append the desired overhangs and send out for synthesis.
4. When the plates come in, run `frag_location_assignment.py` to parse the newly added plate maps and add locations to all of the fragments. This will also determine which genes are now able to be built.
5. Run `resuspend.py` on each synthesis plate that you would like to resuspend.
6. When the necessary fragments arrive, use the `create_assembly_plan.py` script to generate an assembly protocol that coordinates the cloning across all available robots.
7. Once the plan is created, run `assemble.py` on each of the robots, selecting the desired builds for each.
8. After the assembly process is complete, place them in the thermocycler and run the designated program.
9. Once the thermocycle is complete, run `transform.py`.
10. After the transformation run the `plating.py` script.
11. Place plates in the 37C incubator overnight.
12. Retrieve the plates and photograph them within the image box using `image_plates.py`.
13. Fill up a 96 deep well plate with desired media and execute. `colony_picking.py` to inoculate the well plate with the desired colonies.
14. Place the deep well block in a plate shaker at 37C overnight.
15. Either run a 96 well miniprep or send the whole block out for sample preparation and sequencing.
16. If using Sanger sequencing, add the sequencing results to the corresponding build and use `seq_alignment.py` to generate sequence alignments of all the corresponding target sequences and determine the outcomes of the cloning reactions.

Refer to `documentation.md` for extensive protocols and explanation of code.
See the `openfoundry_material_list.csv` for information on the materials used throughout these protocols

## License
MIT License, refer to `LICENSE.md` for more information










#
