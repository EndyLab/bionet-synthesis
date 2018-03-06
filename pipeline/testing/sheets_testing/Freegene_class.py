import json
import yaml
import os
import glob
import re
from collections import namedtuple

class Freegene:
     def __init__(self, gene_id):
         gene_file = gene_id + ".json"
         with open(gene_file,"r") as json_file:
             self.data = json.load(json_file, object_hook=lambda d: namedtuple('X', d.keys())(*d.values()))

def synfragment(self):
    # Open configuration
    with open("fragment_config.yaml","r") as yaml_file:
        yaml_config = yaml.load(yaml_file)
        # Set global variables
        max_length = yaml_config["Synthesis_configuration"]["max_length"]
        min_length = yaml_config["Synthesis_configuration"]["min_length"]
        standard_flank_prefix = yaml_config["Standard_flanks"]["prefix"]
        standard_flank_suffix = yaml_config["Standard_flanks"]["suffix"]
        sequence = self.data.sequence.optimized_sequence




