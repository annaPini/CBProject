import os,  json
import numpy as np
from pathlib import Path

# //////////////////////////////////////////////////////////////////////////////

####################################### PARAMETERS

##### General
DIR_PROJECT = Path(os.getcwd()).parent
RUNS = ["mt1_rep0", "mt1_rep1", "mt2_rep0", "mt2_rep1", "wt1_rep0", "wt1_rep1", "wt2_rep0", "wt2_rep1"]
CURRENT_RUN = "mt2_rep0"

##### Wetness Average Density (WAD)
WAD_REF_FRAME = "0"
WAD_DIVISIONS = 20
WAD_DENSITY_TRESHOLD_LOWER = 0
WAD_DENSITY_TRESHOLD_UPPER = 3

####################################### DIRECTORIES

DIR_DMD = DIR_PROJECT / "data_md"
DIR_DA = DIR_PROJECT / "data_analysis"

DIR_DA_TRAJECTORIES = DIR_DA / "trajectories"
DIR_DA_GENERAL = DIR_DA / "general"
DIR_DA_PYINTERAPH = DIR_DA / "pyinteraph"
DIR_DA_SPECIFIC = DIR_DA / "specific"
DIR_DA_VMD = DIR_DA / "vmd"
DIR_DA_WAD = DIR_DA / "WAD"

####################################### CONSTANTS

WAD_NAME = f"{CURRENT_RUN}-{WAD_REF_FRAME}"
PATH_WAD_INFO = DIR_DA_WAD / f"{WAD_NAME}-info.json"

####################################### OTHER CONSTANTS

RED   = 1, 0, 0, 1
GREEN = 0, 1, 0, 1
BLUE  = 0, 0, 1, 1

####################################### CLASSES

class Info:
    def __init__(self, info_path):
        self.path = Path(info_path)
        if self.path.exists():
            with open(self.path, 'r') as file:
                self.info = json.loads(file.read())
            print(f"--- Loaded {len(self.info)} parameters from '{self.path.name}'.")
        else:
            self.info = {}
            print(f"--- '{self.path.name}' not found, will create new one.")

    def __getattr__(self, key):
        return self.__dict__.get(key, self.info[key])

    def update(self, **data):
        self.info = {**self.info, **data}
        with open(self.path, 'w') as file:
            file.write(json.dumps(self.info))

        print(f"--- Updated '{self.path.name}' with {len(data)} parameters.")

# //////////////////////////////////////////////////////////////////////////////


######### check "solvent accesible surface area"

# extensions --> volmap tool

# 8f - 7f
# 8h - 7h
# 8f - 8h*2
# 7f - 7h*2

# source pyinteraph_venv/bin/activate
# cd /media/sf_CulebraBox/NDProjects/temp_project/8DFN_full/pyinteraph
# pyinteraph -s md_plain.tpr -t md-rottrans.xtc -r md.gro --sb-co 5 -b --sb-graph sb-graph.dat --ff-masses charmm27 -v --sb-cg-file charged_groups.ini

# wt = 7SI9
# mt = 8DFN
