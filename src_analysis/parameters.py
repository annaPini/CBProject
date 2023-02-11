import os,  json
import numpy as np
from pathlib import Path

# //////////////////////////////////////////////////////////////////////////////

####################################### PARAMETERS

##### General
DIR_PROJECT = Path(os.getcwd()).parent
RUN_PREFFIXES = ["mt1", "mt2", "wt1", "wt2"]
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
DIR_DEXTRA = DIR_PROJECT / "data_extra"

DIR_DA_TRAJECTORIES = DIR_DA / "trajectories"
DIR_DA_GENERAL = DIR_DA / "general"
DIR_DA_PYINTERAPH = DIR_DA / "pyinteraph"
DIR_DA_SPECIFIC = DIR_DA / "specific"
DIR_DA_VMD = DIR_DA / "vmd"
DIR_DA_WAD = DIR_DA / "WAD"
DIR_DA_SASA = DIR_DA / "sasa"

####################################### CONSTANTS

WAD_NAME = f"{CURRENT_RUN}-{WAD_REF_FRAME}"
PATH_WAD_INFO = DIR_DA_WAD / f"{WAD_NAME}-info.json"

####################################### OTHER CONSTANTS

RED   = 1, 0, 0, 1
GREEN = 0, 1, 0, 1
BLUE  = 0, 0, 1, 1

HALF_RED   = 1, 0, 0, 0.5
HALF_GREEN = 0, 1, 0, 0.5
HALF_BLUE  = 0, 0, 1, 0.5

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
