# import os
import json
import numpy as np
from pathlib import Path

# //////////////////////////////////////////////////////////////////////////////

####################################### PARAMETERS

# SIMDIR = Path(os.getcwd()).parent
SIMDIR = Path("D:/CulebraExt/ComputationalBiophysicsSimulations/project")

##### General
SYSTEM_NAME = "md"
TRAJECTORY_NAME = "md-MDA"
RUNS = ["8DFN_full", "8DFN_half", "7SI9_full", "7SI9_half"]

##### Wetness Analysis (WA)
WA_CURRENT_RUN = "8DFN_full"
WA_NAME = "test"

WA_DIR = SIMDIR / WA_CURRENT_RUN / "wetness_analysis"
WA_POSTDIR = SIMDIR / WA_CURRENT_RUN / "5_POST"
WA_INFO_PATH = WA_DIR / (WA_NAME + "-info.json")

WA_DIVISIONS = 20
DENSITY_TRESHOLD_LOWER = 0
DENSITY_TRESHOLD_UPPER = 3


####################################### OTHER CONSTANTS

RED = (1, 0, 0)
GREEN = (0, 1, 0)

####################################### FUNCTIONS

def get_POSTDIR(run): return SIMDIR / run / "5_POST"

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
            print(f"--- '{self.path.name}' not found.")



    def __getattr__(self, key):
        return self.__dict__.get(key, self.info[key])

    def write(self):
        with open(self.path, 'w') as file: file.write(json.dumps(self.info))

    def update(self, **data):
        self.info = {**self.info, **data}
        self.write()

        print(f"--- Updating '{self.path.name}' with {len(data)} parameters...")

# //////////////////////////////////////////////////////////////////////////////
