import json
from pathlib import Path

# ////////////////////////////////////////////////////////////////////////////// WAD PARAMETERS
RUN_WAD = "mt2_rep0"
WAD_DIVISIONS = 20
WAD_DENSITY_TRESHOLD_LOWER = 0
WAD_DENSITY_TRESHOLD_UPPER = 3

###########################################################
DIR_DA = Path.cwd().parent / "data_analysis"
DIR_DA_TRAJECTORIES = DIR_DA / "_trajectories"
DIR_DA_WAD = DIR_DA / "wad"

WAD_NAME = RUN_WAD
PATH_WAD_INFO = DIR_DA_WAD / f"{WAD_NAME}-info.json"

# ////////////////////////////////////////////////////////////////////////////// CONTAINER CLASS
class Info:
    def __init__(self, info_path):
        self.path = Path(info_path)
        if self.path.exists():
            with open(self.path, 'r') as file:
                self.info = json.loads(file.read())
            print(f"--- Loaded {len(self.info)} parameters from '{self.path.name}'.")
        else:
            print(f"--- '{self.path.name}' not found, will create new one.")
            self.info = {}


    def __getattr__(self, key):
        return self.__dict__.get(key, self.info[key])

    def update(self, **data):
        self.info = {**self.info, **data}
        with open(self.path, 'w') as file:
            file.write(json.dumps(self.info))

        print(f"--- Updated '{self.path.name}' with {len(data)} parameters.")

# //////////////////////////////////////////////////////////////////////////////
