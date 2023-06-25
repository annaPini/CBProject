from pathlib import Path

# //////////////////////////////////////////////////////////////////////////////
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


# //////////////////////////////////////////////////////////////////////////////
