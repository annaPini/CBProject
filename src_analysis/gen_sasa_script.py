from parameters import *

if __name__ == "__main__":
    PATH_GRO_MT2_0     = DIR_DA_TRAJECTORIES / f"mt2_rep0.gro"
    PATH_XTC_MT2_0     = DIR_DA_TRAJECTORIES / f"mt2_rep0.xtc"
    PATH_TPR_MT2_0     = DIR_DA_TRAJECTORIES / f"mt2_rep0.tpr"
    PATH_AREA_MT2_0     = DIR_DA_SASA / f"area.xvg"
    PATH_SFE_MT2_0     = DIR_DA_SASA / f"sfe.xvg"  # solvation free energy as a function of time
    PATH_OR_MT2_0     = DIR_DA_SASA / f"or.xvg"  # average area per residue
    PATH_OA_MT2_0     = DIR_DA_SASA / f"oa.xvg"  # average area per atom
    PATH_TV_MT2_0     = DIR_DA_SASA / f"volume.xvg"  # total volume and density as a function of time

    #gmx sasa -f PATH_XTC_MT2_0 -s PATH_TPR_MT2_0 -o PATH_AREA_MT2_0 -odg PATH_SFE_MT2_0 -or PATH_OR_MT2_0 -oa PATH_OA_MT2_0 -tv PATH_TV_MT2_0 -surface protein

   # gmx sasa -f mt2_rep0.xtc -s mt2_rep0.tpr -o area.xvg -odg sfe.xvg -or or.xvg -oa oa.xvg -tv volume.xvg -surface protein

    PATH_SASA = "sasa.sh"

    command = f"gmx sasa -f {PATH_XTC_MT2_0} -s {PATH_TPR_MT2_0} -o {PATH_AREA_MT2_0} -odg {PATH_SFE_MT2_0} -or  {PATH_OR_MT2_0} -oa {PATH_OA_MT2_0} -tv {PATH_TV_MT2_0} -surface protein"

    with open(PATH_SASA, 'w') as file: file.write(command.replace("\\", "/"))