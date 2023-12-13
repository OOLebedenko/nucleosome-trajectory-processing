## 1. change the path to SCRIPT_DIR according to your directory tree
SCRIPT_DIR:={YOUR PATH}/nucleosome-trajectory-processing/15N_relaxation_rates/scripts/

## 2. specify MD trajectory parameters
# set trajectory path; for example, md_setup/md_protocol/TIP4P-D_disp/ff99SB-disp/wt/02_production/6_run/
TRAJECTORY_PATH:=""
# set path to reference pdb (of note, chains should be labeled same as in 3LZ0); for example, md_protocol/TIP4P-D_disp/ff99SB-disp/wt/01_equil_histone_tails/0_prepare/wt.pqr
REFERENCE_PDB_PATH:=""
# set trajectory length to be processed
TRAJECTORY_LENGTH:="" # ns
# set type of trajectory files: "dat" - TrjtoolDatFile; "nc" - AmberNetCDF; "xtc" - GromacsXtcFile
FILETYPE:=nc
# set pattern for trajectory filenames: run00001.dat ---> "run%05d".dat
PATTERN:=run%05d
# set timestep used to record trajectory files
DT_NS=0.01 # ns


## 3. set path to Xray reference pdb relative to which RMSD is calculated
XRAY_REF="" # for example 3lz0
