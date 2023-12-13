## 1. change the path to SCRIPT_DIR according to your directory tree
SCRIPT_DIR:=../../scripts/

## 2. specify MD traj parameters
# set trajectory path; for example md_setup/md_protocol/TIP4P-D_disp-replica-1/ff99SB-disp/wt/02_production/6_run/
TRAJECTORY_PATH:=../../../short_example_trajectory/trajectory_GromacsXtc/
# set path to reference pdb (of note, we need ); for example md_protocol/TIP4P-D_disp-replica-1/ff99SB-disp/wt/01_equil_histone_tails/0_prepare/wt.pqr
REFERENCE_PDB_PATH:=../../../short_example_trajectory/reference.pdb
# set trajectory length
TRAJECTORY_LENGTH:=10 # ns
# set type of trajectory files "dat" - TrjtoolDatFile; "nc" - AmberNetCDF, "xtc" - GromacsXtcFile
FILETYPE:=xtc
# set pattern of trajectory files files: run00001.dat ---> "run%05d".dat
PATTERN:=run%05d
# set step of printing out coordinates in trajectory
DT_NS=0.01 # time step of output for frame coordinates


## 3. set path to Xray reference pdb relative to which RMSD is calculated
XRAY_REF=../../../short_example_trajectory/3lz0.pdb