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


## 3. specify fit parameters
FIT_LIMIT_NS=1000 # ns
# you may specify logarithmic resampling of the correlation function (LAG_INDEX="log") and set the number of points (N_LAG_POINTS) to cover the range from 0 ns to FIT_LIMIT_NS.
# this choice may prevent overfitting of the correlation function at long timescales
# if you want to fit the correlation function without the logarithmic resampling, you should delete these parameters (LAG_SPACING and N_LAG_POINTS)
LAG_SPACING="log"
N_LAG_POINTS=1000


## 4. specify experimental parameters needed for calculations
NMR_FREQ=850e6 # Hz
TUMBLING_TIME=163.4 # ns