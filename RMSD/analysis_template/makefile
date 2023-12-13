include common.mk

all: data/all.csv data/rmsd_protein.csv data/rmsd_dna.csv data/rmsd_dna_inner.csv data/rmsd_dna_outer.csv

clean:
	rm -rf data

data/all.csv data/rmsd_protein.csv data/rmsd_dna.csv data/rmsd_dna_inner.csv data/rmsd_dna_outer.csv: $(strip ${SCRIPT_DIR})/calc_rmsd.py
	python "$<" \
        --path-to-trajectory="$(strip ${TRAJECTORY_PATH})" \
        --path-to-reference-pdb="$(strip ${REFERENCE_PDB_PATH})" \
        --path-to-xray-reference-pdb="$(strip ${XRAY_REF})" \
        --filetype="$(strip ${FILETYPE})" \
        --pattern="$(strip ${PATTERN})" \
        --trajectory-length="$(strip ${TRAJECTORY_LENGTH})" \
        --dt-ns="$(strip ${DT_NS})" \
        --output-directory=data/