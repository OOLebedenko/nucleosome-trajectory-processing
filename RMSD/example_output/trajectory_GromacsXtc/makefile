include common.mk

all: figure/rmsd.png

clean:
	rm -rf data
	rm -rf figure

figure/rmsd.png: $(strip ${SCRIPT_DIR})/plot_rmsd.py data/rmsd_*.csv
	python "$<" \
        --path-to-rmsd_dir=data/ \
        --output-directory=figure/

data/rmsd_*.csv: $(strip ${SCRIPT_DIR})/calc_rmsd.py
	python "$<" \
        --path-to-trajectory="$(strip ${TRAJECTORY_PATH})" \
        --path-to-reference-pdb="$(strip ${REFERENCE_PDB_PATH})" \
        --trajectory-stride="$(strip ${TRAJECTORY_STRIDE})" \
        --path-to-xray-reference-pdb="$(strip ${XRAY_REF})" \
        --filetype="$(strip ${FILETYPE})" \
        --pattern="$(strip ${PATTERN})" \
        --trajectory-length="$(strip ${TRAJECTORY_LENGTH})" \
        --dt-ns="$(strip ${DT_NS})" \
        --output-directory=data/