import argparse
import functools
import os
import operator

from tqdm import tqdm
from pyxmolpp2 import PdbFile, Trajectory, TrjtoolDatFile, AmberNetCDF, GromacsXtcFile, mName, rId, aName, UnitCell, Degrees
from pyxmolpp2.pipe import AssembleQuaternaryStructure, Align, Run

from process_utils.select import get_sec_str_residues_predicate
from process_utils.calc import CalcRmsd


class XtcFileReaderWrapper:
    def __init__(self, n_frames):
        self.n_frames = n_frames

    def __call__(self, filename):
        return GromacsXtcFile(filename, self.n_frames)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract vectors')
    parser.add_argument('--path-to-trajectory', required=True)
    parser.add_argument('--path-to-reference-pdb', required=True)
    parser.add_argument('--path-to-xray-reference-pdb', required=True)
    parser.add_argument('--filetype', choices=["dat", "nc", "xtc"], required=True)
    parser.add_argument('--pattern', default="run%05d")
    parser.add_argument('--trajectory-start', default=1, type=int)
    parser.add_argument('--trajectory-length', required=True, type=int)
    parser.add_argument('--trajectory-stride',default=1, type=int)
    parser.add_argument('--frames-per-trajectory-file', type=int, default=100)
    parser.add_argument('--dt-ns', type=float, default=0.01)
    parser.add_argument('--output-directory', default=".")
    args = parser.parse_args()

    trj_reader_dict = {"dat": TrjtoolDatFile,
                       "nc": AmberNetCDF,
                       "xtc": XtcFileReaderWrapper(args.frames_per_trajectory_file)}

    #  load trajectory
    reference = PdbFile(args.path_to_reference_pdb).frames()[0]
    traj = Trajectory(reference)
    for ind in tqdm(range(args.trajectory_start, args.trajectory_length + 1), desc="traj_reading"):
        fname = "{pattern}.{filetype}".format(pattern=args.pattern, filetype=args.filetype)
        traj.extend(trj_reader_dict[args.filetype](os.path.join(args.path_to_trajectory, fname % (ind))))

    # set xray reference to calc RMSD
    xray_ref = PdbFile(args.path_to_xray_reference_pdb).frames()[0]

    # select secondary-strucuture residues
    protein_chains = ["A", "B", "C", "D", "E", "F", "G", "H"]
    ss_residues = get_sec_str_residues_predicate(frame=reference, molnames=protein_chains)
    ss_ca_atoms = ss_residues & (aName == "CA")

    # select inner and outer DNA turns
    dna_chains = ["I", "J"]
    dna_align_pred = aName.is_in("N1", "N9") & mName.is_in(set(dna_chains))
    dna_inner_turn = set(range(-38, 38))
    dna_outer_turn = set(range(-72, -39)) | set(range(39, 72))
    inner_dna_seceletion = ((mName.is_in(set(dna_chains)) & rId.is_in(dna_inner_turn))) & aName.is_in("N1", "N9")
    outer_dna_seceletion = (mName.is_in(set(dna_chains))) & (rId.is_in(dna_outer_turn)) & aName.is_in("N1", "N9")

    # process trajectory
    tqdm(traj[::args.trajectory_stride]) \
    | AssembleQuaternaryStructure(of=(mName.is_in("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")),
                                  by=aName.is_in("P") | ss_ca_atoms,
                                  reference=reference) \
    | Align(by=ss_ca_atoms, reference=xray_ref) \
    | CalcRmsd(reference=xray_ref,
               by_atoms=ss_ca_atoms | dna_align_pred,
               dt_ns=args.dt_ns,
               out_filename="all.csv",
               out_dirname=args.output_directory) \
    | CalcRmsd(reference=xray_ref,
               by_atoms=ss_ca_atoms,
               dt_ns=args.dt_ns,
               out_filename="rmsd_protein.csv",
               out_dirname=args.output_directory) \
    | CalcRmsd(reference=xray_ref,
               by_atoms=dna_align_pred,
               dt_ns=args.dt_ns,
               out_filename="rmsd_dna.csv",
               out_dirname=args.output_directory) \
    | CalcRmsd(reference=xray_ref,
               by_atoms=inner_dna_seceletion,
               dt_ns=args.dt_ns,
               out_filename="rmsd_dna_inner.csv",
               out_dirname=args.output_directory) \
    | CalcRmsd(reference=xray_ref,
               by_atoms=outer_dna_seceletion,
               dt_ns=args.dt_ns,
               out_filename="rmsd_dna_outer.csv",
               out_dirname=args.output_directory) \
    | Run()
