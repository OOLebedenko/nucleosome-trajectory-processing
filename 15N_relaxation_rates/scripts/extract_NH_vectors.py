from process_utils.select import get_sec_str_residues_predicate, atom_pairs_from_one_residue
from pyxmolpp2 import PdbFile, Trajectory, TrjtoolDatFile, AmberNetCDF, GromacsXtcFile, Atom, mName, aName
from pyxmolpp2.pipe import AssembleQuaternaryStructure, WriteVectorsToCsv, Run, Align
from tqdm import tqdm
import os
import argparse


class XtcFileReaderWrapper:
    def __init__(self, n_frames):
        self.n_frames = n_frames

    def __call__(self, filename):
        return GromacsXtcFile(filename, self.n_frames)


class OutputFilenameFormatter:
    def __init__(self, output_directory):
        self.output_directory = output_directory

    def __call__(self, atom1: Atom, atom2: Atom):
        output_directory = os.path.join(self.output_directory)
        os.makedirs(output_directory, exist_ok=True)
        return f"{output_directory}/" \
               f"{atom1.residue.id.serial:02d}_{atom1.name}" \
               f"{atom2.name}" \
               f".csv"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract vectors')
    parser.add_argument('--path-to-trajectory', required=True)
    parser.add_argument('--path-to-reference-pdb', required=True)
    parser.add_argument('--chain-name', required=True)
    parser.add_argument('--residue-of-interest', required=True)
    parser.add_argument('--filetype', choices=["dat", "nc", "xtc"], required=True)
    parser.add_argument('--pattern', default="run%05d")
    parser.add_argument('--trajectory-start', default=1, type=int)
    parser.add_argument('--trajectory-length', required=True, type=int)
    parser.add_argument('--frames-per-trajectory-file', type=int, default=100)
    parser.add_argument('--output-directory', default=".")
    args = parser.parse_args()

    trj_reader_dict = {"dat": TrjtoolDatFile,
                       "nc": AmberNetCDF,
                       "xtc": XtcFileReaderWrapper(args.frames_per_trajectory_file)}

    #  load trajectory
    reference = PdbFile(args.path_to_reference_pdb).frames()[0]
    traj = Trajectory(PdbFile(args.path_to_reference_pdb).frames()[0])
    for ind in tqdm(range(args.trajectory_start, args.trajectory_length + 1), desc="traj_reading"):
        fname = "{pattern}.{filetype}".format(pattern=args.pattern, filetype=args.filetype)
        traj.extend(trj_reader_dict[args.filetype](os.path.join(args.path_to_trajectory, fname % (ind))))

    # select secondary-strucuture residues
    protein_chains = ["A", "B", "C", "D", "E", "F", "G", "H"]
    ss_residues = get_sec_str_residues_predicate(frame=reference, molnames=["A", "B", "C", "D", "E", "F", "G", "H"])

    # set residues of interest
    first_rid, last_rid = args.residue_of_interest.split("-")
    residues_of_interest = set(list(range(int(first_rid), int(last_rid) + 1)))

    # process trajectory
    tqdm(traj) \
    | AssembleQuaternaryStructure(of=(mName.is_in("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")),
                                  by=aName.is_in({"P", "CA"}),
                                  reference=reference) \
    | Align(by=(aName == "CA") & ss_residues, reference=reference) \
    | WriteVectorsToCsv(pair_selector=atom_pairs_from_one_residue(args.chain_name, "N", "H", residues_of_interest),
                        filename_provider=OutputFilenameFormatter(args.output_directory)) \
    | Run()
