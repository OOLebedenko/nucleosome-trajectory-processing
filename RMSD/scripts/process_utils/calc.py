from pyxmolpp2.pipe import TrajectoryProcessor
from pyxmolpp2 import Frame, AtomPredicate, calc_rmsd
import os
import csv


class CalcRmsd(TrajectoryProcessor):

    def __init__(self,
                 reference: Frame,
                 by_atoms: AtomPredicate,
                 dt_ns: float,
                 out_filename: str = "rmsd.csv",
                 out_dirname: str = ".") -> None:
        """
        The idea is to run of the processing trajectory in pipe-like format:
        trajectory | CalcRmsd(dt_ns=0.001) | Run()

        :param reference: reference structure
        :param by: selector of atoms which are used to calculate RMSD
        :param dt_ns: time step between frame in trajectory
        :param out_filename: by default cm.csv
        :param out_dirname: by default the current directory
        """

        self.reference = reference
        self.atoms_selector_for_alignment = by_atoms
        self.dt_ns = dt_ns
        self.out_filename = out_filename
        self.out_dirname = out_dirname
        self.output_file = None
        self.prev_cm = None

        self._reference_atoms_for_alignment = self.reference.atoms.filter(self.atoms_selector_for_alignment)
        self._reference_atoms_for_alignment_coords = self._reference_atoms_for_alignment.coords.values

    def before_first_iteration(self, frame: Frame) -> None:
        """
        this function will be called before the first iteration over trajectory
        :param frame:
        :return:
        """
        self._frame_atoms_for_alignment = frame.atoms.filter(self.atoms_selector_for_alignment)

        # open file to write rmsd values
        os.makedirs(self.out_dirname, exist_ok=True)
        self.out_file = open(os.path.join(self.out_dirname, self.out_filename), "w")
        self.out_csvfile = csv.writer(self.out_file)
        self.out_csvfile.writerow(["time_ns", "rmsd"])

    def after_last_iteration(self, exc_type, exc_value, traceback) -> bool:
        """
        this function will be called after the last iteration over trajectory
        :param exc_type:
        :param exc_value:
        :param traceback:
        :return:
        """
        # close file
        self.out_file.close()
        return False

    def __call__(self, frame: Frame) -> Frame:
        """
        this function will be called for each iteration over trajectory
        :param frame:
        :return:
        """
        alignment = self._frame_atoms_for_alignment.alignment_to(self._reference_atoms_for_alignment)  # alignment
        crd = self._frame_atoms_for_alignment.coords.values.copy()
        crd = crd @ alignment.matrix3d().T + alignment.vector3d().values  # get coordinates

        frame.atoms.coords.apply(alignment)
        rmsd = calc_rmsd(self._reference_atoms_for_alignment_coords, crd)
        self.out_csvfile.writerow([frame.index * self.dt_ns, rmsd])

        return frame
