from typing import Sequence, List, Union, Tuple
from pyxmolpp2 import Frame
from pyxmolpp2 import PdbFile, Trajectory, TrjtoolDatFile, TorsionAngleFactory
from pyxmolpp2.pipe import TrajectoryProcessor, ProcessedTrajectory
import csv
import os


class ConsumerCsv:

    def __init__(self, file):
        self.csv_file = csv.writer(file)

    def writerow(self, value_array):
        self.csv_file.writerow(value_array)

    def header(self, header):
        self.csv_file.writerow(header)


class ConsumerStdout:

    def __init__(self, csv_file):
        self.csv_file = csv_file

    def consume(self, vector):
        print(vector)


class VectorsExtractor(TrajectoryProcessor):

    def copy(self):
        return VectorsExtractor(self.get_selection, self.vector_consumers)

    def __init__(self, get_selection, vector_consumers):
        self.get_selection = get_selection
        self.vector_consumers = vector_consumers

    def __ror__(self, trajectory: Sequence[Frame]):
        return ProcessedTrajectory(trajectory, self)

    def __call__(self, frame: Frame):
        atoms_selection_1, atoms_selection_2 = self.get_selection(frame)
        vectors = atoms_selection_1.coords.values - atoms_selection_2.coords.values
        for vector_consumer, vector in zip(self.vector_consumers, vectors):
            vector_consumer.writerow(vector)
        return frame


class AngleExtractor(TrajectoryProcessor):

    def copy(self):
        return AngleExtractor(self.angle_name, self.consumers)

    def __init__(self, angle_name, angle_consumers):
        self.angle_name = angle_name
        self.consumers = angle_consumers

    def __ror__(self, trajectory: Sequence[Frame]):
        return ProcessedTrajectory(trajectory, self)

    def __call__(self, frame: Frame):
        angles = [TorsionAngleFactory.get(residue=residue, angle_name=self.angle_name) for residue in frame.residues]
        for consumer, angle in zip(self.consumers, angles):
            if angle is None:
                continue

            consumer.writerow([str(angle.value().to_standard_range().degrees)])


class OpenCsvAsVectorsExtractors:
    def __init__(self, ref, get_selection, filename_format="{atom2.residue.id.serial:02d}_{atom2.name}.csv",
                 out_dir="./"):
        self.files = []
        self.vector_consumer_files = []
        self.filename_format = filename_format
        self.get_selection = get_selection
        self.filenames = [
            filename_format.format(atom2=atom2)
            for atom2 in self.get_selection(ref)[1]
        ]
        self.out_dir = out_dir

    def __enter__(self):
        for fname in self.filenames:
            fout = open(os.path.join(self.out_dir, fname), "w")

            vector_consumer = ConsumerCsv(fout)
            vector_consumer.header(["x", "y", "z"])

            self.vector_consumer_files.append(vector_consumer)
            self.files.append(fout)
        return VectorsExtractor(vector_consumers=self.vector_consumer_files, get_selection=self.get_selection)

    def __exit__(self, exc_type, exc_value, exc_traceback):
        for fout in self.files:
            fout.close()


class OpenCsvAsAngleExtractors:

    def __init__(self, ref, angle_name, filename_format="{residue.id.serial:02d}.csv",
                 out_dir="./"):
        self.files = []
        self.angle_consumer_files = []
        self.filename_format = filename_format
        self.angle_name = angle_name
        self.filenames = [
            filename_format.format(residue=residue)
            for residue in ref.residues
        ]
        self.out_dir = out_dir

    def __enter__(self):
        for fname in self.filenames:
            fout = open(os.path.join(self.out_dir, fname), "w")

            vector_consumer = ConsumerCsv(fout)
            vector_consumer.header(["degree"])

            self.angle_consumer_files.append(vector_consumer)
            self.files.append(fout)
        return AngleExtractor(angle_consumers=self.angle_consumer_files, angle_name=self.angle_name)

    def __exit__(self, exc_type, exc_value, exc_traceback):
        for fout, fname in zip(self.files, self.filenames):
            if os.path.getsize(os.path.join(self.out_dir, fname)) == 0:
                fout.close()
                os.remove(os.path.join(self.out_dir, fname))
            fout.close()


class Run:
    def __ror__(self, trajectory: Sequence[Frame]):
        for _ in trajectory:
            pass
