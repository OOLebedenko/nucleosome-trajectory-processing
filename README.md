<div align="justify">

# nucleosome-trajectory-processing

This repository contains the scripts and example data for analysis of MD trajectory of nucleosome core particle (NCP).
You can find our protocol for MD simulation of NCP at https://github.com/OOLebedenko/nucleosome-md-simulation

### System requirements

Key packages and programs:

- Linux platform (tested on ubuntu 20.04 and ubuntu 22.04)
- [python3](https://www.python.org/)
- [gcc compiler](https://gcc.gnu.org/) (7.0 <= version <= 11.0)
- [g++ compiler](https://gcc.gnu.org/) (7.0 <= version <= 11.0)
- [cmake](https://www.gnu.org/software/make/manual/make.html)
- [dssp](https://swift.cmbi.umcn.nl/gv/dssp/) (version == 3.0.0-3build1)

### Installation dependencies

The key package for analysis of MD trajectory is in-house python
library [pyxmolpp2](https://sizmailov.github.io/pyxmolpp2/api/python/install.html)

```code-block:: bash
# install dssp==3.0.0-3build1
sudo add-apt-repository "deb [arch=amd64] http://archive.ubuntu.com/ubuntu focal main universe"
sudo apt-get install dssp=3.0.0-3build1

# install system dependencies
sudo apt-get install g++ gcc cmake python3 python3-dev libnetcdf-dev 

# create virtual enviroment
python3 -m venv ./venv
source ./venv/bin/activate
pip install --upgrade pip
pip install -U setuptools wheel pip

# install python packages
pip install -r requirements.txt
```

### Run MD analysis

To start processing of the MD trajectory, please, see the github page for the relevant type of calculations:

1) [15N relaxation rates](15N_relaxation_rates/README.md)

</div>



