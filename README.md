<div align="justify">

# nucleosome-trajectory-processing

This repository contains the scripts and example data for analysis of MD trajectory of nucleosome core particle (NCP).
You may find our protocol for MD simulation of NCP at https://github.com/OOLebedenko/nucleosome-md-simulation

### System requirements

Key packages and programs:

- Linux platform
- [python3](https://www.python.org/) (version >= 3.8)
- [gcc compiler](https://gcc.gnu.org/) (7.3.0 <= version <= 9.5.0)
- [g++ compiler](https://gcc.gnu.org/) (7.3.0 <= version <= 9.5.0)
- [cmake](https://www.gnu.org/software/make/manual/make.html) (version == 4.3)
- [dssp](https://swift.cmbi.umcn.nl/gv/dssp/) (version == 3.0.0-3build1)

### Installation dependencies

The key package for analysis of MD trajectory is in-house python library [pyxmolpp2](https://sizmailov.github.io/pyxmolpp2/api/python/install.html)

```code-block:: bash
    # install requried packages via
    sudo apt-get install cmake g++-7 gcc-7 python3 python3-dev libnetcdf-dev 
    
    # create virtual enviroment
    python3 -m venv ./venv
    source ./venv/bin/activate
    pip install --upgrade pip
    python -m pip install -U setuptools wheel pip
    
    # export path to C++ compilers needed for pyxmolpp2
    export CC=/usr/bin/gcc-7; export CXX=/usr/bin/g++-7 
    
    # install packages
    pip install -r requirements.txt
```

### Run analysis of MD trajectory

To start processing of the MD trajectory, please, see the github page for the corresponding type of analysis:

1) [15N relaxation rates](15N_relaxation_rates/README.md)


</div>



