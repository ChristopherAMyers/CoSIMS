# CoSIMS
CoSIMS is a Collision Simulator for Ion Mobility Spectrometry. This is a trajectory based method for calculating molecular colision cross-sections for use with IMS-MS experiments. We welcome the submission of bugs, improvements, and suggestions to the code.

If you find this code usefull, we ask that you please cite the following paper:

>Christopher A. Myers, Rebecca J. D’Esposito, Daniele Fabris, Srivathsan V. Ranganathan, and Alan A. Chen. CoSIMS: An Optimized Trajectory Based Collision Simulator for Ion Mobility Spectrometry. *J. Phys. Chem. B.*, May 2019, 10.1021/acs.jpcb.9b01018.
https://pubs.acs.org/doi/10.1021/acs.jpcb.9b01018


## Installation

### Obtaining the Code
To obtain the program from the GitHub repository:
```bash
git clone https://github.com/ChristopherAMyers/CoSIMS.git
```
Alternatively, download the `.zip` file from GitHub and extract the contents manually.

### Building
- **Requirements:**  
  - `g++` version 4.8.1 or higher  
  - OpenMP library for multi-threading  
  - Linux system (Windows not tested)

To build the program:
```bash
make
```
Or compile using multiple cores:
```bash
make -j
```
This creates a statically linked binary at:
```
/path/to/cosims/install/bin/cosims
```

Manual compilation:
```bash
g++ src/*.cpp -o cosims -std=c++11 -fopenmp -O3
```

---

## Running CoSIMS

Basic run command:
```bash
./cosims -m moleculeFile.xxx
```
Where `xxx` = `pdb` or `mfj`.

Additional options:
```bash
./cosims -m moleculeFile.xxx -[option] [value]
```

---

### Molecule File Types

#### PDB Files
- Uses `ATOM` and `HETATM` entries.  
- Assumes all entries belong to one structure.  
- Columns 78–79 specify element types; if blank, atom name columns 13–16 are used.  

Example:
```
ATOM      1  OH5’ C5       1       3.429  -7.861   3.641  1.00  0.00  H
```

#### MFJ Files
- Follows [MOBCAL](https://www.sciencedirect.com/science/article/abs/pii/0009261496006388) syntax with added integer charge support.  

---

### Random Seed
Set with `-s` option:
```bash
./cosims -m file.pdb -s 12345
```
Uses 64-bit Mersenne Twister RNG.  
Same seed may yield different results with multiple threads due to dynamic allocation.

---

### Number of Threads
Set with `-n` option:
```bash
./cosims -m file.pdb -n 4
```
- Should be ≤ number of physical CPU cores.  
- Identical runs may produce different results due to thread-safe RNG.

---

### Atom Types
Set with `-a` option to provide custom parameters.  
Format:
```
#name  symbol  mass   sigma    epsilon    intMass
Carbon C       12.010 2.043    0.03090114 12
Hydrogen H      1.008 2.043    0.01498936 1
```

---

### Input File
Use `-i` option for runtime configuration without recompiling.

Example input file:
```
# CoSIMS input file
name myName
dt 0.01
temp 298
traj 50
iter 15
threads 1
seed 12345
dispersion cutoff true
dispersion radius 40
multipole true
multipole order 1
```

---

## Available Program Options

| Option               | Description                                                   | Type      | Default                     |
|----------------------|---------------------------------------------------------------|-----------|------------------------------|
| `name`               | Project name / output file prefix                             | string    | `log charge`                 |
| `dt`                 | Verlet integration time step (ps)                             | float     | 0.01                         |
| `temp`               | Temperature (K)                                               | float     | 298                          |
| `traj`               | Trajectories per CCS integral (×10³)                           | int       | 50                           |
| `iter`               | Number of CCS integrals to average                             | int       | 10                           |
| `print rate`          | Print CCS progress every N trajectories                        | int       | 20000                        |
| `print mol`           | Print rotated molecule to XYZ                                 | bool      | false                        |
| `dir`                | Output directory                                              | string    | Current directory            |
| `charge`             | Uniform total charge (e)                                       | int       | 0                            |
| `threads`            | OpenMP threads (overrides `-n`)                                | int       | 1                            |
| `seed`               | RNG seed (overrides `-s`)                                      | int64     | random device()               |
| `dispersion cutoff`   | LJ cutoff enabled                                             | bool      | true                         |
| `dispersion radius`   | LJ cutoff radius (Å)                                          | float     | Based on vdW radii            |
| `multipole`           | Enable multipole approximation                                | bool      | false                        |
| `multipole radius`    | Radius for multipole approximation (Å)                        | float     | 20.0                         |
| `multipole order`     | 0 = monopole, 1 = dipole, 2 = quadrupole                      | int       | 1                            |
| `max cluster size`    | Max radial distance for atom clusters                          | float     | 4.0                          |

---

### References
1. M. F. Mesleh *et al.*, *The Journal of Physical Chemistry*, 100(40):16082–16086, 1996.  
2. A. A. Shvartsburg, M. F. Jarrold, *Chemical Physics Letters*, 261(1):86–91, 1996.  
