#  Tauola-bbb 

[![GitHub release (latest by date)](https://img.shields.io/github/v/release/jimmjohn/TAUOLA-bbb-25-05-2020?include_prereleases)](https://github.com/jimmjohn/TAUOLA-bbb-25-05-2020)
[![GitHub](https://img.shields.io/github/license/jothepro/doxygen-awesome-css)![GitHub Repo stars](https://img.shields.io/github/stars/jimmjohn/tauola-bbb)



The Tauola-bbb package includes the Fortran version of the Tauola package. The documentation provides a detailed explanation of the purpose of each subroutine.



## Installation

Tauola is a standalone package, but for more sophisticated exampes may require the following packages.

To run default examples simply execute:

"make"

"./taumain.exe"

in any of the tauola-bbb/demo-* directories.



To use tauola-bbb content in other projects follow instructions (README files) from dedicated tauola-bbb/patch-* folders

NOTE: if instructions of patch-babar-validation are followed, changes of files of the installation are introduced.

​      Be sure to backup tauola-bbb folder as well as photos/photos.f file for easy restoration of the original files. 

The compilation configurations are set in the adequate Makefiles.



For `Tauola` with c++ interface installation, you need the `HepMC2`, `HepMC3`, `PYTHIA8`, `LHAPDF` and `MC-TESTER`.

For `HepMC2` and `HepMC3` package, you download the source code and create a build directory. Then run the following command from the build directory.

```sh
cmake -DCMAKE_INSTALL_PREFIX=<path>/HepMC2/HepMC-2.06.11/ <path>/HepMC2/HepMC-2.06.11/ \
-Dmomentum:STRING=GEV -Dlength:STRING=MM
make 
make install
```

For `HepMC3`, we have to just change the folder and the version of `HepMC3`

For `MC-TESTER` installation, go to the the unzipped location and run following command

```sh
./configure --prefix=<path>/MC-TESTER/mc-tester/ \
--with-HepMC=<path>/HepMC2/HepMC-2.06.11/ \
--with-HepMC3=<path>/HepMC3/HepMC3-3.3.0
make
make install
```

For the installation of `Tauola`, you can do

```sh
./configure --prefix=<path>/Tauola --with-hepmc=<path>/HepMC2/HepMC-2.06.11/ \
--with-hepmc3=<path>/HepMC3/HepMC3-3.3.0 --with-lhapdf=<path>/LHAPDF/LHAPDF-6.5.5 \
--with-pythia8=<path>/Pythia8/pythia8313 --with-mc-tester=<path>/MC-TESTER/mc-tester
make
```



## Directory tree

The fortran code of the library is strored in the main directory, some hadronic currents are stored in sub-directories:

   tauola-bbb/new-currents

the C++ code is stored in directory:

   tauola-bbb/tauola-c

If it would not be used truncated version is stored in:

   tauola-bbb/tauola-no-c

See README files of these directories for more details.

\---

Other folders contain supplementary libraries allowing standalone run of prepared examples.



## Examples

The following examples (sub-directories) are prepared with corresponding README files:

- tauola-bbb/demo-babar
- tauola-bbb/demo-lfv
- tauola-bbb/demo-redefine
- tauola-bbb/demo-pairs

You can got to the particular directory and execute the make command and run the program. The makefile links the `glibk` library for plotting, `photos` library for the radiative corrections and the `jetset` library for printing the event records in pythia format.



## Patches for other projects

There are some patches for installation the code into other projects. Again, each sub-directory features its own README 

- tauola-bbb/patch-KK-face
- tauola-bbb/patch-tauolapp
- tauola-bbb/patch-babar-validation




## Written by

Zbigniew Was, 

**[zbigniew.was@ifj.edu.pl](mailto:zbigniew.was@ifj.edu.pl)**


|                        Read Next |
| -------------------------------: |
| [taumain](docs/theme/taumain.md) |

