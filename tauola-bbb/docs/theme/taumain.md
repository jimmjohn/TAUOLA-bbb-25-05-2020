#  taumain

[![GitHub release (latest by date)](https://img.shields.io/github/v/release/jimmjohn/tauola-bbb?include_prereleases)](https://github.com/jothepro/doxygen-awesome-css/releases/latest)
[![GitHub](https://img.shields.io/github/license/jothepro/doxygen-awesome-css)![GitHub Repo stars](https://img.shields.io/github/stars/jimmjohn/tauola-bbb)







## Installation

For `Tauola` installation, you need the `HepMC2`, `HepMC3`, `PYTHIA8`, `LHAPDF` and `MC-TESTER`.

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

The examples are given in the folder `demo-babar` ,  `demo-lfv` ,  `demo-redifine`  folders. You can got to the particular directory and execute the make command and run the program. The makefile links the `glibk` library for plotting, `photos` library for the radiative corrections and the `jetset` library for printing the event records in pythia format.




## Written by

Zbigniew Was, 

**[zbigniew.was@ifj.edu.pl](mailto:zbigniew.was@ifj.edu.pl)**


|                        Read Next |
|---------------------------------:|
| [Extensions](docs/extensions.md) |

</div>
