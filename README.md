# cwipiFoam
Establishes CWIPI coupling of OpenFOAM-v2212 and Nektar++ Acoustic Solver version 5.3.0

Installation:  
1) Compile Nektar++ with AcousticSolver and CWIPI; install into $HOME/opt  
2) Compile OpenFOAM-v2212
3) Create a file named loadNektar.sh
4) Add the following lines to loadNektar.sh:

export PATH=$HOME/opt/bin:$PATH
export LIBRARY_PATH=$HOME/opt/lib64:$LIBRARY_PATH
export LIBRARY_PATH=$HOME/opt/lib64/nektar++:$LIBRARY_PATH
export LD_LIBRARY_PATH=$HOME/opt/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$HOME/opt/lib64/nektar++:$LD_LIBRARY_PATH

5) Run source loadNektar.sh
7) Run wmake cwipiPstream
8) Run wmake rhoCentralFoam_cwipi
