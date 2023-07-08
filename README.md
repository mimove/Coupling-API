# cwipiFoam
Establishes CWIPI coupling of OpenFOAM-v2212 and Nektar++ Acoustic Solver version 5.3.0

# Currently the entropy material derivative in the APE-4 continuity equation is not working

Installation:  
1) Compile Nektar++ with AcousticSolver and CWIPI; install into $HOME/opt  
2) Compile OpenFOAM-v2212
3) Run the following commands (or add to the bottom of your .bashrc file):

export PATH=$HOME/opt/bin:$PATH

export LIBRARY_PATH=$HOME/opt/lib64:$LIBRARY_PATH

export LIBRARY_PATH=$HOME/opt/lib64/nektar++:$LIBRARY_PATH

export LD_LIBRARY_PATH=$HOME/opt/lib64:$LD_LIBRARY_PATH

export LD_LIBRARY_PATH=$HOME/opt/lib64/nektar++:$LD_LIBRARY_PATH

6) Run wmake in cwipiPstream
7) Run wmake in rhoCentralFoam_cwipi
