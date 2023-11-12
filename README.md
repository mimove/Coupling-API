# cwipiFoam
Establishes CWIPI coupling of OpenFOAM-v2212 and Nektar++ Acoustic Solver version 5.3.0

The aim of this library is to enable memory-based source coupling of compressible OpenFOAM solvers and Nektar++'s AcousticSolver.  The OpenFOAM solver is used to calculate the unsteady hydrodynamic source fields, which are then used as input sources by AcousticSolver to compute the propagation of noise into the far-field.  The method of memory-based coupling avoids the I/O bottleneck associated with file-based data transfer since the data is exchanged through system memory via CWIPI.

An edited version of the rhoCentralFoam solver has been provided as an example of how to edit the source code of compressible OpenFOAM solvers to use this CWIPI coupling.

To use cwipiRhoCentralFoam, one must first run the solver in decoupled mode by setting the entry cwipiSwitch to false in system/controlDict:

cwipiSwitch       false;

The solver should be run in this mode until the start-up transient period has passed.  The user should then add time-averaging for the base flow fields, namely for U, c, T, rho, s and L to the bottom of the controlDict file:

{

    functions
    {
      fieldAverage0
      {
        type            fieldAverage;
        libs            (fieldFunctionObjects);
        fields
        (
          U
          {
            mean          true;
            prime2Mean    false;
            base          time;
            allowRestart  true;
          }
          c
          {
            mean          true;
            prime2Mean    false;
            base          time;
            allowRestart  true;
          }
          T
          {
            mean          true;
            prime2Mean    false;
            base          time;
            allowRestart  true;
          }
          rho
          {
            mean          true;
            prime2Mean    false;
            base          time;
            allowRestart  true;
          }
          s
          {
            mean          true;
            prime2Mean    false;
            base          time;
            allowRestart  true;
          }
          L
          {
            mean          true;
            prime2Mean    false;
            base          time;
            allowRestart  true;
          }
        );
        writeControl      writeTime;
      }
    }
    
}

The solver should then be run, once again in decoupled mode, until the time-averaged flow fields are statistically stationary.  Once this is the case, the solver can then be run in coupled mode with the following switches:

cwipiSwitch       true;
cwipiDim          2;
cwipiStep         1;
cwipiLambVector   true;
cwipiEntropy      true;
cwipiDsDt         true;

The above combination of parameters enables CWIPI coupling for a Nektar solution grid with 2 spatial dimensions, sending all of the acoustic source fields at each time step.  The individual acoustic sources can be turned on/off with their respective named entries (cwipiLambVector, cwipiEntropy and cwipiDsDt).

The cwipiStep parameter governs the frequency with which the source fields are sent to AcousticSolver; this parameter must be defined such that cwipiStep(openfoam) * timeStep(openfoam) = ReceiveSteps(Nektar++) * TimeStep(Nektar++), i.e. that the two solvers step through the same amount of simulated time for every invokation of the send() method.  As an example, a value of timeStep for the OpenFOAM solver of 1e-6 and a TimeStep value pf 1e-7 in Nektar++ necessitates cwipiStep = 1 and ReceiveStep = 10, since the ratio of time steps is 1/10.

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
