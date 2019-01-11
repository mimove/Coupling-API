# Instructions to couple OpenFOAM with CWIPI

cwipiIcoFoam is the a solver based on icoFoam, with a few additional lines with calls to CWIPI functions:
```
cwipiCoupling(mesh);
...
cwipiSend(mesh, U);
...
cwipiWait();
```
These functions are defined in the library cwipiPstream.

The cavity folder is a tutorial case for coupling cwipiIcoFoam with the AcousticSolver of Nektar++.



## Compiling the OpenFOAM library and solver

To compile the cwipiPstream library, the option file cwipiPstream/Make/options needs to be modified with the correct cwipi installation path.
```
wmake cwipiPstream
wmake cwipiIcoFoam
```


## Running the tutorial case

Within the cavity turorial folder, CFD is the subfolder for the OpenFOAM side, while CAA is the subfolder for the Nektar++ side.
In the following example, the case is run on 5 cores, 2 of which are assigned to OpenFOAM and 3 to Nektar++.

```
cd cavity
decomposePar -case CFD
mpirun --output-filename log -wd CFD -np 2 cwipiIcoFoam -parallel : -wd CAA -np 3 AcousticSolver --verbose --cwipi FOAM_APE test_icofoam.xml Conditions_APE_CWIPI.xml
```
