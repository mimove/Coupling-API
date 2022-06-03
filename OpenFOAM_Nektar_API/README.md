# Instructions to couple OpenFOAM and Nektar++ via CWIPI

cwipiIcoFoam is a solver based on icoFoam, with a few additional lines with calls to CWIPI functions:
```
cwipiCoupling(mesh);
...
cwipiSend(mesh, LPrime);
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


## Cavity tutorial case

Within the cavity turorial folder, CFD is the subfolder for the OpenFOAM side, while CAA is the subfolder for the Nektar++ side.
In the following example, the case is run on 5 cores, 2 of which are assigned to OpenFOAM and 3 to Nektar++. This case was run on versions 5.0 of Nektar++ and 1712 of OpenFOAM.

This case is not intended to generate any meaningful aeroacoustic result, but to show how the sources are sent from OpenFOAM and Nektar++.


To run the case, the following commands should be run in a terminal:

```
cd cavity
decomposePar -case CFD
mpirun --output-filename log -wd CFD -np 2 cwipiIcoFoam -parallel : -wd CAA -np 3 AcousticSolver --verbose --cwipi FOAM_APE test_icofoam.xml Conditions_APE_CWIPI.xml
```

Once the simulation is finished, two check points would have been saved. Those check points should be postprocessed using the appropriate OpenFOAM and Nektar++ commands to obtain the .plt or .vtk files. These files should then be opened with a visualisation software (such as ParaView or Tecplot) to verify that the fields have been correctly interpolated. The three components of the OpenFOAM field LPrime should be the same as the fields F_0_u, F_0_v, F_0_w of Nektar++. If this is the case, the one-way coupling between the two solvers is working correctly.

-----------------------
If you use publish any work based on this repository, we kindly ask you to cite our publication: 

"Moratilla-Vega, M.A., Angelino, M., Xia, H. and Page, G.J., 2022. An open-source coupled method for aeroacoustics modelling. Computer Physics Communications, p.108420. DOI: https://doi.org/10.1016/j.cpc.2022.108420"

-----------------------

Any question/doubt about the coupling can be sent to mimove14@gmail.com
