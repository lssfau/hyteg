#### Geophysical Benchmarks

This section contains code for the geophysical benchmarks which are categorised into three parts,

* Compressible simulation with TALA formulation on an unit square identical to the TALA case from King 2010 is done in `CompressibleTALA.cpp`. Some parameters such as the Rayleigh number and the Dissipation number can be controlled from the parameter file `CompressibleTALA.prm`

* A single convection cell with nonlinear rheology on an unit square identical to Tosi 2015 case 4 is done in `TosiNonlinear.cpp` and the parameters can be adjusted in the corresponding parameter file

* Symmetrical plume benchmarks on the Spherical shell with the incompressible variable viscous Stokes system identical to Euen 2023 is done in `SphericalShellBench.cpp`
