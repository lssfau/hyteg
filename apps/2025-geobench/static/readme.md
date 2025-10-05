### Analytical convergence studies

The results for the three cases can be obtained by inputting the right parameters in the parameter file `3DSphShellConvStudies.prm` and adjusting `cRotPenalty` value for the freeslip-freeslip cases. In addition for the fs-fs cases, one needs to use
`sphNTan 2;` and `sphNRad 2;`, while for every other cases we use `sphNTan 5;`, `sphNRad 3;`.

* Delta forcing - Mixed BC
    - `mixed true;`
    - `delta true;`
    - `freeslip false;`

* Smooth forcing - Mixed BC
    - `mixed true;`
    - `delta false;`
    - `freeslip false;`

* Smooth forcing - Freeslip BCs
    - `mixed false;`
    - `delta false;`
    - `freeslip true;`

##### Kernel generation

The `rDash` value can be altered to place the delta function and the topography kernel value for a given `rDash` is also printed in the output logs. All kernels were generated with `sphNTan 5;`, `sphNRad 3;` and `maxLevel 4;`