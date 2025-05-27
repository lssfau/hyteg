#### Analytical convergence studies

The results for the three cases can be obtained by inputting the right parameters in the parameter file `3DSphShellConvStudies.prm` and adjusting `cRotPenalty` value for the freeslip-freeslip cases,

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
