### 2025 HyTeG/TerraNeo Benchmarks

This folder contains the code used to produce benchmark results. The codes are split into three folders corresponding to the analytical, gephysical benchmarks and scalability studies,

* `static`
    - Contains code for order of convergence studies wrt the analytical solution

* `geobench`
    - Contains code for compressible and nonlinear benchmark on the unit square and the symmetrical plume benchmarks on the SPH shell

* `scalability`
    - Contains a simple app which runs the FGMRES solver for the Stoke system and computes the runtime
