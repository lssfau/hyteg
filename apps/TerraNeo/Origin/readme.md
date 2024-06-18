# TerraNeo

(View this file via [HyTeG's documentation](https://hyteg.pages.i10git.cs.fau.de/hyteg/index.html) (🚧 add correct link to this page) for proper math rendering.)

A full mantle circulation modeling application.

Core features:

* mantle convection modeled via the truncated anelastic liquid approximation (TALA)
* surface velocity boundary conditions through plate reconstruction data
* free-slip boundary conditions at the core-mantle boundary (CMB)
* matrix-free implementation for maximum scalability

---

## Getting started

### Compiling

Make sure you have enabled the TerraNeo module via CMake by setting
```
HYTEG_TERRANEO_MODULE=ON
```

It's optional (but recommended) to install ADIOS2 for checkpointing and efficient parallel output. Enable it via
```
HYTEG_BUILD_WITH_ADIOS2=ON
```

Then build the `TerraNeo` target, e.g., via:
```
make TerraNeo -j8
```
(`-j` triggers a parallel build).

### Running the app

The application parameters are set through a parameter file. Running the app via

```
mpirun -np 8 ./TerraNeo
```

will default to using the default (`parameters.prm`) parameter file. You can create your own file and pass it via the 
command line instead:

```
mpirun -np 8 ./TerraNeo my_parameter_file.prm
```

---

## Documentation

### Model and formulation

> 🚧 references and REVIEW!

The TALA models the mantle via a coupled system of the compressible full Stokes equation and the compressible 
energy equation:

\f[
\begin{align}
    - \nabla \cdot \tau + \nabla p' &= - \alpha \bar{\rho} T' g, \\
    \nabla \cdot u &= - \frac{1}{\bar{\rho}} \nabla \bar{\rho} \cdot u^*, \\
    \bar{\rho} C_p \left( \frac{\partial T}{\partial t} + u \cdot \nabla T \right) - \nabla \cdot (k \nabla T) &= \bar{\rho} H + \tau : \xi + \alpha \bar{\rho} T(u \cdot g)
\end{align}
\f]

with 

\f[
\begin{align}
    \epsilon &:= \frac{1}{2} \left( \nabla u + (\nabla u)^T \right), \\
    \xi &:= \epsilon - \frac{1}{3} (\nabla \cdot u) I, \\
    \tau &:= 2 \eta \xi
\end{align}
\f]

and \f$ I \f$ being the identity.

> 🚧 properly document non-dimensionalization

### Symbol and parameter legend

Barred quantities (e.g., \f$ \bar{\rho} \f$) refer to a reference or background state while primed quantities (e.g.,
\f$ \rho' \f$) are deviations from the reference (i.e., \f$ \rho = \bar{\rho} + \rho' \f$).

| ***symbol***             | ***meaning***                           |
|--------------------------|-----------------------------------------|
| \f$ u \f$                | velocity                                |
| \f$ p \f$                | pressure                                |
| \f$ \rho \f$             | density                                 |
| \f$ T \f$                | temperature                             |
| \f$ \eta \f$             | dynamic viscosity                       |
| \f$ g \f$                | gravitational acceleration              |
| \f$ \mathrm{Ra} \f$      | Rayleigh number                         |
| \f$ \mathrm{Di} \f$      | dissipation number                      |
| \f$ C_p \f$              | heat capacity                           |
| \f$ k \f$                | diffusivity                             |
| \f$ \alpha \f$           | thermal expansion                       |
| \f$ T_\mathrm{S} \f$     | surface temperature                     |
| \f$ T_{\mathrm{CMB}} \f$ | temperature at the CMB                  |
| \f$ \Delta T \f$         | \f$ T_{\mathrm{CMB}} - T_\mathrm{S} \f$ |

> 🚧 this is certainly not a complete list

### Boundary conditions

Both, surface and CMB boundary conditions can be selected via the parameter file. The general approach the present 
mantle circulation model is to prescribe the surface velocity boundary conditions through plate reconstruction data and
to prescribe free-slip boundary conditions at the core-mantle boundary. However, for testing purposes, other boundary 
conditions can be chosen.

> 🚧 more details on plates

### Initial conditions

Several options are available for the initial temperature field. They include random initialization and smooth spherical 
harmonics. 

### Viscosity

The depth- and temperature-dependent viscosity can be selected in two steps. The (constant) background viscosity
\f$ \eta_0 \f$ can be set using one of the provided radial viscosity profiles under `data/terraneo/viscosityProfiles/`.

The viscosity is then computed via incorporation of an additional temperature-dependent term. For details revisit the 
documentation in `data/terraneo/viscosityProfiles/`.

### Discretization

#### Domain and mesh

The thick spherical shell is discretized via an [icosahedral meshing approach](https://hyteg.pages.i10git.cs.fau.de/hyteg/classhyteg_1_1MeshInfo.html)
that generates a tetrahedral coarse mesh. The structure of this coarse mesh can be steered through parameters that define
the number of radial shells and the number of vertices in tangential direction.

The initial coarse mesh is the uniformly refined into a block-structured tetrahedral mesh (see [Kohl 2023](https://doi.org/10.1080/17445760.2023.2266875)).

The elements are projected onto the phycial domain using the tranfinite element approach 
(see e.g., [Gordon 1973](https://doi.org/10.1007/BF01436298) and [Bauer 2017](https://doi.org/10.1016/j.apnum.2017.07.006)).

#### Finite element spaces

Velocity and pressure are discretized using the standard \f$ \mathbb{P}_2-\mathbb{P}_1 \f$ Taylor-Hood element pairing 
on the underlying tetrahedral grid. The temperature is discretized with \f$ \mathbb{P}_2 \f$ and the viscosity is 
embedded into a \f$ \mathbb{P}_2 \f$ space for efficiency.

### Solvers

#### Stokes

The discrete Stokes system is solved using either a block-diagonally multigrid precinditioned FGMRES solver or a 
monolithic multigrid method with inexact Uzawa relaxation (see [Kohl 2022a](https://doi.org/10.1137/20M1376005)).
Optimal strategies for different viscosity formulations are still under investigation. 

All multigrid methods are geometric multigrid methods that exploit the block-structure of the underlying mesh.
All compute kernels are matrix-free and have been generated and optimized using the 
[HyTeG Operator Generator](https://i10git.cs.fau.de/hyteg/hog) (see [Böhm 2024](https://arxiv.org/abs/2404.08371)).

#### Energy equation

The energy equation is treated with an Eulerian-Lagrangian method that treats the advective term with a modified method
of characteristics and the diffusive term with standard preconditioned Krylov solvers (see [Kohl 2022b](https://doi.org/10.1137/21M1402510)).

### I/O

> 🚧 ADIOS2/VTK
> 🚧 checkpointing

