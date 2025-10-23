/*
 * Copyright (c) 2025 Andreas Burkhart.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

··········································································································
: __  __             _   _         ____                          _   _                  _                :
:|  \/  | __ _ _ __ | |_| | ___   / ___|___  _ ____   _____  ___| |_(_) ___  _ __      / \   _ __  _ __  :
:| |\/| |/ _` | '_ \| __| |/ _ \ | |   / _ \| '_ \ \ / / _ \/ __| __| |/ _ \| '_ \    / _ \ | '_ \| '_ \ :
:| |  | | (_| | | | | |_| |  __/ | |__| (_) | | | \ V /  __/ (__| |_| | (_) | | | |  / ___ \| |_) | |_) |:
:|_|  |_|\__,_|_| |_|\__|_|\___|  \____\___/|_| |_|\_/ \___|\___|\__|_|\___/|_| |_| /_/   \_\ .__/| .__/ :
:                                                                                           |_|   |_|    :
··········································································································

The large-scale mantle convection app as presented in the preprint (https://arxiv.org/abs/2506.04157) (currently in review).

Requires the Boost library, which is a header-only library as the app uses the TerraNeo plate provider headers.

In case the tools from src/terraneo at a later time switch to not being header only, you also might need to set HYTEG_TERRANEO_MODULE=ON.

··················································
: ____  _          _       _                     :
:|  _ \(_)___  ___| | __ _(_)_ __ ___   ___ _ __ :
:| | | | / __|/ __| |/ _` | | '_ ` _ \ / _ \ '__|:
:| |_| | \__ \ (__| | (_| | | | | | | |  __/ |   :
:|____/|_|___/\___|_|\__,_|_|_| |_| |_|\___|_|   :
··················································

This is a revised version of the code that had to be changed and updated such that it can run in the current master branch.

Due to the aforementioned past and possible future changes, the app might produce slightly different numerical results.

The original code pertaining to the preprint / publication above can be found at: https://www.doi.org/10.5281/zenodo.15497635

·····················································
: _   _                 _                           :
:| | | | _____      __ | |_ ___    _   _ ___  ___ _ :
:| |_| |/ _ \ \ /\ / / | __/ _ \  | | | / __|/ _ (_):
:|  _  | (_) \ V  V /  | || (_) | | |_| \__ \  __/_ :
:|_| |_|\___/ \_/\_/    \__\___/   \__,_|___/\___(_):
·····················································

The models in the Reproducibility folder are tied to the preprint / publication and often times use specific operators generated only for specific viscosities.

Files featuring a "MCMAIN_" tag use a general implementation usable for arbitrary choices of rho, eta and initial temperature. The multitude of models in 

the main folder are to be interpreted as blueprints on how to implement certain features within the context of the mantle convection app.

Of course all of these features could also have been implemented in a templated / if constexpr / ifdef / macro fashion but that would have lead to a completely

unreadable code file containing all features. The multitude of files provided here are also not a perfect solution but I wanted to avoid providing a codegenerator

as part of the app. I genuinely hope that this compromise works well enough.

In order to choose a model file answer the following questions:

1) Should the model use the frozen velocity approximation (compressibility term on the RHS) or solve a generalised asymmetric saddle point system (compressibility term on the RHS)?

2) Should the model use a advection-diffusion operator splitting in which the advection is solved via a modified method of characteristics (MMOC) or use SUPG stabilisation without operator splitting?

3) Should the model use the inverse viscosity scaled Mass, VCycle BFBT or WBFBT Schur complement approximation?

4) Should the model domain be 2D Annulus or a 3D spherical shell?

5) Should the model use a CG + GAMG PETSc coarse grid solver for the A block?

When in doubt: Use Asymmetric, MMOC, VCycle BFBT, non PETSc.

After you have chosen a model file you can select a parameter file (usually parameters_MC). You can either pass through the parameterfile name as a command line parameter or put

the name of your chosen parameter in "ScenarioSelector.prm". For 2D you should use the specified 2D parameter files.


Within the parameter file you can control all other parts of the model.

A hopefully somewhat complete documentation of how you can control the model via the parameter file:

- minLevel and maxLevel control the minimum and maximum refinement level of the Grid

- lowMemoryMode controls if the model uses the temporary function manager to save memory.
     If alwaysDestroyTemporaryFunctions is set to true then temporary functions are immediately deleted after use (slower but saves the highest amount of memory).

- vtk and bp4 control whether the model outputs vtk or bp4 files during the model run. These files are viewable in e.g. Paraview. Bp4 requires ADIOS2.

- loadCheckpointOnStart controls whether a the models tries to start from an existing checkpoint.
     loadCheckpointNumber allows you to control the checkpoint which the model will try to load (checkpoints are matched by filename and timestep).

- writeCheckpoints enables or disables the writing of checkpoints during the model run.
     defaultToUsingAdiosCheckpoints defines whether the model should try to save and load ADIOS2 checkpoints before trying the alternative FileVector (ADIOS2 not required) checkpointing.

- useMyrsInsteadOfTimeStepsAsOutputFrequency defines whether the model tries to save checkpoints and vtk/bp4 files every N timesteps or Myrs of model run time.
     writeFrequencyOutput and writeFrequencyCheckpoint define the number of timesteps or Myrs (integer value!) between checkpoint or vtk/bp4 outputs.

- boundaryTolerance is uses to set the boundary flags during storage generation. All DoF which are at most boundaryTolerance away from the analytical boundary are set as boundary DoFs.

- temperatureExtrapolationOrder and velocityExtrapolationOrder define the order of extrapolation used to predict the temperature and velocity at the next time step.
     For example order velocityExtrapolationOrder = 1 uses the current and last velocity together with the respective time step sizes to calculate a linear extrapolation.
     If the model does not have enough past states and timesteps to calculate the extrapolation up to the specified order then a extrapolation of the maximum possible is used.
     For example an extrapolation of second order requires the current and two past states. In the first time step the model would therefore use an extrapolation of order 0
     for the first time step. Past states and timesteps are also saved within checkpoints so this only affects the first few timesteps.

- BDFOrder sets the order of the BDF time discretisation used in the diffusion solve (MMOC) or advection-diffusion solve (SUPG). Currently 1 and 2 are supported.

- densityDerivativeOrder and viscosityDerivativeOrder are currently unused parameters that would allow you to approximate e.g. an approximation of the density derivative
  accurate up to the given order based on the density in past time steps. This could for example be used for the projected density approximation as described in Gassmöller 2020.
  The functionality to calculate these approximations is technically already finished. Leave at 0 for now.

- maxSteps defines a maximum number of time steps after which the model should end.

- minTimestepMyrs and maxTimestepMyrs allow you to constrict the timestep size to a given interval. Time steps falling outside of the interval are replaced by the closes value inside the interval.
  Set both to the same value for fixed time step size.

- CFL allows you to control the CFL constant used to determine the timestep size.
     useGlobalCFL allows you to control whether the maximal time step size is calculated for each element and afterwards the largest possible time step size is selected (useGlobalCFL false)
     or if  we should just use HyTeGs global mesh qualities like hMax to determine an approximation (useGlobalCFL true). When in doubt use false.

- randomInitialGuessU and randomInitialGuessP allow you to control whether the model should start with a random or zero initial guess in velocity and pressure.

- maximumSaddlePointSolverIterations allows you to set a limit how often the saddle point solver can be called per time step. Usually the solver only stops once it has reached
  its tolerance (e.g. for FGMRES) so you can leave this at 1. If your saddle point solver does not enforce a given tolerance (e.g. when it consists of a single VCycle) you can use this parameter to set
  a maximum number of solver calls.
     absoluteResidualToleranceOuterSaddlePointSolverLoop and relativeResidualToleranceOuterSaddlePointSolverLoop allow you to set the tolerances you aim for. Usually this is set to the same value as
     the tolerances set for the Saddle Point FGMRES Solver.

- blockpreconditionerType allows you to define the type of blockpreconditioner that you want to use. 
  0 = Inexact Uzawa, 1 = Adjoint Inexact Uzawa, 2 = Block Factorisation Approximation, 3 = Symmetric Uzawa

- shearHeatingCutoff allows you to set a cutoff for the shear heating in meters (if used, shear heating is set to zero within shearHeatingCutoff meters from the surface)

- agglomerateACoarseSolver defines whether the coarse grid solver of the A block should be solved on a lower amount of ranks to avoid excessive communication.
  This only makes sense for a very high amount of ranks. The number of ranks on which the coarse grid problem should be solved can be set in the "ABlock Agglomerated CG Coarse Grid Solver"
  options in the parameter file. This can only be used if your file does not use the PETSc coarse grid solver.

- disableOldA allows you to disable the usage of a different A block operator on lower levels (see preprint / publication). If disables the same A operator is used on all levels (specifically the NewAType one).

- OldAMaxLevel allows you to control the refinement level up to which the different A block operator is used (specifically the OldAType one). You can also use negative values.
  Negative value means OldAMaxLevel = min( minLevel + abs(OldAMaxLevel), maxLevel )

- chebyshevUpdateRateMyrs allows you to control how often (in terms of Myrs) the spectral radii of the Chebyshev smoothers are updated during the run.
  Value <= 0 disables the feature and always updates the spectral radii after every time step. Note: Since the chebyshev spectral radii are always
  updated after loading a checkpoint using this option might create slightly numerically different results as if the original had continued ( e.g. the
  orginal model would not have updated the spectral radii after a time step and after reloading they were updated anyway ).
  If useMyrsInsteadOfTimeStepsAsOutputFrequency is true and writeFrequencyOutput or writeFrequencyCheckpoint are a multiple of chebyshevUpdateRateMyrs,
  this problem should not occur!

- AsymmetricPreconditionerMode allows you to define whether the block preconditioner is evaluated for the symmetrical system (even if you solve a generalised asymmetric saddle point problem)
  or if the B block should be replaced by B+C during block preconditioner evaluation (see preprint / publication). Only applicable if you have chosen an AsymmetricSystem model file.
  0 = Asymmetric system is solved with the preconditioner for the symmetric system,
  unmodified (divergence) B Block gets used for the block preconditioner and BFBT Schur Complement approx.
  1 = Asymmetric system is solved with the preconditioner for the asymmetric system, 
  modified B Block gets used for the block preconditioner and BFBT Schur Complement approx.  
  Note: For mode 1 the BFBT suboperator should be solved via a GMRES solver. This is already supported (see Utility/Solver/Schur/SchurBFBTSolver.hpp and choose FGMRES instead of CG).

- WBFBTType allows you to set the type of WBFBT Schur complement approximation. 1 = Approx Schur Op with Poisson Neumann Prec, 2 = Poisson Neumann Only.  Only applicable if you have chosen a WBFBT model file.
  When in doubt: Leave this at type 2.
     WBFBT_ar and WBFBT_al allow you to set an asymmetric scaling for the WBFBT Schur complement approximation (see preprint / publication and Rudi2017). When in doubt: Set to 1.0.

- SUPG_scaling allows you to set a constant scaling for the SUPG terms. See the docstrings of the SUPG operators for hints on how SUPG is implemented specifically.
  Only applicable if you have chosen a SUPG model file. When in doubt: Set to 1.0.

- minVisc and maxVisc allow you to set relative bounds for the viscosity. For this to take effect you need to use the BoundedViscosity (see Utility/Viscosity/BoundedViscosity.hpp)
  class as part of your viscosity setup. minVisc and maxVisc are relative to the nondimensionalisation constant for the viscosity.

- usePlates allows you to control whether the TerraNeo plate oracle should be used to set the time dependent surface Dirichlet boundary conditions for the velocity between time steps.
  Only applicable if surfaceBoundaryType = 1.
     fileTopologies and fileReconstructions allow you to set a plate topology and reconstruction file that should be loaded.
     plateVelocityScaling allows you to set a constant inverse scaling of the plate velocities at the surface. For a discussion of plate scaling see Colli2020.
     loopPlateAge allows you to continue a simulation after the plate reconstructions given by the loaded files are over. In this case the plates just start 
     at the beginning again (implemented via a modulo with the maximum plate age).

- surfaceBoundaryType and CMBBoundaryType allow you to set the type of boundary condition you want to use at the surface and CMB. 1 = Dirichlet (Noslip or PlateData), 2 = Neumann, 3 = FreeSlip.
  Currently only 1 and 3 make sense.

- nTan2D, nRad2D, nTan3D and nRad3D allow you to control how the coarse mesh is generated in 2D and 3D, specifically how the elements are generated in tangential and radial direction.
  For 3D (ntan - 1) must be a power of 2. See "src/hyteg/mesh/MeshInfo.hpp" or the generated hyteg documentation at "https://hyteg.pages.i10git.cs.fau.de/hyteg/classhyteg_1_1MeshInfo.html"
  for more explanations of the nTan and nRad parameters.

- temperatureSurface, temperatureCMB, radiusSurface, radiusCMB, etaRef, rhoRef, C_pRef, alphaRef, gRef, GammaRef, uRef and kRef allow you to set
  the (nondimensionalisation) reference constants for the model as detailed in the appendix of the preprint / publication.

- const_H, const_alpha, const_K_T, const_k and const_C_p allow you to define the nondimensional internal heating, thermal expansivity, isothermal bulk modulus, specific heat capacity
  as constants used in the model

- rockChemicalCompositionParameter, depthDependency and additiveOffSet allow you to configure the Frank–Kamenetskii type viscosity (if used), see "Utility/Viscosity/ExponentialViscosity.hpp"
  for more details.

- ArrheniusC1 and ArrheniusC2 allow you to configure the Arrhenius type viscosity (if used), see "Utility/Viscosity/ArrheniusViscosity.hpp" for more details.
  Note that the ArrheniusViscosity class also uses the rockChemicalCompositionParameter, depthDependency and additiveOffSet parameters.

- initialTemperatureFilter, initialTemperatureSteepness, degreeMinSH, degreeMaxSH and buoyancyFactor allow you to configure the spherical harmonics (initial) temperature model (if used).
  See "Utility/Temperature/SphericalHarmonicsTemperature.hpp" for more details.

- relativeTemperatureNoiseFactor allow you to control the maximum relative temperature change of the RelativeRandomTemperature class (if used). See "Utility/Temperature/RelativeRandomTemperature.hpp"
  for more details.

- rhoSurf allows you to set the surface density for the ExponentialDensity class (if used). See "Utility/Density/ExponentialDensity.hpp" for more details.

All remaining parameters pertain to the solvers that can potentially be used as part of your model and should hopefully be self explanatory. Note that you can define a "prefix" for any solver
found in "Utility/Solver". This allows you to define multiple solvers of the same type that load a different set of parameters. See, e.g., the "Inexact ABlock CG Outer Loop Solver" parameters
will be loaded if you create an instance of the "ABlockCGOuterLoopSolver" solver class with prefix "Inexact".


··················································································································································································································································
: _   _                 _               _                                  _                _ _                   _                   _ _            ___     _       _ _   _       _   _                                      _                  :
:| | | | _____      __ | |_ ___     ___| |__   __ _ _ __   __ _  ___    __| | ___ _ __  ___(_) |_ _   _    __   _(_)___  ___ ___  ___(_) |_ _   _   ( _ )   (_)_ __ (_) |_(_) __ _| | | |_ ___ _ __ ___  _ __   ___ _ __ __ _| |_ _   _ _ __ ___ :
:| |_| |/ _ \ \ /\ / / | __/ _ \   / __| '_ \ / _` | '_ \ / _` |/ _ \  / _` |/ _ \ '_ \/ __| | __| | | |   \ \ / / / __|/ __/ _ \/ __| | __| | | |  / _ \/\ | | '_ \| | __| |/ _` | | | __/ _ \ '_ ` _ \| '_ \ / _ \ '__/ _` | __| | | | '__/ _ \:
:|  _  | (_) \ V  V /  | || (_) | | (__| | | | (_| | | | | (_| |  __/ | (_| |  __/ | | \__ \ | |_| |_| |_   \ V /| \__ \ (_| (_) \__ \ | |_| |_| | | (_>  < | | | | | | |_| | (_| | | | ||  __/ | | | | | |_) |  __/ | | (_| | |_| |_| | | |  __/:
:|_| |_|\___/ \_/\_/    \__\___/   \___|_| |_|\__,_|_| |_|\__, |\___|  \__,_|\___|_| |_|___/_|\__|\__, ( )   \_/ |_|___/\___\___/|___/_|\__|\__, |  \___/\/ |_|_| |_|_|\__|_|\__,_|_|  \__\___|_| |_| |_| .__/ \___|_|  \__,_|\__|\__,_|_|  \___|:
:                                                         |___/                                   |___/|/                                   |___/                                                       |_|                                      :
··················································································································································································································································

As an example, lets investigate how one would implement a new viscosity given a fictitious (space and temperature dependent) viscosity map.

Lets assume you want to implement the viscosity

    eta(x,T) := eta0(x) * f(x,T).

A viscosity model can be created as a combination of multiple "TemperatureDependentViscosityModel" child classes. First lets assume that
eta0(x) only depends on the radius and that you have a given depth dependent viscosity profile saved in a "myCustomViscosity.csv" file.
This csv file can be loaded using the "Utility/Data/DataLoader.hpp" class:

    auto dataVisc = MantleConvection::loadCSV( "Utility/Data/ViscosityProfiles/myCustomViscosity.csv" );

The data in the csv file can consists of two sorted vectors, one representing the radius in meters and the other the viscosity in Pa s.
You also need to specify the data type in the first column, see e.g., "Utility/Data/ViscosityProfiles/viscosityStotz2018Cutoff_Arrhenius.csv"
as an example how to do this.

Now we can nondimensionalise the viscosity data

    MantleConvection::nondimensionaliseCSV( parameters_, dataVisc );

where "parameters_" is the "walberla::config::Config::BlockHandle" of your loaded parameterfile ( you can get this by calling Model->getParameters() ).

Now we can create our eta0 viscosity model based on the data:

    auto ViscProfile = std::make_shared< MantleConvection::LinearViscosityInterpolation >( ND_, parameters_, dataVisc.range.front(), dataVisc.values.front() );

To implement eta(x,T), you can now create a child class of TemperatureDependentViscosityModel ( "Utility/Viscosity/TemperatureDependentViscosityModel.hpp" ) called "CustomViscosity"
and override the "evaluate" method taking both the point in space and the nondimensional temperature ( not necessarily between 0 and 1, see e.g. ExponentialViscosity.hpp ).
This viscosity model class should now take another "TemperatureDependentViscosityModel" object called "ModifiedModel" during construction that it can use during
evaluation. In our case "ModifiedModel" would refer to eta0(x). See, e.g., "ExponentialViscosity.hpp" for an example. Of course our child class also needs to calculate
f(x,T) as part of the overriden "evaluate" method.

    real_t evaluate( const hyteg::Point3D& x, real_t temp ) override
    {
      // get the temperature to be between 0 and 1
      temp -= temperatureSurface_;

      return ModifiedModel_->evaluate( x, temp ) * f( x, temp );
    }

Now we can create the eta(x,T) that we want:

    auto CustomVisc = std::make_shared< CustomViscosity >( ND_, parameters_, ViscProfile );

Lets assume that we now want to also bound the viscosity to an interval. We can do this by using the "BoundedViscosity" class:

    auto BoundedVisc = std::make_shared< BoundedViscosity >( ND_, parameters_, CustomVisc );

The actual viscosity bounds (relative to the viscosity reference constant) can be configured in the parameter file.

Now our viscosity is ready to be used and we can set it as part of our model initialisation:

    auto& viscosityModel = BoundedVisc;

The initial temperature, reference tempature, density, viscosity and potentially pressure model classes all follow this design pattern and can be
easily combined with eachother.

Note: The "DataLoader" class can also read in two dimensional tables like thermodynamic pressure and temperature dependent lookup tables.

··········································································································
: _____                          _                     _   _                                   _ _       :
:|  __ \                        | |                   | | | |                                 | | |      :
:| |__) |___ _ __  _ __ ___   __| |_   _  ___ ___     | |_| |__   ___      _ __ ___  ___ _   _| | |_ ___ :
:|  _  // _ \ '_ \| '__/ _ \ / _` | | | |/ __/ _ \    | __| '_ \ / _ \    | '__/ _ \/ __| | | | | __/ __|:
:| | \ \  __/ |_) | | | (_) | (_| | |_| | (_|  __/    | |_| | | |  __/    | | |  __/\__ \ |_| | | |_\__ \:
:|_|  \_\___| .__/|_|  \___/ \__,_|\__,_|\___\___|     \__|_| |_|\___|    |_|  \___||___/\__,_|_|\__|___/:
:           | |                                                                                          :
:           |_|                                                                                          :
··········································································································

For an explanation how to reproduce the results in the preprint / publication download https://www.doi.org/10.5281/zenodo.15497635 and follow
the instructions included in the README.txt file. The instructions should be valid for both the orginal and the master branch version found here.