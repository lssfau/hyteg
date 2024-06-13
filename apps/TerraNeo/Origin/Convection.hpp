/*
 * Copyright (c) 2024 Eugenio D'Ascoli, Andreas Burkhart,
 * Nils Kohl, Hamish Brown, Ponsuganth Ilangovan, Marcus Mohr.
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

#pragma once

#include <core/Environment.h>
#include <core/config/Config.h>
#include <core/mpi/MPIManager.h>

#ifdef HYTEG_BUILD_WITH_ADIOS2
#include "hyteg/checkpointrestore/ADIOS2/AdiosCheckpointExporter.hpp"
#include "hyteg/checkpointrestore/ADIOS2/AdiosCheckpointImporter.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosWriter.hpp"
#endif

#include "hyteg/MeshQuality.hpp"
#include "hyteg/dataexport/SQL.hpp"
#include "hyteg/dataexport/TimingOutput.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/gridtransferoperators/P2toP1Conversion.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticVectorProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticVectorRestriction.hpp"
#include "hyteg/memory/MemoryAllocation.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/numerictools/CFDHelpers.hpp"
#include "hyteg/numerictools/L2Space.hpp"
#include "hyteg/numerictools/SpectrumEstimation.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/p2functionspace/P2TransferOperators.hpp"
#include "hyteg/petsc/PETScBlockPreconditionedStokesSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/ChebyshevSmoother.hpp"
#include "hyteg/solvers/FGMRESSolver.hpp"
#include "hyteg/solvers/FullMultigridSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/SubstitutePreconditioner.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/controlflow/SolverLoop.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesBlockPreconditioners.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesVelocityBlockBlockDiagonalPreconditioner.hpp"
#include "hyteg/solvers/solvertemplates/StokesFSGMGSolverTemplate.hpp"
// HOG generated HyTeG operator
#include "hyteg_operators/operators/div_k_grad/P2ElementwiseDivKGradIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/k_mass/P1ElementwiseKMass.hpp"
#include "hyteg_operators/operators/k_mass/P1ElementwiseKMassIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/k_mass/P2ToP1ElementwiseKMassIcosahedralShellMap.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesFullOperator.hpp"
#include "hyteg_operators_composites/viscousblock/P2ViscousBlockLaplaceOperator.hpp"
// Particle transport
#include "convection_particles/vtk/ParticleVtkOutput.h"
#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"
#include "mixed_operator/ScalarToVectorOperator.hpp"
#include "mixed_operator/VectorMassOperator.hpp"
// TerraNeo helpers for initialisation and data structures
#include "terraneo/helpers/RadialProfiles.hpp"
#include "terraneo/helpers/TerraNeoParameters.hpp"
#include "terraneo/initialisation/TemperatureInitialisation.hpp"
#include "terraneo/operators/P2P1StokesOperatorWithProjection.hpp"
#include "terraneo/operators/P2TransportTALAOperator.hpp"

namespace terraneo {

inline TerraNeoParameters TN;

class ConvectionSimulation
{
 public:
   // no default constructor allowed
   ConvectionSimulation() = delete;

   ConvectionSimulation( const walberla::Config::BlockHandle& mainConf );
   ~ConvectionSimulation(){};

   void init();

   void step();

   void                        outputTimingTree();
   real_t                      viscosityFunction( const Point3D& x, real_t Temperature );
   real_t                      densityFunction( const Point3D& x );
   real_t                      diffPreFactorFunction( const Point3D& x );
   real_t                      calculateStokesResidual( uint_t level );
   real_t                      referenceTemperatureFunction( const Point3D& x );
   const SimulationParameters& getSimulationParams();

 private:
   typedef P2P1TaylorHoodFunction< real_t >                                                  StokesFunction;
   typedef P2Function< real_t >                                                              ScalarFunction;
   typedef hyteg::operatorgeneration::P2P1StokesFullIcosahedralShellMapOperator              StokesOperator;
   typedef P2P1StokesFullIcosahedralShellMapOperatorFS                                       StokesOperatorFS;
   typedef hyteg::operatorgeneration::P2ViscousBlockLaplaceIcosahedralShellMapOperator       BlockLaplaceOperator;
   typedef hyteg::operatorgeneration::P1ElementwiseKMassIcosahedralShellMap                  SchurOperator;
   typedef hyteg::operatorgeneration::P2ElementwiseDivKGradIcosahedralShellMap               DiffusionOperator;
   typedef VectorMassOperator< real_t, P1VectorFunction, P1ElementwiseBlendingMassOperator > P1MassOperatorVelocity;
   typedef hyteg::operatorgeneration::P2ToP1ElementwiseKMassIcosahedralShellMap              FrozenVelocityOperator;

   void setupDomain();
   void setupBoundaryConditions();
   void setupFunctions();
   void initialiseFunctions();
   void setupSolversAndOperators();

   void setupOutput();
   void setupStokesRHS();
   void setupEnergyRHS();
   void updateRefViscosity();
   void updatePlateVelocities( StokesFunction& );
   void solveStokes();
   void solveEnergy();

   void dataOutput();
   void normalFunc( const Point3D& p, Point3D& n );

   // Declare boundary condition objects for unknown functions for temperature and velocity

   BoundaryCondition bcVelocity;
   BoundaryCondition bcTemperature;

   // Corresponding Boundary UIDs for velocity and temperature to setup

   BoundaryUID idSurface;
   BoundaryUID idCMB;

   BoundaryUID idSurfaceT;
   BoundaryUID idCMBT;

   // function and solver initialization

   // Scalar fields for temperature
   std::shared_ptr< ScalarFunction > temperature;
   std::shared_ptr< ScalarFunction > temperaturePrev;
   std::shared_ptr< ScalarFunction > temperatureDev;
   std::shared_ptr< ScalarFunction > temperatureTmp;
   std::shared_ptr< ScalarFunction > temperatureReference;

   // scalar field for viscosity
   std::shared_ptr< ScalarFunction > viscosityFE;
   std::shared_ptr< ScalarFunction > viscosityFEInv;

   // Scalar field for RHS of energy equation
   std::shared_ptr< ScalarFunction > energyRHS;
   std::shared_ptr< ScalarFunction > energyRHSWeak;

   //field containing ones
   std::shared_ptr< ScalarFunction > onesFE;

   // scalar field for density
   std::shared_ptr< ScalarFunction > densityFE;

   // scalar field for Diffusion pre factor
   std::shared_ptr< ScalarFunction > diffusionFE;
   std::shared_ptr< ScalarFunction > adiabaticTermCoeff;
   std::shared_ptr< ScalarFunction > shearHeatingTermCoeff;
   std::shared_ptr< ScalarFunction > constEnergyCoeff;
   std::shared_ptr< ScalarFunction > surfTempCoeff;

   // Vector and scalar field storing velocity and pressure fields (stokes equation)
   std::shared_ptr< StokesFunction >             stokesLHS;
   std::shared_ptr< P1VectorFunction< real_t > > stokesLHSP1;
   std::shared_ptr< P1VectorFunction< real_t > > stokesLHSP1Weak;
   // Vector and scalar field storing velocity and pressure fields (of the Prev timestep) (stokes equation)
   std::shared_ptr< StokesFunction > stokesLHSPrev;
   // Vector and scalar field storing right hand side of stokes equation
   std::shared_ptr< StokesFunction > stokesRHS;
   // Temporary vector and scalar filed storing RHS of stokes equation
   std::shared_ptr< StokesFunction > stokesTmp;
   // Vector field storing the inward normal of the sphere
   std::shared_ptr< P2VectorFunction< real_t > > inwardNormal;
   //Vector field for Grad(Rho)/Rho in RHS mass
   std::shared_ptr< P2VectorFunction< real_t > > gradRhoOverRho;
   // Scalar field storing the magnitude of the velocity
   std::shared_ptr< ScalarFunction > velocityMagnitudeSquared;

   std::shared_ptr< P2VectorFunction< real_t > > oppositeGravityField;
   // Temporary scalar fields
   std::shared_ptr< ScalarFunction > scalarTmp;

   // Storage for primitives (includes functionality for distributed computing)
   std::shared_ptr< PrimitiveStorage > storage;

   // Solvers

   std::shared_ptr< CGSolver< DiffusionOperator > > diffusionSolver;
   // std::shared_ptr< CGSolver< P2DiffusionOperatorWrapper > > diffusionSolverTest;
   std::shared_ptr< FGMRESSolver< StokesOperator > >                     stokesSolver;
   std::shared_ptr< Solver< StokesOperatorFS > >                         stokesSolverFS;
   std::shared_ptr< CGSolver< P2TransportIcosahedralShellMapOperator > > transportSolverTALA;

   // Operators
   std::shared_ptr< StokesOperator >                         stokesOperator;
   std::shared_ptr< StokesOperatorFS >                       stokesOperatorFS;
   std::shared_ptr< SchurOperator >                          schurOperator;
   std::shared_ptr< MMOCTransport< ScalarFunction > >        transportOperator;
   std::shared_ptr< P2TransportIcosahedralShellMapOperator > transportOperatorTALA;
   std::shared_ptr< DiffusionOperator >                      diffusionOperator;
   std::shared_ptr< P2ElementwiseBlendingMassOperator >      P2MassOperator;
   std::shared_ptr< P1MassOperatorVelocity >                 MassOperatorVelocityP1;
   std::shared_ptr< P2ProjectNormalOperator >                projectionOperator;

   std::shared_ptr< FrozenVelocityOperator > frozenVelocityRHSX;
   std::shared_ptr< FrozenVelocityOperator > frozenVelocityRHSY;
   std::shared_ptr< FrozenVelocityOperator > frozenVelocityRHSZ;

   std::shared_ptr< InexactUzawaPreconditioner< StokesOperator, StokesOperator::ViscousOperator_T, SchurOperator > >
                                               uzawaSmoother;
   std::shared_ptr< Solver< StokesOperator > > coarseGridSolver;
   std::shared_ptr< hyteg::VTKOutput >         vtkOutput;

   // ADIOS2 data output
#ifdef HYTEG_BUILD_WITH_ADIOS2
   std::shared_ptr< AdiosWriter > _output;

   std::shared_ptr< AdiosCheckpointExporter > checkpointExporter;
   std::shared_ptr< AdiosCheckpointImporter > checkpointImporter;
#endif

   // std::functions for various functionalities
   std::function< real_t( const Point3D& ) >         densityFunc;
   std::function< real_t( const Point3D& ) >         diffFactorFunc;
   std::function< real_t( const Point3D&, real_t ) > viscosityFunc;
   std::function< real_t( const Point3D& ) >         referenceTemperatureFct;

   std::function< real_t( const Point3D& ) > normalX = []( const Point3D& x ) { return -x[0] / x.norm(); };
   std::function< real_t( const Point3D& ) > normalY = []( const Point3D& x ) { return -x[1] / x.norm(); };
   std::function< real_t( const Point3D& ) > normalZ = []( const Point3D& x ) { return -x[2] / x.norm(); };

   std::shared_ptr< RadialProfile >                             viscosityProfiles;
   std::shared_ptr< RadialProfile >                             temperatureProfiles;
   std::shared_ptr< TemperatureInitializationParameters >       temperatureInitParams;
   std::shared_ptr< std::function< real_t( const Point3D& ) > > temperatureReferenceFct;
};

} // namespace terraneo
