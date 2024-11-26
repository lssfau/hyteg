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

#include "hyteg/HytegDefinitions.hpp"

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
#include "hyteg/gridtransferoperators/P2toP2QuadraticInjection.hpp"
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
#include "hyteg/solvers/solvertemplates/StokesFSGMGUzawaSolverTemplate.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"
// HOG generated HyTeG operator
#include "hyteg_operators/operators/div_k_grad/P2ElementwiseDivKGradIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/k_mass/P1ElementwiseKMass.hpp"
#include "hyteg_operators/operators/k_mass/P1ElementwiseKMassIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/k_mass/P2ToP1ElementwiseKMassIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/terraneo/P2VectorToP1ElementwiseFrozenVelocityIcosahedralShellMap.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesEpsilonOperator.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesFullOperator.hpp"
#include "hyteg_operators_composites/viscousblock/P2ViscousBlockLaplaceOperator.hpp"
// Particle transport
#include "convection_particles/vtk/ParticleVtkOutput.h"
#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"
#include "mixed_operator/ScalarToVectorOperator.hpp"
#include "mixed_operator/VectorMassOperator.hpp"
// TerraNeo helpers for initialisation and data structures
#include "terraneo/helpers/InterpolateProfile.hpp"
#include "terraneo/helpers/RadialProfiles.hpp"
#include "terraneo/helpers/TerraNeoParameters.hpp"
#include "terraneo/helpers/Viscosity.hpp"
#include "terraneo/initialisation/TemperatureInitialisation.hpp"
#include "terraneo/operators/P2P1StokesOperatorWithProjection.hpp"
#include "terraneo/operators/P2TransportRHSOperator.hpp"
// #include "terraneo/operators/P2TransportTALAOperator.hpp"
#include "terraneo/operators/P2TransportTALAOperatorStd.hpp"

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
   real_t                      densityFunction( const Point3D& x );
   real_t                      diffPreFactorFunction( const Point3D& x );
   real_t                      calculateStokesResidual( uint_t level );
   real_t                      calculateEnergyResidual( uint_t level );
   real_t                      referenceTemperatureFunction( const Point3D& x );
   const SimulationParameters& getSimulationParams();

   void updateNonDimParameters( const Point3D& x );

 private:
   typedef P2P1TaylorHoodFunction< real_t >                                            StokesFunction;
   typedef P2Function< real_t >                                                        ScalarFunction;
   typedef P2P1TaylorHoodFunction< real_t >                                            StokesFunctionP2P1;
   typedef P2Function< real_t >                                                        ScalarFunctionP2;
   typedef P2VectorFunction< real_t >                                                  VectorFunctionP2;
   typedef P1Function< real_t >                                                        ScalarFunctionP1;
   typedef hyteg::operatorgeneration::P2P1StokesFullIcosahedralShellMapOperator        StokesOperator;
   typedef P2P1StokesP1ViscosityFullIcosahedralShellMapOperatorFS                      StokesOperatorFS;
   typedef hyteg::operatorgeneration::P2ViscousBlockLaplaceIcosahedralShellMapOperator BlockLaplaceOperator;
   typedef hyteg::operatorgeneration::P1ElementwiseKMassIcosahedralShellMap            SchurOperator;
   typedef hyteg::operatorgeneration::P2ElementwiseDivKGradIcosahedralShellMap         DiffusionOperator;
   typedef hyteg::operatorgeneration::P2ToP1ElementwiseKMassIcosahedralShellMap        FrozenVelocityOperator;

   typedef hyteg::operatorgeneration::P2VectorToP1ElementwiseFrozenVelocityIcosahedralShellMap FrozenVelocityFullOperator;

   // typedef P2P1StokesP1ViscosityFullIcosahedralShellMapOperatorFS StokesOperatorP1Visc;

   void setupDomain();
   void setupBoundaryConditions();
   void setupFunctions();
   void initialiseFunctions();
   void setupSolversAndOperators();

   void setupOutput();
   void setupStokesRHS();
   void setupEnergyRHS();
   void updateViscosity();
   void updatePlateVelocities( StokesFunction& );
   void solveStokes();
   void solveEnergy();

   void dataOutput();
   void normalFunc( const Point3D& p, Point3D& n );
   void oppositeGravity( const Point3D& p, Point3D& n );

   real_t adiabaticCoefficientFunction( const Point3D& x );
   real_t constantEnergyCoefficientFunction( const Point3D& x );
   real_t surfaceTempCoefficientFunction( const Point3D& x );

   // Declare boundary condition objects for unknown functions for temperature and velocity

   BoundaryCondition bcVelocity;
   BoundaryCondition bcTemperature;

   // Corresponding Boundary UIDs for velocity and temperature to setup

   BoundaryUID idSurface;
   BoundaryUID idCMB;

   BoundaryUID idSurfaceT;
   BoundaryUID idCMBT;

   // function and solver initialization

   enum class BoundaryConditionType
   {
      VELOCITY_BOUNDARY_CONDITION,
      TEMPERATURE_BOUNDARY_CONDITION,
      NO_BOUNDARY_CONDITION
   };

   /** Initialise functions
    * 
    * To create a new function, just add into the correct dict and it can be used from the container anywhere in the app
    * The tuple denotes < 
    * std::string           = (Name of the function),
    * uint_t                = (Minimum level),
    * uint_t                = (Maximum level),
    * BoundaryConditionType = (Velocity or Temperature or No Boundary condition)
    * >
    */
   std::vector< std::tuple< std::string, uint_t, uint_t, BoundaryConditionType > > p2p1StokesFunctionDict = {
       { "VelocityFE", 0u, 0u, BoundaryConditionType::VELOCITY_BOUNDARY_CONDITION },
       { "VelocityFEPrev", 0u, 0u, BoundaryConditionType::VELOCITY_BOUNDARY_CONDITION },
       { "StokesRHS", 0u, 0u, BoundaryConditionType::VELOCITY_BOUNDARY_CONDITION },
       { "StokesTmp1", 0u, 0u, BoundaryConditionType::VELOCITY_BOUNDARY_CONDITION },
       { "StokesTmp2", 0u, 0u, BoundaryConditionType::VELOCITY_BOUNDARY_CONDITION },
       { "StokesTmpProlongation", 0u, 0u, BoundaryConditionType::VELOCITY_BOUNDARY_CONDITION } };
   std::map< std::string, std::shared_ptr< StokesFunctionP2P1 > > p2p1StokesFunctionContainer;

   std::vector< std::tuple< std::string, uint_t, uint_t, BoundaryConditionType > > p2ScalarFunctionDict = {
       { "TemperatureFE", 0u, 0u, BoundaryConditionType::TEMPERATURE_BOUNDARY_CONDITION },
       { "TemperaturePrev", 0u, 0u, BoundaryConditionType::TEMPERATURE_BOUNDARY_CONDITION },
       { "TemperatureDev", 0u, 0u, BoundaryConditionType::TEMPERATURE_BOUNDARY_CONDITION },
       { "Temperature[K]", 0u, 0u, BoundaryConditionType::TEMPERATURE_BOUNDARY_CONDITION },
       { "ReferenceTemperature[K]", 0u, 0u, BoundaryConditionType::TEMPERATURE_BOUNDARY_CONDITION },
       { "ViscosityFE", 0u, 0u, BoundaryConditionType::NO_BOUNDARY_CONDITION },
       { "ViscosityFEInv", 0u, 0u, BoundaryConditionType::NO_BOUNDARY_CONDITION },
       { "Viscosity[Pas]", 0u, 0u, BoundaryConditionType::NO_BOUNDARY_CONDITION },
       { "EnergyRHSWeak", 0u, 0u, BoundaryConditionType::TEMPERATURE_BOUNDARY_CONDITION },
       { "DensityFE", 0u, 0u, BoundaryConditionType::NO_BOUNDARY_CONDITION },
       { "ShearHeatingTermCoeff", 0u, 0u, BoundaryConditionType::NO_BOUNDARY_CONDITION },
       { "VelocityMagnitudeSquared", 0u, 0u, BoundaryConditionType::NO_BOUNDARY_CONDITION } };
   std::map< std::string, std::shared_ptr< ScalarFunctionP2 > > p2ScalarFunctionContainer;

   std::vector< std::tuple< std::string, uint_t, uint_t, BoundaryConditionType > > p2VectorFunctionDict = {};
   std::map< std::string, std::shared_ptr< VectorFunctionP2 > >                    p2VectorFunctionContainer;

   // Storage for primitives (includes functionality for distributed computing)
   std::shared_ptr< PrimitiveStorage > storage;

   // Solvers

   std::shared_ptr< CGSolver< DiffusionOperator > >                      diffusionSolver;
   std::shared_ptr< FGMRESSolver< StokesOperator > >                     stokesSolver;
   std::shared_ptr< Solver< StokesOperatorFS > >                         stokesSolverFS;
   std::shared_ptr< CGSolver< P2TransportIcosahedralShellMapOperator > > transportSolverTALA;

   std::shared_ptr< Solver< StokesOperatorFS::ViscousOperatorFS_T > > stokesABlockSmoother;

   // Operators
   std::shared_ptr< StokesOperator >                         stokesOperator;
   std::shared_ptr< StokesOperatorFS >                       stokesOperatorFS;
   std::shared_ptr< SchurOperator >                          schurOperator;
   std::shared_ptr< MMOCTransport< ScalarFunction > >        transportOperator;
   std::shared_ptr< P2TransportIcosahedralShellMapOperator > transportOperatorTALA;
   // std::shared_ptr< P2TransportRHSIcosahedralShellMapOperator > transportOperatorRHS;
   std::shared_ptr< DiffusionOperator >                 diffusionOperator;
   std::shared_ptr< P2ElementwiseBlendingMassOperator > P2MassOperator;
   std::shared_ptr< P2ProjectNormalOperator >           projectionOperator;
   std::shared_ptr< P2toP2QuadraticProlongation >       p2ProlongationOperator;

   std::shared_ptr< P2toP2QuadraticInjection > p2InjectionOperator;

   std::shared_ptr< P2TransportP1CoefficientsIcosahedralShellMapOperator > transportOperatorP1Coefficients;

   std::shared_ptr< FrozenVelocityOperator > frozenVelocityRHSX;
   std::shared_ptr< FrozenVelocityOperator > frozenVelocityRHSY;
   std::shared_ptr< FrozenVelocityOperator > frozenVelocityRHSZ;

   std::shared_ptr< FrozenVelocityFullOperator > frozenVelocityRHS;

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
   std::function< real_t( const Point3D& ) > densityFunc;
   std::function< real_t( const Point3D& ) > diffFactorFunc;
   std::function< real_t( const Point3D& ) > adiabaticCoeffFunc;
   std::function< real_t( const Point3D& ) > constEnergyCoeffFunc;
   std::function< real_t( const Point3D& ) > surfTempCoeffFunc;
   std::function< real_t( const Point3D& ) > referenceTemperatureFct;
   std::function< real_t( const Point3D& ) > oppositeGravityX;
   std::function< real_t( const Point3D& ) > oppositeGravityY;
   std::function< real_t( const Point3D& ) > oppositeGravityZ;

   std::function< real_t( const Point3D& ) > normalX = []( const Point3D& x ) { return -x[0] / x.norm(); };
   std::function< real_t( const Point3D& ) > normalY = []( const Point3D& x ) { return -x[1] / x.norm(); };
   std::function< real_t( const Point3D& ) > normalZ = []( const Point3D& x ) { return -x[2] / x.norm(); };

   std::shared_ptr< RadialProfile >                             viscosityProfiles;
   std::shared_ptr< RadialProfile >                             temperatureProfiles;
   std::shared_ptr< TemperatureInitializationParameters >       temperatureInitParams;
   std::shared_ptr< std::function< real_t( const Point3D& ) > > temperatureReferenceFct;
};

} // namespace terraneo
