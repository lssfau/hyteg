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

#include "hyteg/Git.hpp"
#include "hyteg/MeshQuality.hpp"
#include "hyteg/dataexport/SQL.hpp"
#include "hyteg/dataexport/TimingOutput.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/gridtransferoperators/P0toP0AveragedInjection.hpp"
#include "hyteg/gridtransferoperators/P1toP0Conversion.hpp"
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

#include "sqlite/SQLite.h"
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
#include "terraneo/operators/P2P1StokesOperatorRotationOpgen.hpp"
#include "terraneo/operators/P2P1StokesOperatorWithProjection.hpp"
#include "terraneo/operators/P2TransportRHSOperator.hpp"
#include "terraneo/operators/TransportOperatorStd.hpp"
#include "terraneo/solvers/MCSolverBase.hpp"
#include "terraneo/solvers/StokesMCFGMRESSolver.hpp"
#include "terraneo/solvers/StokesMCFGMRESSolverWithProjection.hpp"
#include "terraneo/solvers/StokesMCUzawaSolver.hpp"
#include "terraneo/types/types.hpp"
#include "terraneo/utils/NusseltNumberOperator.hpp"

namespace terraneo {

inline TerraNeoParameters TN;

template < typename TemperatureFunction_T, typename ViscosityFunction_T >
class ConvectionSimulation
{
 public:
   using P2P1StokesFunction_T = P2P1TaylorHoodFunction< real_t >;

   using P2VectorFunction_T = P2VectorFunction< real_t >;
   using P1VectorFunction_T = P1VectorFunction< real_t >;

   using P2ScalarFunction_T = P2Function< real_t >;
   using P1ScalarFunction_T = P1Function< real_t >;
   using P0ScalarFunction_T = P0Function< real_t >;

   // no default constructor allowed
   ConvectionSimulation() = delete;

   ConvectionSimulation( const walberla::Config::BlockHandle& mainConf );
   ~ConvectionSimulation() {};

   void init();

   void step();

   void                        outputTimingTree();
   real_t                      densityFunction( const Point3D& x );
   real_t                      diffPreFactorFunction( const Point3D& x );
   real_t                      calculateStokesResidual( uint_t level );
   real_t                      calculateEnergyResidual( uint_t level );
   real_t                      calculatePressureResidual( uint_t level );
   real_t                      referenceTemperatureFunction( const Point3D& x );
   const SimulationParameters& getSimulationParams();

 private:
   typedef hyteg::operatorgeneration::P2P1StokesFullIcosahedralShellMapOperator        StokesOperator;
   typedef P2P1StokesP1ViscosityFullIcosahedralShellMapOperatorFS                      StokesOperatorFS;
   typedef hyteg::operatorgeneration::P2ViscousBlockLaplaceIcosahedralShellMapOperator BlockLaplaceOperator;
   typedef hyteg::operatorgeneration::P1ElementwiseKMassIcosahedralShellMap            SchurOperator;
   typedef hyteg::operatorgeneration::P2ElementwiseDivKGradIcosahedralShellMap         DiffusionOperator;

   typedef hyteg::operatorgeneration::P2VectorToP1ElementwiseFrozenVelocityIcosahedralShellMap FrozenVelocityFullOperator;

   using TransportOperator_T =
       std::conditional_t< std::is_same_v< TemperatureFunction_T, hyteg::P2Function< real_t > >,
                           terraneo::P2TransportIcosahedralShellMapOperator,
                           std::conditional_t< std::is_same_v< TemperatureFunction_T, hyteg::P1Function< real_t > >,
                                               terraneo::P1TransportIcosahedralShellMapOperator,
                                               std::nullptr_t > >;

   using StokesOperator_T =
       std::conditional_t< std::is_same_v< ViscosityFunction_T, hyteg::P0Function< real_t > >,
                           P2P1StokesP0ViscosityFullIcosahedralShellMapOperatorFS,
                           std::conditional_t< std::is_same_v< ViscosityFunction_T, hyteg::P1Function< real_t > >,
                                               P2P1StokesP1ViscosityFullIcosahedralShellMapOperatorFS,
                                               P2P1StokesFullIcosahedralShellMapOperatorFS > >;

   using StokesOperatorRotationOpgen_T =
       std::conditional_t< std::is_same_v< ViscosityFunction_T, hyteg::P0Function< real_t > >,
                           P2P1StokesP0ViscousOpgenRotationWrapper,
                           std::conditional_t< std::is_same_v< ViscosityFunction_T, hyteg::P1Function< real_t > >,
                                               P2P1StokesOpgenRotationWrapper,
                                               std::nullptr_t > >;

   using ProjectionOperator_T = P2ProjectNormalOperator;
   using RotationOperator_T   = P2RotationOperator;

   using P1ScalarMassOperator_T = operatorgeneration::P1ElementwiseMassIcosahedralShellMap;
   using P2ScalarMassOperator_T = operatorgeneration::P2ElementwiseMassIcosahedralShellMap;

   using P1VectorMassOperator_T = P1ElementwiseBlendingVectorMassOperator;
   using P2VectorMassOperator_T = P2ElementwiseBlendingVectorMassOperator;

   void setupDomain();
   void setupBoundaryConditions();
   void setupFunctions();
   void initialiseFunctions();
   void setupSolversAndOperators();

   void setupOutput();
   void setupStokesRHS();
   void setupEnergyRHS();
   void updateViscosity();
   void updatePlateVelocities( P2P1StokesFunction_T& );
   void solveStokes();
   void solveEnergy();
   void calculateHeatflow( const std::shared_ptr< RadialProfile >& );
   void calculateHeatflowIntegral( const std::shared_ptr< RadialProfile >& );
   void dataOutput();
   void outputCheckpoint();
   void initTimingDB();
   void normalFunc( const Point3D& p, Point3D& n );
   void oppositeGravity( const Point3D& p, Point3D& n );

   real_t adiabaticCoefficientFunction( const Point3D& x );
   real_t constantEnergyCoefficientFunction( const Point3D& x );
   real_t surfaceTempCoefficientFunction( const Point3D& x );

   // Declare boundary condition objects for unknown functions for temperature and velocity

   BoundaryCondition bcVelocityR;
   BoundaryCondition bcVelocityThetaPhi;

   BoundaryCondition bcVelocity;
   BoundaryCondition bcTemperature;

   // Corresponding Boundary UIDs for velocity and temperature to setup

   BoundaryUID idSurface;
   BoundaryUID idCMB;

   BoundaryUID idSurfaceT;
   BoundaryUID idCMBT;

   // Timing Analysis SQL db
   std::shared_ptr< FixedSizeSQLDB > db;

   // function and solver initialization

   enum class BCs_T
   {
      VELOCITY_BC,
      VELOCITYROTATION_BC,
      TEMPERATURE_BC,
      NO_BC
   };

   enum class Level_T
   {
      MINLEVEL,
      MAXLEVEL,
      MAXLEVELPLUSONE
   };

   /** Initialise functions
     * 
     * To create a new function, just add into the correct dict and it can be used from the container anywhere in the app
     * The tuple denotes < 
     * std::string           = (Name of the function),
     * uint_t                = (Minimum level),
     * uint_t                = (Maximum level),
     * BCs_T = (Velocity or Temperature or No Boundary condition)
     * >
     */
   //
   // TODO: Implement a seperate FE Function container for TerraNeo when the crashes get annoying
   //
   std::vector< std::tuple< std::string, Level_T, Level_T, BCs_T > > p2p1StokesFunctionDict = {
       { "VelocityFE", Level_T::MINLEVEL, Level_T::MAXLEVEL, BCs_T::VELOCITY_BC },
       { "VelocityFEPrev", Level_T::MINLEVEL, Level_T::MAXLEVEL, BCs_T::VELOCITY_BC },
       { "StokesRHS", Level_T::MINLEVEL, Level_T::MAXLEVEL, BCs_T::VELOCITY_BC },
       { "VelocityFERotated", Level_T::MINLEVEL, Level_T::MAXLEVEL, BCs_T::VELOCITYROTATION_BC },
       { "StokesRHSRotated", Level_T::MINLEVEL, Level_T::MAXLEVEL, BCs_T::VELOCITYROTATION_BC },
       { "StokesTmp1", Level_T::MINLEVEL, Level_T::MAXLEVEL, BCs_T::VELOCITY_BC },
       { "StokesTmp2", Level_T::MINLEVEL, Level_T::MAXLEVEL, BCs_T::VELOCITY_BC },
       { "StokesTmp3", Level_T::MINLEVEL, Level_T::MAXLEVEL, BCs_T::VELOCITY_BC },
       { "StokesResidual", Level_T::MINLEVEL, Level_T::MAXLEVEL, BCs_T::VELOCITY_BC },
       { "StokesTmpProlongation", Level_T::MINLEVEL, Level_T::MAXLEVEL, BCs_T::VELOCITY_BC } };
   std::map< std::string, std::shared_ptr< P2P1StokesFunction_T > > p2p1StokesFunctionContainer;

   std::vector< std::tuple< std::string, Level_T, Level_T, BCs_T > > p2ScalarFunctionDict = {
       { "TemperatureFE", Level_T::MINLEVEL, Level_T::MAXLEVEL, BCs_T::TEMPERATURE_BC },
       { "TemperaturePrev", Level_T::MINLEVEL, Level_T::MAXLEVEL, BCs_T::TEMPERATURE_BC },
       { "TemperatureDev", Level_T::MINLEVEL, Level_T::MAXLEVEL, BCs_T::TEMPERATURE_BC },
       { "Temperature[K]", Level_T::MINLEVEL, Level_T::MAXLEVEL, BCs_T::TEMPERATURE_BC },
       { "ReferenceTemperature", Level_T::MINLEVEL, Level_T::MAXLEVEL, BCs_T::TEMPERATURE_BC },
       { "ViscosityFE", Level_T::MINLEVEL, Level_T::MAXLEVEL, BCs_T::NO_BC },
       { "ViscosityFEInv", Level_T::MINLEVEL, Level_T::MAXLEVEL, BCs_T::NO_BC },
       { "Viscosity[Pas]", Level_T::MINLEVEL, Level_T::MAXLEVEL, BCs_T::NO_BC },
       { "EnergyRHSWeak", Level_T::MINLEVEL, Level_T::MAXLEVEL, BCs_T::TEMPERATURE_BC },
       { "DensityFE", Level_T::MINLEVEL, Level_T::MAXLEVEL, BCs_T::NO_BC },
       { "ShearHeatingTermCoeff", Level_T::MINLEVEL, Level_T::MAXLEVEL, BCs_T::NO_BC },
       { "ShearHeatingTermCoeffDebug", Level_T::MINLEVEL, Level_T::MAXLEVEL, BCs_T::NO_BC },
       { "VelocityMagnitudeSquared", Level_T::MINLEVEL, Level_T::MAXLEVEL, BCs_T::NO_BC } };
   std::map< std::string, std::shared_ptr< P2ScalarFunction_T > > p2ScalarFunctionContainer;

   std::vector< std::tuple< std::string, Level_T, Level_T, BCs_T > > p2VectorFunctionDict = {
       { "NormalsFS", Level_T::MINLEVEL, Level_T::MAXLEVEL, BCs_T::VELOCITY_BC } };
   std::map< std::string, std::shared_ptr< P2VectorFunction_T > > p2VectorFunctionContainer;

   std::vector< std::tuple< std::string, Level_T, Level_T, BCs_T > > p1ScalarFunctionDict = {
       { "ViscosityFEP1", Level_T::MINLEVEL, Level_T::MAXLEVELPLUSONE, BCs_T::NO_BC },
       { "ShearHeatingTermCoeffP1", Level_T::MINLEVEL, Level_T::MAXLEVELPLUSONE, BCs_T::NO_BC },
       { "TemperatureFEP1", Level_T::MINLEVEL, Level_T::MAXLEVELPLUSONE, BCs_T::TEMPERATURE_BC },
       { "TemperatureFEPrevP1", Level_T::MINLEVEL, Level_T::MAXLEVELPLUSONE, BCs_T::TEMPERATURE_BC },
       { "EnergyRHSP1", Level_T::MINLEVEL, Level_T::MAXLEVELPLUSONE, BCs_T::TEMPERATURE_BC } };
   std::map< std::string, std::shared_ptr< P1ScalarFunction_T > > p1ScalarFunctionContainer;

   std::vector< std::tuple< std::string, Level_T, Level_T, BCs_T > > p1VectorFunctionDict = {
       { "VelocityFEP1", Level_T::MINLEVEL, Level_T::MAXLEVELPLUSONE, BCs_T::VELOCITY_BC },
       { "VelocityFEPrevP1", Level_T::MINLEVEL, Level_T::MAXLEVELPLUSONE, BCs_T::VELOCITY_BC } };
   std::map< std::string, std::shared_ptr< P1VectorFunction_T > > p1VectorFunctionContainer;

   std::vector< std::tuple< std::string, Level_T, Level_T, BCs_T > > p0ScalarFunctionDict = {
       { "ViscosityFEP0", Level_T::MINLEVEL, Level_T::MAXLEVELPLUSONE, BCs_T::NO_BC } };
   std::map< std::string, std::shared_ptr< P0ScalarFunction_T > > p0ScalarFunctionContainer;

   // Storage for primitives (includes functionality for distributed computing)
   std::shared_ptr< PrimitiveStorage > storage;

   // Solvers
   std::shared_ptr< MCSolverBase< StokesOperator_T > > stokesMCSolver_;
   std::shared_ptr< Solver< StokesOperator_T > >       stokesSolver_;

   std::shared_ptr< MCSolverBase< StokesOperatorRotationOpgen_T > > stokesRotationOpgenMCSolver_;
   std::shared_ptr< Solver< StokesOperatorRotationOpgen_T > >       stokesRotationOpgenSolver_;

   std::shared_ptr< CGSolver< TransportOperator_T > > temperatureTransportSolver_;

   // Operators
   std::shared_ptr< StokesOperator_T >              stokesOperator_;
   std::shared_ptr< StokesOperatorRotationOpgen_T > stokesOperatorRotationOpgen_;

   std::shared_ptr< MMOCTransport< P1Function< real_t > > > temperatureMMOCOperatorP1_;

   std::shared_ptr< MMOCTransport< TemperatureFunction_T > > temperatureMMOCOperator_;
   std::shared_ptr< TransportOperator_T >                    temperatureTransportOperator_;

   std::shared_ptr< P2ScalarMassOperator_T > p2ScalarMassOperator_;

   std::shared_ptr< ProjectionOperator_T > projectionOperator_;
   std::shared_ptr< RotationOperator_T >   rotationOperator_;

   std::shared_ptr< P2toP2QuadraticProlongation > p2ProlongationOperator_;
   std::shared_ptr< P2toP2QuadraticInjection >    p2InjectionOperator_;
   std::shared_ptr< FrozenVelocityFullOperator >  frozenVelocityRHS_;

   bool outputDirectoriesCreated = false;

   std::string modelPath;
   std::string modelParaviewPath;
   std::string modelCheckpointPath;
   std::string modelRadialProfilesPath;

   std::shared_ptr< hyteg::VTKOutput > vtkOutput_;

   // ADIOS2 data output
#ifdef HYTEG_BUILD_WITH_ADIOS2
   std::shared_ptr< AdiosWriter > adiosOutput_;

   std::shared_ptr< AdiosCheckpointExporter > checkpointExporter;
   std::shared_ptr< AdiosCheckpointImporter > checkpointImporter;

   std::map< std::string, adiosHelpers::adiostype_t > attributeList_;
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

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > shearHeatingCoeffCalc;

   std::function< real_t( const Point3D& ) > normalX = []( const Point3D& x ) { return -x[0] / x.norm(); };
   std::function< real_t( const Point3D& ) > normalY = []( const Point3D& x ) { return -x[1] / x.norm(); };
   std::function< real_t( const Point3D& ) > normalZ = []( const Point3D& x ) { return -x[2] / x.norm(); };

   std::shared_ptr< RadialProfile >                       viscosityProfiles;
   std::shared_ptr< RadialProfile >                       temperatureProfiles;
   std::shared_ptr< RadialProfile >                       velocityProfiles;
   std::shared_ptr< TemperatureInitializationParameters > temperatureInitParams;

   std::function< real_t( const Point3D& ) > randFunc;
   std::function< real_t( const Point3D& ) > temperatureInitialCondition;
};

} // namespace terraneo
