/*
 * Copyright (c) 2017-2020 Nils Kohl.
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

#include <algorithm>
#include <cmath>
#include <core/Environment.h>
#include <core/math/Constants.h>

#include "core/DataTypes.h"
#include "core/config/Config.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/MeshQuality.hpp"
#include "hyteg/composites/StrongFreeSlipWrapper.hpp"
#include "hyteg/composites/UnsteadyDiffusion.hpp"
#include "hyteg/dataexport/SQL.hpp"
#include "hyteg/dataexport/TimingOutput.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/geometry/SphereTools.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/memory/MemoryAllocation.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/numerictools/CFDHelpers.hpp"
#include "hyteg/numerictools/SphericalHarmonicsTool.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/petsc/PETScBlockPreconditionedStokesSolver.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/FullMultigridSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"
#include "hyteg/solvers/controlflow/SolverLoop.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesBlockDiagonalPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesVelocityBlockBlockDiagonalPreconditioner.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"

#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"
#include "sqlite/SQLite.h"

namespace hyteg {

using walberla::int_c;
using walberla::real_t;
using walberla::math::pi;

#if 0
static std::string getDateTimeID()
{
   std::vector< char > cTimeString( 64 );
   WALBERLA_ROOT_SECTION()
   {
      std::time_t t;
      std::time( &t );
      std::strftime( cTimeString.data(), 64, "%F_%H-%M-%S", std::localtime( &t ) );
   }

   walberla::mpi::broadcastObject( cTimeString );

   std::string timeString( cTimeString.data() );
   return timeString;
}
#endif

enum class StokesSolverType : int
{
   PETSC_MUMPS         = 0,
   PETSC_MINRES_JACOBI = 1,
   PETSC_MINRES_BOOMER = 2,
   HYTEG_MINRES        = 3,
   HYTEG_MINRES_GMG    = 4,
   HYTEG_UZAWA_V       = 5,
   HYTEG_UZAWA_FMG     = 6
};

enum class DiffusionSolverType : int
{
   PETSC_MINRES = 0,
   HYTEG_CG     = 1,
   HYTEG_GMG    = 2
};

struct DomainInfo
{
   bool threeDim = false;

   real_t rMin = 0;
   real_t rMax = 0;
   uint_t nTan = 0;
   uint_t nRad = 0;

   real_t domainVolume() const
   {
      if ( !threeDim )
      {
         return walberla::math::pi * rMax * rMax - walberla::math::pi * rMin * rMin;
      }
      else
      {
         return ( 4. / 3. ) * walberla::math::pi * rMax * rMax * rMax - ( 4. / 3. ) * walberla::math::pi * rMin * rMin * rMin;
      }
   }
};

struct SolverInfo
{
   StokesSolverType stokesSolverType = StokesSolverType::PETSC_MUMPS;

   uint_t stokesMaxNumIterations           = 10;
   real_t stokesAbsoluteResidualUTolerance = 0;
   real_t stokesRelativeResidualUTolerance = 0;
   uint_t uzawaInnerIterations             = 10;
   uint_t uzawaPreSmooth                   = 6;
   uint_t uzawaPostSmooth                  = 6;

   DiffusionSolverType diffusionSolverType = DiffusionSolverType::PETSC_MINRES;

   uint_t diffusionMaxNumIterations           = 10000;
   real_t diffusionAbsoluteResidualUTolerance = 10000;

   real_t uzawaOmega = 0.3;
};

struct OutputInfo
{
   bool        vtk;
   bool        vtkOutputVelocity;
   bool        sphericalTemperatureSlice;
   int         sphericalTemperatureSliceNumMeridians;
   int         sphericalTemperatureSliceNumParallels;
   int         sphericalTemperatureSliceIcoRefinements;
   uint_t      printInterval;
   uint_t      vtkInterval;
   uint_t      sphericalTemperatureSliceInterval;
   bool        vtkOutputVertexDoFs;
   std::string outputDirectory;
   std::string outputBaseName;
};

/// Calculates and returns
///
///     ||u||_L2 = sqrt( u^T M u )
///
template < typename FunctionType, typename MassOperator >
real_t normL2( const FunctionType& u, const FunctionType& tmp, const MassOperator& M, const uint_t& level, const DoFType& flag )
{
   tmp.interpolate( 0, level );
   M.apply( u, tmp, level, flag );
   return std::sqrt( u.dotGlobal( tmp, level, flag ) );
}

template < typename FunctionType, typename MassOperator >
real_t normL2Squared( const FunctionType& u,
                      const FunctionType& tmp,
                      const MassOperator& M,
                      const uint_t&       level,
                      const DoFType&      flag )
{
   tmp.interpolate( 0, level );
   M.apply( u, tmp, level, flag );
   return u.dotGlobal( tmp, level, flag );
}

template < typename StokesFunction, typename StokesOperator, typename MassOperatorVelocity, typename MassOperatorPressure >
void calculateStokesResiduals( const StokesOperator&       A,
                               const MassOperatorVelocity& Mu,
                               const MassOperatorPressure& Mp,
                               const StokesFunction&       x,
                               const StokesFunction&       b,
                               uint_t                      level,
                               StokesFunction&             r,
                               StokesFunction&             tmp,
                               real_t&                     residual )
{
   tmp.interpolate( 0, level, All );
   r.interpolate( 0, level, All );
   A.apply( x, tmp, level, Inner | NeumannBoundary | FreeslipBoundary );
   r.assign( { 1.0, -1.0 }, { b, tmp }, level, Inner | NeumannBoundary | FreeslipBoundary );
   residual = normL2Squared( r.uvw()[0], tmp.uvw()[0], Mu, level, Inner | NeumannBoundary | FreeslipBoundary );
   residual += normL2Squared( r.uvw()[1], tmp.uvw()[1], Mu, level, Inner | NeumannBoundary | FreeslipBoundary );
   if ( x.getStorage()->hasGlobalCells() )
   {
      residual += normL2Squared( r.uvw()[2], tmp.uvw()[2], Mu, level, Inner | NeumannBoundary | FreeslipBoundary );
   }
   residual += normL2Squared( r.p(), tmp.p(), Mp, level, Inner | NeumannBoundary | FreeslipBoundary );
   residual = std::sqrt( residual );
}

template < typename UnsteadyDiffusion,
           typename UnsteadyDiffusionOperator,
           typename LaplaceOperator,
           typename MassOperatorVelocity,
           typename ScalarFuncionType >
void calculateDiffusionResidual( UnsteadyDiffusion&               unsteadyDiffusion,
                                 const UnsteadyDiffusionOperator& unsteadyDiffusionOperator,
                                 const LaplaceOperator&           laplacian,
                                 const MassOperatorVelocity&      Mu,
                                 const ScalarFuncionType&         c,
                                 const ScalarFuncionType&         cOld,
                                 const ScalarFuncionType&         f,
                                 ScalarFuncionType&               r,
                                 ScalarFuncionType&               tmp,
                                 uint_t                           level,
                                 real_t&                          residual )
{
   unsteadyDiffusion.calculateResidual(
       unsteadyDiffusionOperator, laplacian, Mu, c, cOld, f, f, r, level, Inner | NeumannBoundary | FreeslipBoundary );
   residual = normL2( r, tmp, Mu, level, Inner | NeumannBoundary | FreeslipBoundary );
}

template < typename StokesFunction, typename VelocityMass >
real_t velocityRMS( const StokesFunction& u, const StokesFunction& tmp, const VelocityMass& M, real_t domainVolume, uint_t level )
{
   auto norm = std::pow( normL2( u.uvw()[0], tmp.uvw()[0], M, level, All ), 2.0 ) +
               std::pow( normL2( u.uvw()[1], tmp.uvw()[1], M, level, All ), 2.0 );

   if ( u.getStorage()->hasGlobalCells() )
   {
      norm += std::pow( normL2( u.uvw()[2], tmp.uvw()[2], M, level, All ), 2.0 );
   }

   return std::sqrt( norm / domainVolume );
}

std::shared_ptr< SetupPrimitiveStorage > createSetupStorage( DomainInfo domainInfo )
{
   MeshInfo meshInfo = MeshInfo::emptyMeshInfo();
   if ( domainInfo.threeDim )
   {
      meshInfo = MeshInfo::meshSphericalShell( domainInfo.nTan, domainInfo.nRad, domainInfo.rMin, domainInfo.rMax );
   }
   else
   {
      meshInfo = MeshInfo::meshAnnulus( domainInfo.rMin, domainInfo.rMax, MeshInfo::CRISS, domainInfo.nTan, domainInfo.nRad );
   }

   auto setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   loadbalancing::roundRobinVolume( *setupStorage );

   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   if ( domainInfo.threeDim )
   {
      IcosahedralShellMap::setMap( *setupStorage );
   }
   else
   {
      AnnulusMap::setMap( *setupStorage );
   }

   return setupStorage;
}

void runBenchmark( real_t     cflMax,
                   real_t     rayleighNumber,
                   bool       fixedTimeStep,
                   real_t     dtConstant,
                   uint_t     maxNumTimeSteps,
                   uint_t     minLevel,
                   uint_t     level,
                   SolverInfo solverInfo,
                   DomainInfo domainInfo,
                   bool       predictorCorrector,
                   real_t     simulationTime,
                   OutputInfo outputInfo,
                   bool       verbose )
{
   walberla::WcTimer localTimer;

   const bool outputTimingJSON = true;

   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Building setup storage ..." );
   }

   auto setupStorage = createSetupStorage( domainInfo );

   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Building distributed storage ..." );
   }

   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage, 1 );

   if ( outputInfo.vtk )
   {
      if ( verbose )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Writing domain partitioning VTK ..." );
      }

      writeDomainPartitioningVTK( storage, outputInfo.outputDirectory, outputInfo.outputBaseName + "_domain" );
   }

   auto timer = storage->getTimingTree();

   timer->start( "Setup" );

   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Gathering domain information ..." );
   }

   const uint_t unknownsTemperature = numberOfGlobalDoFs< P2FunctionTag >( *storage, level );
   const uint_t unknownsStokes      = numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, level );
   const real_t hMin                = MeshQuality::getMinimalEdgeLength( storage, level );
   const real_t hMax                = MeshQuality::getMaximalEdgeLength( storage, level );

   const real_t diffusivity     = 1.0;
   const real_t internalHeating = 0.0;

   const real_t sliceEvaluationRadius = domainInfo.rMin + 0.5 * ( domainInfo.rMax - domainInfo.rMin );

   WALBERLA_LOG_INFO_ON_ROOT( "" )
   WALBERLA_LOG_INFO_ON_ROOT( "Benchmark name: " << outputInfo.outputBaseName )
   WALBERLA_LOG_INFO_ON_ROOT( " - time discretization: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + fixed time step:                              " << fixedTimeStep )
   if ( fixedTimeStep )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "   + dt (constant):                                " << dtConstant )
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( "   + max CFL:                                      " << cflMax )
   }
   WALBERLA_LOG_INFO_ON_ROOT( "   + simulation time:                              " << simulationTime )
   WALBERLA_LOG_INFO_ON_ROOT( " - space discretization: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + rMin:                                         " << domainInfo.rMin )
   WALBERLA_LOG_INFO_ON_ROOT( "   + rMax:                                         " << domainInfo.rMax )
   WALBERLA_LOG_INFO_ON_ROOT( "   + nTan:                                         " << domainInfo.nTan )
   WALBERLA_LOG_INFO_ON_ROOT( "   + nRad:                                         " << domainInfo.nRad )
   WALBERLA_LOG_INFO_ON_ROOT( "   + dimensions:                                   " << ( storage->hasGlobalCells() ? "3" : "2" ) )
   WALBERLA_LOG_INFO_ON_ROOT( "   + level:                                        " << level )
   WALBERLA_LOG_INFO_ON_ROOT( "   + unknowns temperature, including boundary:     " << unknownsTemperature )
   for ( uint_t l = minLevel; l <= level; l++ )
   {
      const uint_t unknownsStokesLevel = numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, l );
      WALBERLA_LOG_INFO_ON_ROOT( "   + unknowns Stokes, including boundary, level " << l << ": " << unknownsStokesLevel )
   }
   WALBERLA_LOG_INFO_ON_ROOT( "   + h_min:                                        " << hMin )
   WALBERLA_LOG_INFO_ON_ROOT( "   + h_max:                                        " << hMax )
   WALBERLA_LOG_INFO_ON_ROOT( " - benchmark settings: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + Rayleigh number:                              " << rayleighNumber )
   WALBERLA_LOG_INFO_ON_ROOT( "   + internal heating:                             " << internalHeating )
   WALBERLA_LOG_INFO_ON_ROOT( " - app settings: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + Stokes solver:                                " << int_c( solverInfo.stokesSolverType ) )
   WALBERLA_LOG_INFO_ON_ROOT( "   + Uzawa omega:                                  " << solverInfo.uzawaOmega )
   WALBERLA_LOG_INFO_ON_ROOT( "   + diffusion solver:                             " << int_c( solverInfo.diffusionSolverType ) )

   if ( outputInfo.vtk )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "   + VTK interval:                                 " << outputInfo.vtkInterval )
   }
   if ( outputInfo.sphericalTemperatureSlice )
   {
      WALBERLA_LOG_INFO_ON_ROOT(
          "   + spherical slice output interval:              " << outputInfo.sphericalTemperatureSliceInterval )
   }
   WALBERLA_LOG_INFO_ON_ROOT( "   + print interval:                               " << outputInfo.printInterval )
   WALBERLA_LOG_INFO_ON_ROOT( "   + output directory:                             " << outputInfo.outputDirectory )
   WALBERLA_LOG_INFO_ON_ROOT( "   + output base name:                             " << outputInfo.outputBaseName )
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   auto storageInfo = storage->getGlobalInfo();
   WALBERLA_LOG_INFO_ON_ROOT( storageInfo );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   FixedSizeSQLDB db( outputInfo.outputDirectory + "/" + outputInfo.outputBaseName + ".db" );

   db.setConstantEntry( "ra", rayleighNumber );
   db.setConstantEntry( "ntan", domainInfo.nTan );
   db.setConstantEntry( "nrad", domainInfo.nRad );
   db.setConstantEntry( "rmin", domainInfo.rMin );
   db.setConstantEntry( "rmax", domainInfo.rMax );
   db.setConstantEntry( "fixed_time_step", fixedTimeStep );
   db.setConstantEntry( "cfl_max", cflMax );
   db.setConstantEntry( "dt_constant", dtConstant );
   db.setConstantEntry( "simulation_time", simulationTime );
   db.setConstantEntry( "level", uint_c( level ) );
   db.setConstantEntry( "unknowns_temperature", uint_c( unknownsTemperature ) );
   db.setConstantEntry( "unknowns_stokes", uint_c( unknownsStokes ) );
   db.setConstantEntry( "h_min", hMin );
   db.setConstantEntry( "h_max", hMax );
   db.setConstantEntry( "num_macro_cells", uint_c( setupStorage->getNumberOfCells() ) );
   db.setConstantEntry( "num_macro_faces", uint_c( setupStorage->getNumberOfFaces() ) );
   db.setConstantEntry( "num_macro_edges", uint_c( setupStorage->getNumberOfEdges() ) );
   db.setConstantEntry( "num_macro_vertices", uint_c( setupStorage->getNumberOfVertices() ) );
   db.setConstantEntry( "num_macro_primitives", uint_c( setupStorage->getNumberOfPrimitives() ) );
   db.setConstantEntry( "diffusivity", diffusivity );

   // db.setConstantEntry( "stokes_solver_type", uint_c( solverInfo.stokesSolverType ) );
   db.setConstantEntry( "uzawa_omega", solverInfo.uzawaOmega );
   db.setConstantEntry( "uzawa_inner_smooth", solverInfo.uzawaInnerIterations );
   db.setConstantEntry( "uzawa_pre_smooth", solverInfo.uzawaPreSmooth );
   db.setConstantEntry( "uzawa_post_smooth", solverInfo.uzawaPostSmooth );

   db.setConstantEntry( "stokes_absolute_residual_threshold", solverInfo.stokesAbsoluteResidualUTolerance );
   db.setConstantEntry( "stokes_relative_residual_threshold", solverInfo.stokesRelativeResidualUTolerance );

   typedef P2P1TaylorHoodFunction< real_t >       StokesFunction;
   typedef P2Function< real_t >                   ScalarFunction;
   typedef P2P1ElementwiseBlendingStokesOperator  StokesOperator;
   typedef P2ElementwiseBlendingLaplaceOperator   LaplaceOperator;
   typedef P2ElementwiseBlendingMassOperator      MassOperatorVelocity;
   typedef P1ConstantMassOperator                 MassOperatorPressure;
   typedef P2ElementwiseUnsteadyDiffusionOperator UnsteadyDiffusionOperator;

   BoundaryCondition bcVelocity;
   BoundaryCondition bcTemperature;

   bcVelocity.createDirichletBC( "outerBoundaryVelocity", 1 );
   bcVelocity.createDirichletBC( "innerBoundaryVelocity", 1 );
   bcTemperature.createDirichletBC( "outerBoundaryTemperature", 1 );
   bcTemperature.createDirichletBC( "innerBoundaryTemperature", 1 );

   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Allocating functions ..." );
   }

   ScalarFunction c( "c", storage, minLevel, level, bcTemperature );
   ScalarFunction cPr( "cPr", storage, minLevel, level, bcTemperature );
   ScalarFunction cOld( "cOld", storage, minLevel, level, bcTemperature );
   ScalarFunction cTmp( "cTmp", storage, minLevel, level, bcTemperature );
   ScalarFunction cTmp2( "cTmp2", storage, minLevel, level, bcTemperature );
   ScalarFunction q( "q", storage, minLevel, level, bcTemperature );

   StokesFunction u( "u", storage, minLevel, level, bcVelocity );
   StokesFunction uLast( "uLast", storage, minLevel, level, bcVelocity );
   StokesFunction f( "f", storage, minLevel, level, bcVelocity );
   StokesFunction outwardNormal( "outwardNormal", storage, minLevel, level, bcVelocity );
   StokesFunction stokesTmp( "stokesTmp", storage, minLevel, level, bcVelocity );
   StokesFunction stokesResidual( "stokesResidual", storage, minLevel, level, bcVelocity );

   ScalarFunction uMagnitudeSquared( "uMagnitudeSquared", storage, minLevel, level, bcVelocity );

   ScalarFunction uTmp( "uTmp", storage, minLevel, level, bcVelocity );
   ScalarFunction uTmp2( "uTmp2", storage, minLevel, level, bcVelocity );

   SphericalHarmonicsTool sphTool( 6 );

   std::function< real_t( const Point3D& ) > initialTemperature = [&]( const Point3D& x ) {
      auto radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );

      const auto linearSlope =
          std::max( real_c( 0 ), 1 - ( ( radius - domainInfo.rMin ) / ( domainInfo.rMax - domainInfo.rMin ) ) );

      const auto xn = x / x.norm();

      // sph ~ in (-2, 2)
      const auto powerPertubation = sphTool.shconvert_eval( 4, 3, xn[0], xn[1], xn[2] );

      const auto powerSlope = std::pow( linearSlope, 4 + powerPertubation );
      return powerSlope;
   };

   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Interpolating initial temperature and internal heating ..." );
   }

   for ( uint_t l = minLevel; l <= level; l++ )
   {
      c.interpolate( initialTemperature, l, All );
      q.interpolate( internalHeating, l, All );
   }

   std::function< real_t( const Point3D& ) > normalX = []( const Point3D& x ) { return x[0] / x.norm(); };
   std::function< real_t( const Point3D& ) > normalY = []( const Point3D& x ) { return x[1] / x.norm(); };
   std::function< real_t( const Point3D& ) > normalZ = []( const Point3D& x ) { return x[2] / x.norm(); };

   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Interpolating normals ..." );
   }

   for ( uint_t l = minLevel; l <= level; l++ )
   {
      outwardNormal.uvw().interpolate( { normalX, normalY, normalZ }, l );
   }

   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Preparing operators ..." );
   }

   UnsteadyDiffusionOperator diffusionOperator(
       storage, minLevel, level, 1.0, diffusivity, DiffusionTimeIntegrator::ImplicitEuler );
   auto                            A = std::make_shared< StokesOperator >( storage, minLevel, level );
   LaplaceOperator                 L( storage, minLevel, level );
   MassOperatorVelocity            MVelocity( storage, minLevel, level );
   MassOperatorPressure            MPressure( storage, minLevel, level );
   MMOCTransport< ScalarFunction > transport( storage, minLevel, level, TimeSteppingScheme::RK4 );

   if ( storage->hasGlobalCells() )
   {
      A->computeAndStoreLocalElementMatrices();
      L.computeAndStoreLocalElementMatrices();
      MVelocity.computeAndStoreLocalElementMatrices();
   }

   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Preparing solvers ..." );
   }

   std::shared_ptr< Solver< StokesOperator > > stokesSolver;
   std::shared_ptr< Solver< StokesOperator > > stokesSolverBlockPrecMinRes;

   real_t initialResiudalU          = 0;
   real_t vCycleResidualULast       = 0;
   real_t vCycleResidualCLast       = 0;
   uint_t numVCycles                = 0;
   real_t averageResidualReductionU = 0;

   if ( solverInfo.stokesSolverType == StokesSolverType::PETSC_MUMPS )
   {
      stokesSolver = std::make_shared< PETScLUSolver< StokesOperator > >( storage, level );
   }

   else if ( solverInfo.stokesSolverType == StokesSolverType::PETSC_MINRES_JACOBI ||
             solverInfo.stokesSolverType == StokesSolverType::PETSC_MINRES_BOOMER )
   {
      auto velocityPreconditionerType = 1;
      if ( solverInfo.stokesSolverType == StokesSolverType::PETSC_MINRES_BOOMER )
      {
         velocityPreconditionerType = 3;
      }
      auto stokesSolverTmp =
          std::make_shared< PETScBlockPreconditionedStokesSolver< StokesOperator > >( storage,
                                                                                      level,
                                                                                      solverInfo.stokesAbsoluteResidualUTolerance,
                                                                                      solverInfo.stokesMaxNumIterations,
                                                                                      velocityPreconditionerType,
                                                                                      1 );
      stokesSolverTmp->reassembleMatrix( false );
      stokesSolverTmp->setVerbose( verbose );
      stokesSolver = stokesSolverTmp;
   }

   else if ( solverInfo.stokesSolverType == StokesSolverType::HYTEG_MINRES )
   {
      stokesSolver = solvertemplates::stokesMinResSolver< StokesOperator >(
          storage, level, solverInfo.stokesAbsoluteResidualUTolerance, solverInfo.stokesMaxNumIterations, verbose );
   }

   else if ( solverInfo.stokesSolverType == StokesSolverType::HYTEG_MINRES_GMG )
   {
      typedef CGSolver< P2ElementwiseBlendingLaplaceOperator >                 CoarseGridSolver_T;
      typedef GeometricMultigridSolver< P2ElementwiseBlendingLaplaceOperator > GMGSolver_T;

      auto coarseGridSolver = std::make_shared< CoarseGridSolver_T >( storage, minLevel, level );
      auto smoother =
          std::make_shared< WeightedJacobiSmoother< P2ElementwiseBlendingLaplaceOperator > >( storage, minLevel, level, 0.66 );
      auto prolongationOperator = std::make_shared< P2toP2QuadraticProlongation >();
      auto restrictionOperator  = std::make_shared< P2toP2QuadraticRestriction >();
      auto gmgSolver            = std::make_shared< GMGSolver_T >(
          storage, smoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, level, 3, 3 );

      auto lumpedInvPressureMass = std::make_shared< P1LumpedInvMassOperator >( storage, minLevel, level );

      typedef StokesBlockDiagonalPreconditioner< StokesOperator, P1LumpedInvMassOperator > Preconditioner_T;
      auto preconditioner = std::make_shared< Preconditioner_T >( storage, minLevel, level, 3, gmgSolver );

      auto stokesSolverTmp = std::make_shared< MinResSolver< StokesOperator > >(
          storage, minLevel, level, solverInfo.stokesMaxNumIterations, 1e-30, preconditioner );
      stokesSolverTmp->setPrintInfo( verbose );
      stokesSolverTmp->setAbsoluteTolerance( solverInfo.stokesAbsoluteResidualUTolerance );
      if ( solverInfo.stokesSolverType == StokesSolverType::HYTEG_MINRES_GMG )
      {
         stokesSolver = stokesSolverTmp;
      }
      else
      {
         stokesSolverBlockPrecMinRes = stokesSolverTmp;
      }
   }

   else if ( solverInfo.stokesSolverType == StokesSolverType::HYTEG_UZAWA_V ||
             solverInfo.stokesSolverType == StokesSolverType::HYTEG_UZAWA_FMG )
   {
      auto prolongationOperator = std::make_shared< P2P1StokesToP2P1StokesProlongation >();
      auto restrictionOperator  = std::make_shared< P2P1StokesToP2P1StokesRestriction >( true );

      std::shared_ptr< Solver< typename StokesOperator::VelocityOperator_T > > smoother =
          std::make_shared< WeightedJacobiSmoother< typename StokesOperator::VelocityOperator_T > >(
              storage, minLevel, level, 0.66 );

      auto uzawaVelocityPreconditioner =
          std::make_shared< StokesVelocityBlockBlockDiagonalPreconditioner< StokesOperator > >( storage, smoother );

      auto uzawaSmoother = std::make_shared< UzawaSmoother< StokesOperator > >( storage,
                                                                                uzawaVelocityPreconditioner,
                                                                                minLevel,
                                                                                level,
                                                                                solverInfo.uzawaOmega,
                                                                                Inner | NeumannBoundary,
                                                                                solverInfo.uzawaInnerIterations );

      std::shared_ptr< Solver< StokesOperator > > coarseGridSolverInternal;

      auto petscSolverInternalTmp = std::make_shared< PETScBlockPreconditionedStokesSolver< StokesOperator > >(
          storage, minLevel, solverInfo.stokesAbsoluteResidualUTolerance, 1000, 1 );
      petscSolverInternalTmp->setVerbose( verbose );
      auto coarseGridSolver = petscSolverInternalTmp;

      auto multigridSolver = std::make_shared< GeometricMultigridSolver< StokesOperator > >( storage,
                                                                                             uzawaSmoother,
                                                                                             coarseGridSolver,
                                                                                             restrictionOperator,
                                                                                             prolongationOperator,
                                                                                             minLevel,
                                                                                             level,
                                                                                             solverInfo.uzawaPreSmooth,
                                                                                             solverInfo.uzawaPostSmooth,
                                                                                             2,
                                                                                             CycleType::VCYCLE );

      if ( solverInfo.stokesSolverType == StokesSolverType::HYTEG_UZAWA_V )
      {
         auto stopIterationCallback =
             [&]( const StokesOperator& _A, const StokesFunction& _u, const StokesFunction& _b, const uint_t _level ) {
                real_t r_u;

                calculateStokesResiduals( _A, MVelocity, MPressure, _u, _b, _level, stokesResidual, stokesTmp, r_u );

                if ( numVCycles == 0 )
                {
                   WALBERLA_LOG_INFO_ON_ROOT(
                       walberla::format( "[Uzawa] iter %3d | residual: %10.5e | initial ", 0, vCycleResidualULast ) );
                }

                auto reductionRateU = r_u / vCycleResidualULast;

                vCycleResidualULast = r_u;

                numVCycles++;
                averageResidualReductionU += reductionRateU;

                if ( verbose )
                {
                   WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
                       "[Uzawa] iter %3d | residual: %10.5e | reduction: %10.5e ", numVCycles, r_u, reductionRateU ) );
                }

                if ( r_u / initialResiudalU < solverInfo.stokesRelativeResidualUTolerance )
                {
                   WALBERLA_LOG_INFO_ON_ROOT( "[Uzawa] reached relative residual threshold" )
                   return true;
                }

                if ( r_u < solverInfo.stokesAbsoluteResidualUTolerance )
                {
                   WALBERLA_LOG_INFO_ON_ROOT( "[Uzawa] reached absolute residual threshold" )
                   return true;
                }

                if ( reductionRateU > 0.8 )
                {
                   WALBERLA_LOG_INFO_ON_ROOT( "[Uzawa] reached convergence rate threshold" )
                   return true;
                }

                return false;
             };

         stokesSolver = std::make_shared< SolverLoop< StokesOperator > >(
             multigridSolver, solverInfo.stokesMaxNumIterations, stopIterationCallback );
      }
      else
      {
         auto fmgSolver = std::make_shared< FullMultigridSolver< StokesOperator > >(
             storage, multigridSolver, prolongationOperator, minLevel, level, 3 );
         stokesSolver = fmgSolver;
      }
   }
   else
   {
      WALBERLA_ABORT( "Invalid solver type." );
   }

   std::shared_ptr< Solver< UnsteadyDiffusionOperator > > diffusionLinearSolver;

   UnsteadyDiffusion< ScalarFunction, UnsteadyDiffusionOperator, LaplaceOperator, MassOperatorVelocity > diffusionSolver(
       storage, minLevel, level, diffusionLinearSolver );

   if ( solverInfo.diffusionSolverType == DiffusionSolverType::PETSC_MINRES )
   {
      auto internalDiffusionSolver = std::make_shared< PETScMinResSolver< UnsteadyDiffusionOperator > >(
          storage, level, 1e-30, solverInfo.diffusionAbsoluteResidualUTolerance, solverInfo.diffusionMaxNumIterations );
      internalDiffusionSolver->reassembleMatrix( true );
      diffusionLinearSolver = internalDiffusionSolver;
   }
   else if ( solverInfo.diffusionSolverType == DiffusionSolverType::HYTEG_CG )
   {
      auto internalDiffusionSolver = std::make_shared< CGSolver< UnsteadyDiffusionOperator > >(
          storage, minLevel, level, solverInfo.diffusionMaxNumIterations, solverInfo.diffusionAbsoluteResidualUTolerance );
      internalDiffusionSolver->setPrintInfo( verbose );
      diffusionLinearSolver = internalDiffusionSolver;
   }
   else if ( solverInfo.diffusionSolverType == DiffusionSolverType::HYTEG_GMG )
   {
      WALBERLA_ABORT( "Somethings not working here for 3D" );
      typedef CGSolver< UnsteadyDiffusionOperator >                 CoarseGridSolver_T;
      typedef GeometricMultigridSolver< UnsteadyDiffusionOperator > GMGSolver_T;

      auto coarseGridSolver = std::make_shared< CoarseGridSolver_T >( storage, minLevel, level );
      auto smoother = std::make_shared< WeightedJacobiSmoother< UnsteadyDiffusionOperator > >( storage, minLevel, level, 0.66 );
      auto prolongationOperator = std::make_shared< P2toP2QuadraticProlongation >();
      auto restrictionOperator  = std::make_shared< P2toP2QuadraticRestriction >();
      auto gmgSolver            = std::make_shared< GMGSolver_T >(
          storage, smoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, level, 3, 3 );

      auto stopIterationCallback = [&]( const UnsteadyDiffusionOperator&,
                                        const ScalarFunction& _u,
                                        const ScalarFunction& _b,
                                        const uint_t ) {
         real_t r_c;

         calculateDiffusionResidual( diffusionSolver, diffusionOperator, L, MVelocity, _u, cOld, _b, cTmp, cTmp2, level, r_c );

         auto reductionRateC = r_c / vCycleResidualCLast;

         vCycleResidualCLast = r_c;

         if ( verbose )
         {
            WALBERLA_LOG_INFO_ON_ROOT(
                walberla::format( "[Diffusion GMG] residual: %10.5e | reduction: %10.5e", r_c, reductionRateC ) );
         }

         if ( r_c < solverInfo.diffusionAbsoluteResidualUTolerance )
         {
            return true;
         }

         if ( reductionRateC > 0.8 )
         {
            return true;
         }

         return false;
      };

      diffusionLinearSolver = std::make_shared< SolverLoop< UnsteadyDiffusionOperator > >( gmgSolver, 20, stopIterationCallback );
   }

   diffusionSolver.setSolver( diffusionLinearSolver );

   real_t timeTotal = 0;
   real_t vMax      = 0;
   real_t vRms      = 0;
   real_t residualU = 0;

   real_t timeStepTotal = 0;
   real_t timeStokes    = 0;
   real_t timeMMOC      = 0;
   real_t timeDiffusion = 0;
   real_t timeVTK       = 0;

   hyteg::VTKOutput vtkOutput( outputInfo.outputDirectory, outputInfo.outputBaseName, storage, outputInfo.vtkInterval );
   vtkOutput.setVTKDataFormat( vtk::DataFormat::BINARY );

   if ( outputInfo.vtkOutputVertexDoFs )
   {
      vtkOutput.add( c.getVertexDoFFunction() );
   }
   else
   {
      vtkOutput.add( c );
   }

   if ( outputInfo.vtkOutputVelocity )
   {
      if ( outputInfo.vtkOutputVertexDoFs )
      {
         vtkOutput.add( u.uvw()[0].getVertexDoFFunction() );
         vtkOutput.add( u.uvw()[1].getVertexDoFFunction() );
         vtkOutput.add( u.uvw()[2].getVertexDoFFunction() );
      }
      else
      {
         vtkOutput.add( u );
      }
   }
   else
   {
      if ( outputInfo.vtkOutputVertexDoFs )
      {
         vtkOutput.add( uMagnitudeSquared.getVertexDoFFunction() );
      }
      else
      {
         vtkOutput.add( uMagnitudeSquared );
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   printFunctionAllocationInfo( *storage, 1 );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   printCurrentMemoryUsage();
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   timer->stop( "Setup" );

   timer->start( "Simulation" );

   uint_t            timeStep = 0;
   walberla::WcTimer timeStepTimer;

   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Initial Stokes solve ..." );
   }

   timeStepTimer.start();

   for ( uint_t l = minLevel; l <= level; l++ )
   {
      MVelocity.apply( c, f.uvw()[0], l, All );
      MVelocity.apply( c, f.uvw()[1], l, All );
      if ( storage->hasGlobalCells() )
      {
         MVelocity.apply( c, f.uvw()[2], l, All );
      }

      f.uvw().multElementwise( { f.uvw(), outwardNormal.uvw() }, l );
      f.uvw().assign( { rayleighNumber }, { f.uvw() }, l, All );
   }

   calculateStokesResiduals( *A, MVelocity, MPressure, u, f, level, stokesResidual, stokesTmp, residualU );

   initialResiudalU    = residualU;
   vCycleResidualULast = residualU;

   localTimer.start();
   stokesSolver->solve( *A, u, f, level );
   localTimer.end();
   timeStokes = localTimer.last();

   calculateStokesResiduals( *A, MVelocity, MPressure, u, f, level, stokesResidual, stokesTmp, residualU );

   vMax = velocityMaxMagnitude( u.uvw(), uTmp, uMagnitudeSquared, level, All );

   localTimer.start();
   if ( outputInfo.vtk )
   {
      vtkOutput.write( level );
   }
   localTimer.end();
   timeVTK = localTimer.last();

   if ( outputInfo.sphericalTemperatureSlice && timeStep % outputInfo.sphericalTemperatureSliceInterval == 0 )
   {
      if ( verbose )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Evaluating spherical slices ..." );
      }

      evaluateSphericalSliceUV( sliceEvaluationRadius,
                                outputInfo.sphericalTemperatureSliceNumMeridians,
                                outputInfo.sphericalTemperatureSliceNumParallels,
                                c,
                                level,
                                outputInfo.outputDirectory + "/" + outputInfo.outputBaseName + "_temp_slice_uv_" +
                                    std::to_string( timeStep ) + ".csv" );

      evaluateSphericalSliceIco( sliceEvaluationRadius,
                                 outputInfo.sphericalTemperatureSliceIcoRefinements,
                                 c,
                                 level,
                                 outputInfo.outputDirectory + "/" + outputInfo.outputBaseName + "_temp_slice_ico_" +
                                     std::to_string( timeStep ) + ".csv" );
      if ( verbose )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Done." );
      }
   }

   timeStepTimer.end();
   timeStepTotal = timeStepTimer.last();

   WALBERLA_LOG_INFO_ON_ROOT(
       " timestep |           dt |   time total | velocity RMS | velocity max magnitude |   residual u |  total | Stokes |   MMOC |   diff |    VTK |" )
   WALBERLA_LOG_INFO_ON_ROOT(
       "----------+--------------+--------------+--------------+------------------------+--------------+--------+--------+--------+--------+--------+" )
   WALBERLA_LOG_INFO_ON_ROOT(
       walberla::format( " %8s | %12s | %12.8f | %12.4f | %22.4f | %12.5e | %6.2f | %6.2f | %6.2f | %6.2f | %6.2f |",
                         "initial",
                         "-",
                         timeTotal,
                         vRms,
                         vMax,
                         residualU,
                         timeStepTotal,
                         timeStokes,
                         timeMMOC,
                         timeDiffusion,
                         timeVTK ) )

   db.setVariableEntry( "ts", uint_c( timeStep ) );
   db.setVariableEntry( "sim_time", timeTotal );
   db.setVariableEntry( "run_time_advection", timeMMOC );
   db.setVariableEntry( "run_time_diffusion", timeDiffusion );
   db.setVariableEntry( "run_time_stokes", timeStokes );
   db.setVariableEntry( "run_time_time_step", timeStepTotal );
   db.setVariableEntry( "run_time_vtk", timeVTK );
   db.setVariableEntry( "v_max", vMax );
   db.setVariableEntry( "v_rms", vRms );
   db.setVariableEntry( "dt", real_c( 0 ) );

   db.setVariableEntry( "initial_residual_u_predictor", initialResiudalU );
   db.setVariableEntry( "num_v_cycles_predictor", numVCycles );
   db.setVariableEntry( "avg_residual_reduction_u_predictor", averageResidualReductionU / real_c( numVCycles ) );
   db.setVariableEntry( "final_residual_u_predictor", vCycleResidualULast );

   db.setVariableEntry( "initial_residual_u_corrector", real_c( 0 ) );
   db.setVariableEntry( "num_v_cycles_corrector", real_c( 0 ) );
   db.setVariableEntry( "avg_residual_reduction_u_corrector", real_c( 0 ) );
   db.setVariableEntry( "final_residual_u_corrector", real_c( 0 ) );

   db.writeRowOnRoot();

   timer->stop( "Simulation" );

   while ( timeTotal < simulationTime && timeStep < maxNumTimeSteps )
   {
      timer->start( "Simulation" );

      timeStepTimer.start();

      timeStep++;

      // new time step size
      vMax = velocityMaxMagnitude( u.uvw(), uTmp, uTmp2, level, All );

      real_t dt;
      if ( fixedTimeStep )
      {
         dt = dtConstant;
      }
      else
      {
         dt = ( cflMax / vMax ) * hMin;
      }

      // energy

      // advection

      // start value for predictor
      cPr.assign( { 1.0 }, { c }, level, All );

      // let's just use the current velocity for the prediction
      uLast.assign( { 1.0 }, { u }, level, All );

      // diffusion

      cPr.interpolate( initialTemperature, level, DirichletBoundary );

      cOld.assign( { 1.0 }, { cPr }, level, All );

      diffusionOperator.setDt( 0.5 * dt );

      if ( storage->hasGlobalCells() )
      {
         diffusionOperator.getOperator().computeAndStoreLocalElementMatrices();
      }

      calculateDiffusionResidual(
          diffusionSolver, diffusionOperator, L, MVelocity, cPr, cOld, q, cTmp, cTmp2, level, vCycleResidualCLast );

      localTimer.start();
      diffusionSolver.step( diffusionOperator, L, MVelocity, cPr, cOld, q, q, level, Inner | NeumannBoundary | FreeslipBoundary );
      localTimer.end();
      timeDiffusion = localTimer.last();

      localTimer.start();
      transport.step( cPr, u.uvw(), uLast.uvw(), level, All, dt, 1, true );
      localTimer.end();
      timeMMOC = localTimer.last();

      // diffusion

      cPr.interpolate( initialTemperature, level, DirichletBoundary );

      cOld.assign( { 1.0 }, { cPr }, level, All );

      calculateDiffusionResidual(
          diffusionSolver, diffusionOperator, L, MVelocity, cPr, cOld, q, cTmp, cTmp2, level, vCycleResidualCLast );

      localTimer.start();
      diffusionSolver.step( diffusionOperator, L, MVelocity, cPr, cOld, q, q, level, Inner | NeumannBoundary | FreeslipBoundary );
      localTimer.end();
      timeDiffusion += localTimer.last();

      // Stokes

      for ( uint_t l = minLevel; l <= level; l++ )
      {
         MVelocity.apply( cPr, f.uvw()[0], l, All );
         MVelocity.apply( cPr, f.uvw()[1], l, All );
         if ( storage->hasGlobalCells() )
         {
            MVelocity.apply( cPr, f.uvw()[2], l, All );
         }

         f.uvw().multElementwise( { f.uvw(), outwardNormal.uvw() }, l );
         f.uvw().assign( { rayleighNumber }, { f.uvw() }, l, All );
      }

      calculateStokesResiduals( *A, MVelocity, MPressure, u, f, level, stokesResidual, stokesTmp, residualU );

      vCycleResidualULast       = residualU;
      initialResiudalU          = residualU;
      averageResidualReductionU = real_c( 0 );
      numVCycles                = 0;

      localTimer.start();
      stokesSolver->solve( *A, u, f, level );
      localTimer.end();
      timeStokes = localTimer.last();

      calculateStokesResiduals( *A, MVelocity, MPressure, u, f, level, stokesResidual, stokesTmp, residualU );

      db.setVariableEntry( "initial_residual_u_predictor", initialResiudalU );
      db.setVariableEntry( "num_v_cycles_predictor", numVCycles );
      db.setVariableEntry( "avg_residual_reduction_u_predictor", averageResidualReductionU / real_c( numVCycles ) );
      db.setVariableEntry( "final_residual_u_predictor", vCycleResidualULast );

      if ( predictorCorrector )
      {
         // energy

         // diffusion

         c.interpolate( initialTemperature, level, DirichletBoundary );

         cOld.assign( { 1.0 }, { c }, level, All );

         calculateDiffusionResidual(
             diffusionSolver, diffusionOperator, L, MVelocity, c, cOld, q, cTmp, cTmp2, level, vCycleResidualCLast );

         localTimer.start();
         diffusionSolver.step(
             diffusionOperator, L, MVelocity, c, cOld, q, q, level, Inner | NeumannBoundary | FreeslipBoundary );
         localTimer.end();
         timeDiffusion += localTimer.last();

         // advection

         localTimer.start();
         transport.step( c, u.uvw(), uLast.uvw(), level, All, dt, 1, true );
         localTimer.end();
         timeMMOC += localTimer.last();

         // diffusion

         c.interpolate( initialTemperature, level, DirichletBoundary );

         cOld.assign( { 1.0 }, { c }, level, All );

         calculateDiffusionResidual(
             diffusionSolver, diffusionOperator, L, MVelocity, c, cOld, q, cTmp, cTmp2, level, vCycleResidualCLast );

         localTimer.start();
         diffusionSolver.step(
             diffusionOperator, L, MVelocity, c, cOld, q, q, level, Inner | NeumannBoundary | FreeslipBoundary );
         localTimer.end();
         timeDiffusion += localTimer.last();

         // Stokes

         for ( uint_t l = minLevel; l <= level; l++ )
         {
            MVelocity.apply( c, f.uvw()[0], l, All );
            MVelocity.apply( c, f.uvw()[1], l, All );
            if ( storage->hasGlobalCells() )
            {
               MVelocity.apply( c, f.uvw()[2], l, All );
            }

            f.uvw().multElementwise( { f.uvw(), outwardNormal.uvw() }, l );
            f.uvw().assign( { rayleighNumber }, { f.uvw() }, l, All );
         }

         calculateStokesResiduals( *A, MVelocity, MPressure, u, f, level, stokesResidual, stokesTmp, residualU );

         vCycleResidualULast       = residualU;
         initialResiudalU          = residualU;
         averageResidualReductionU = real_c( 0 );
         numVCycles                = 0;

         localTimer.start();
         stokesSolver->solve( *A, u, f, level );
         localTimer.end();
         timeStokes += localTimer.last();

         calculateStokesResiduals( *A, MVelocity, MPressure, u, f, level, stokesResidual, stokesTmp, residualU );

         db.setVariableEntry( "initial_residual_u_corrector", initialResiudalU );
         db.setVariableEntry( "num_v_cycles_corrector", numVCycles );
         db.setVariableEntry( "avg_residual_reduction_u_corrector", averageResidualReductionU / real_c( numVCycles ) );
         db.setVariableEntry( "final_residual_u_corrector", vCycleResidualULast );
      }
      else
      {
         // use predicted value
         c.assign( { 1.0 }, { cPr }, level, All );

         db.setVariableEntry( "initial_residual_u_corrector", real_c( 0 ) );
         db.setVariableEntry( "num_v_cycles_corrector", real_c( 0 ) );
         db.setVariableEntry( "avg_residual_reduction_u_corrector", real_c( 0 ) );
         db.setVariableEntry( "final_residual_u_corrector", real_c( 0 ) );
      }

      timeTotal += dt;

      vRms = velocityRMS( u, stokesTmp, MVelocity, domainInfo.domainVolume(), level );
      vMax = velocityMaxMagnitude( u.uvw(), uTmp, uMagnitudeSquared, level, All );

      localTimer.start();
      if ( outputInfo.vtk )
      {
         vtkOutput.write( level, timeStep );
      }
      localTimer.end();
      timeVTK = localTimer.last();

      if ( outputInfo.sphericalTemperatureSlice && timeStep % outputInfo.sphericalTemperatureSliceInterval == 0 )
      {
         if ( verbose )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "Evaluating spherical slices ..." );
         }

         evaluateSphericalSliceUV( sliceEvaluationRadius,
                                   outputInfo.sphericalTemperatureSliceNumMeridians,
                                   outputInfo.sphericalTemperatureSliceNumParallels,
                                   c,
                                   level,
                                   outputInfo.outputDirectory + "/" + outputInfo.outputBaseName + "_temp_slice_uv_" +
                                       std::to_string( timeStep ) + ".csv" );

         evaluateSphericalSliceIco( sliceEvaluationRadius,
                                    outputInfo.sphericalTemperatureSliceIcoRefinements,
                                    c,
                                    level,
                                    outputInfo.outputDirectory + "/" + outputInfo.outputBaseName + "_temp_slice_ico_" +
                                        std::to_string( timeStep ) + ".csv" );
         if ( verbose )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "Done." );
         }
      }

      timeStepTimer.end();
      timeStepTotal = timeStepTimer.last();

      db.setVariableEntry( "ts", uint_c( timeStep ) );
      db.setVariableEntry( "sim_time", timeTotal );
      db.setVariableEntry( "run_time_advection", timeMMOC );
      db.setVariableEntry( "run_time_diffusion", timeDiffusion );
      db.setVariableEntry( "run_time_stokes", timeStokes );
      db.setVariableEntry( "run_time_time_step", timeStepTotal );
      db.setVariableEntry( "run_time_vtk", timeVTK );
      db.setVariableEntry( "v_max", vMax );
      db.setVariableEntry( "v_rms", vRms );
      db.setVariableEntry( "dt", dt );

      db.writeRowOnRoot();

      if ( outputInfo.printInterval > 0 && timeStep % outputInfo.printInterval == 0 )
      {
         WALBERLA_LOG_INFO_ON_ROOT(
             walberla::format( " %8d | %12.5e | %12.8f | %12.4f | %22.4f | %12.5e | %6.2f | %6.2f | %6.2f | %6.2f | %6.2f |",
                               timeStep,
                               dt,
                               timeTotal,
                               vRms,
                               vMax,
                               residualU,
                               timeStepTotal,
                               timeStokes,
                               timeMMOC,
                               timeDiffusion,
                               timeVTK ) )
      }

      timer->stop( "Simulation" );

      if ( outputTimingJSON )
      {
         if ( verbose )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "Writing timing tree to .json ..." );
         }

         writeTimingTreeJSON( *timer,
                              outputInfo.outputDirectory + "/" + outputInfo.outputBaseName + "_up_to_ts_" +
                                  std::to_string( timeStep ) + "_timing.json" );
      }
   }
}

} // namespace hyteg

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

#ifdef HYTEG_BUILD_WITH_PETSC
   hyteg::PETScManager manager( &argc, &argv );
#endif

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./MantleConvection.prm";
      WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   hyteg::DomainInfo domainInfo;
   hyteg::SolverInfo solverInfo;
   hyteg::OutputInfo outputInfo;

   domainInfo.threeDim = mainConf.getParameter< bool >( "threeDim" );
   domainInfo.rMin     = mainConf.getParameter< hyteg::real_t >( "rMin" );
   domainInfo.rMax     = mainConf.getParameter< hyteg::real_t >( "rMax" );
   domainInfo.nTan     = mainConf.getParameter< hyteg::uint_t >( "nTan" );
   domainInfo.nRad     = mainConf.getParameter< hyteg::uint_t >( "nRad" );

   const hyteg::real_t cflMax          = mainConf.getParameter< hyteg::real_t >( "cflMax" );
   const hyteg::real_t rayleighNumber  = mainConf.getParameter< hyteg::real_t >( "rayleighNumber" );
   const bool          fixedTimeStep   = mainConf.getParameter< bool >( "fixedTimeStep" );
   const hyteg::real_t dtConstant      = mainConf.getParameter< hyteg::real_t >( "dtConstant" );
   const uint_t        minLevel        = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t        level           = mainConf.getParameter< uint_t >( "level" );
   const hyteg::real_t simulationTime  = mainConf.getParameter< hyteg::real_t >( "simulationTime" );
   const hyteg::uint_t maxNumTimeSteps = mainConf.getParameter< hyteg::uint_t >( "maxNumTimeSteps" );

   const int stokesSolverTypeInt    = mainConf.getParameter< int >( "stokesSolverType" );
   const int diffusionSolverTypeInt = mainConf.getParameter< int >( "diffusionSolverType" );

   solverInfo.stokesSolverType                 = static_cast< hyteg::StokesSolverType >( stokesSolverTypeInt );
   solverInfo.stokesMaxNumIterations           = mainConf.getParameter< uint_t >( "stokesMaxNumIterations" );
   solverInfo.stokesAbsoluteResidualUTolerance = mainConf.getParameter< real_t >( "stokesAbsoluteResidualUTolerance" );
   solverInfo.stokesRelativeResidualUTolerance = mainConf.getParameter< real_t >( "stokesRelativeResidualUTolerance" );
   solverInfo.uzawaOmega                       = mainConf.getParameter< real_t >( "uzawaOmega" );
   solverInfo.uzawaInnerIterations             = mainConf.getParameter< uint_t >( "uzawaInnerIterations" );
   solverInfo.uzawaPreSmooth                   = mainConf.getParameter< uint_t >( "uzawaPreSmooth" );
   solverInfo.uzawaPostSmooth                  = mainConf.getParameter< uint_t >( "uzawaPostSmooth" );

   solverInfo.diffusionSolverType                 = static_cast< hyteg::DiffusionSolverType >( diffusionSolverTypeInt );
   solverInfo.diffusionMaxNumIterations           = mainConf.getParameter< uint_t >( "diffusionMaxNumIterations" );
   solverInfo.diffusionAbsoluteResidualUTolerance = mainConf.getParameter< real_t >( "diffusionAbsoluteResidualUTolerance" );

   outputInfo.outputDirectory                         = mainConf.getParameter< std::string >( "outputDirectory" );
   outputInfo.outputBaseName                          = mainConf.getParameter< std::string >( "outputBaseName" );
   outputInfo.vtk                                     = mainConf.getParameter< bool >( "vtk" );
   outputInfo.vtkOutputVelocity                       = mainConf.getParameter< bool >( "vtkOutputVelocity" );
   outputInfo.vtkInterval                             = mainConf.getParameter< uint_t >( "vtkOutputInterval" );
   outputInfo.vtkOutputVertexDoFs                     = mainConf.getParameter< bool >( "vtkOutputVertexDoFs" );
   outputInfo.printInterval                           = 1;
   outputInfo.sphericalTemperatureSlice               = mainConf.getParameter< bool >( "sphericalTemperatureSlice" );
   outputInfo.sphericalTemperatureSliceInterval       = mainConf.getParameter< uint_t >( "sphericalTemperatureSliceInterval" );
   outputInfo.sphericalTemperatureSliceNumMeridians   = mainConf.getParameter< int >( "sphericalTemperatureSliceNumMeridians" );
   outputInfo.sphericalTemperatureSliceNumParallels   = mainConf.getParameter< int >( "sphericalTemperatureSliceNumParallels" );
   outputInfo.sphericalTemperatureSliceIcoRefinements = mainConf.getParameter< int >( "sphericalTemperatureSliceIcoRefinements" );

   const bool verbose = mainConf.getParameter< bool >( "verbose" );

   hyteg::runBenchmark( cflMax,
                        rayleighNumber,
                        fixedTimeStep,
                        dtConstant,
                        maxNumTimeSteps,
                        minLevel,
                        level,
                        solverInfo,
                        domainInfo,
                        true,
                        simulationTime,
                        outputInfo,
                        verbose );
}
