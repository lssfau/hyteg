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

#include <cmath>
#include <core/Environment.h>
#include <core/math/Constants.h>

#include "core/DataTypes.h"
#include "core/config/Config.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/FunctionProperties.hpp"
#include "hyteg/MeshQuality.hpp"
#include "hyteg/composites/StrongFreeSlipWrapper.hpp"
#include "hyteg/composites/UnsteadyDiffusion.hpp"
#include "hyteg/dataexport/SQL.hpp"
#include "hyteg/dataexport/TimingOutput.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/numerictools/CFDHelpers.hpp"
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

struct DomainInfo
{
   bool threeDim = false;

   real_t rMin = 0;
   real_t rMax = 0;
   uint_t nTan = 0;
   uint_t nRad = 0;

   real_t domainArea()
   {
      if ( !threeDim )
      {
         return walberla::math::pi * rMax * rMax - walberla::math::pi * rMin * rMin;
      }
      else
      {
         WALBERLA_ABORT( "3D not implemented." );
      }
   }
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
                               real_t&                     residualU,
                               real_t&                     residualP )
{
   tmp.interpolate( 0, level, All );
   r.interpolate( 0, level, All );
   A.apply( x, tmp, level, Inner | NeumannBoundary | FreeslipBoundary );
   r.assign( {1.0, -1.0}, {b, tmp}, level, Inner | NeumannBoundary | FreeslipBoundary );
   residualU = normL2Squared( r.uvw.u, tmp.uvw.u, Mu, level, Inner | NeumannBoundary | FreeslipBoundary );
   residualU += normL2Squared( r.uvw.v, tmp.uvw.v, Mu, level, Inner | NeumannBoundary | FreeslipBoundary );
   if ( x.getStorage()->hasGlobalCells() )
   {
      residualU += normL2Squared( r.uvw.w, tmp.uvw.w, Mu, level, Inner | NeumannBoundary | FreeslipBoundary );
   }
   residualU = std::sqrt( residualU );
   residualP = normL2( r.p, tmp.p, Mp, level, Inner | NeumannBoundary | FreeslipBoundary );
}

template < typename StokesFunction, typename VelocityMass >
real_t velocityRMS( const StokesFunction& u, const StokesFunction& tmp, const VelocityMass& M, real_t domainVolume, uint_t level )
{
   auto norm = std::pow( normL2( u.uvw.u, tmp.uvw.u, M, level, All ), 2.0 ) +
               std::pow( normL2( u.uvw.v, tmp.uvw.v, M, level, All ), 2.0 );

   if ( u.getStorage()->hasGlobalCells() )
   {
      norm += std::pow( normL2( u.uvw.w, tmp.uvw.w, M, level, All ), 2.0 );
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
   HYTEG_CG     = 1
};

void runBenchmark( real_t              cflMax,
                   real_t              rayleighNumber,
                   bool                fixedTimeStep,
                   real_t              dtConstant,
                   uint_t              minLevel,
                   uint_t              level,
                   StokesSolverType    stokesSolverType,
                   DiffusionSolverType diffusionSolverType,
                   DomainInfo          domainInfo,
                   bool                predictorCorrector,
                   real_t              simulationTime,
                   bool                vtk,
                   uint_t              printInterval,
                   uint_t              vtkInterval,
                   std::string         outputBaseName,
                   bool                verbose )
{
   walberla::WcTimer localTimer;

   const bool outputTimingJSON = true;

   auto setupStorage = createSetupStorage( domainInfo );
   auto storage      = std::make_shared< PrimitiveStorage >( *setupStorage );

   if ( vtk )
   {
      writeDomainPartitioningVTK( storage, "vtk/", outputBaseName + "_domain" );
   }

   auto timer = storage->getTimingTree();
   timer->start( "Total" );
   timer->start( "Setup" );

   const uint_t unknowns = numberOfGlobalDoFs< P2FunctionTag >( *storage, level );
   const real_t hMin     = MeshQuality::getMinimalEdgeLength( storage, level );
   const real_t hMax     = MeshQuality::getMaximalEdgeLength( storage, level );

   const real_t diffusivity                 = 1.0;
   const real_t internalHeating             = 0.0;
   const real_t initialTemperatureSteepness = 10.0;

   const uint_t stokesMaxIter                   = 100000;
   const real_t stokesAbsoluteResidualTolerance = 1e-10;

   WALBERLA_LOG_INFO_ON_ROOT( "Benchmark name: " << outputBaseName )
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
   WALBERLA_LOG_INFO_ON_ROOT( "   + dimensions:                                   " << ( storage->hasGlobalCells() ? "3" : "2" ) )
   WALBERLA_LOG_INFO_ON_ROOT( "   + level:                                        " << level )
   WALBERLA_LOG_INFO_ON_ROOT( "   + unknowns (== particles), including boundary:  " << unknowns )
   WALBERLA_LOG_INFO_ON_ROOT( "   + h_min:                                        " << hMin )
   WALBERLA_LOG_INFO_ON_ROOT( "   + h_max:                                        " << hMax )
   WALBERLA_LOG_INFO_ON_ROOT( " - benchmark settings: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + Rayleigh number:                              " << rayleighNumber )
   WALBERLA_LOG_INFO_ON_ROOT( "   + internal heating:                             " << internalHeating )
   WALBERLA_LOG_INFO_ON_ROOT( " - app settings: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + Stokes solver:                                " << int_c( stokesSolverType ) )
   WALBERLA_LOG_INFO_ON_ROOT( "   + diffusion solver:                             " << int_c( diffusionSolverType ) )
   if ( vtk )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "   + VTK interval:                                 " << vtkInterval )
   }
   WALBERLA_LOG_INFO_ON_ROOT( "   + print interval:                               " << printInterval )
   WALBERLA_LOG_INFO_ON_ROOT( "   + output base name:                             " << outputBaseName )
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   FixedSizeSQLDB db( outputBaseName + ".db" );

   db.setConstantEntry( "fixed_time_step", fixedTimeStep );
   db.setConstantEntry( "cfl_max", cflMax );
   db.setConstantEntry( "dt_constant", dtConstant );
   db.setConstantEntry( "simulation_time", simulationTime );
   db.setConstantEntry( "level", uint_c( level ) );
   db.setConstantEntry( "unknowns", uint_c( unknowns ) );
   db.setConstantEntry( "h_min", hMin );
   db.setConstantEntry( "h_max", hMax );
   db.setConstantEntry( "num_macro_cells", uint_c( setupStorage->getNumberOfCells() ) );
   db.setConstantEntry( "num_macro_faces", uint_c( setupStorage->getNumberOfFaces() ) );
   db.setConstantEntry( "num_macro_edges", uint_c( setupStorage->getNumberOfEdges() ) );
   db.setConstantEntry( "num_macro_vertices", uint_c( setupStorage->getNumberOfVertices() ) );
   db.setConstantEntry( "num_macro_primitives", uint_c( setupStorage->getNumberOfPrimitives() ) );
   db.setConstantEntry( "diffusivity", diffusivity );

   typedef P2P1TaylorHoodFunction< real_t >       StokesFunction;
   typedef P2Function< real_t >                   ScalarFunction;
   typedef P2P1ElementwiseBlendingStokesOperator  StokesOperator;
   typedef P2ElementwiseBlendingLaplaceOperator   LaplaceOperator;
   typedef P2ElementwiseBlendingMassOperator      MassOperatorVelocity;
   typedef P1ElementwiseBlendingMassOperator      MassOperatorPressure;
   typedef P2ElementwiseUnsteadyDiffusionOperator UnsteadyDiffusionOperator;

   BoundaryCondition bcVelocity;
   BoundaryCondition bcTemperature;

   bcVelocity.createDirichletBC( "outerBoundaryVelocity", 1 );
   bcVelocity.createDirichletBC( "innerBoundaryVelocity", 1 );
   bcTemperature.createDirichletBC( "outerBoundaryTemperature", 1 );
   bcTemperature.createDirichletBC( "innerBoundaryTemperature", 1 );

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

   ScalarFunction uTmp( "uTmp", storage, minLevel, level, bcVelocity );
   ScalarFunction uTmp2( "uTmp2", storage, minLevel, level, bcVelocity );

   std::function< real_t( const Point3D& ) > initialTemperature = [&]( const Point3D& x ) {
      auto radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
      return std::exp( -initialTemperatureSteepness * ( ( radius - domainInfo.rMin ) / ( domainInfo.rMax - domainInfo.rMin ) ) );
   };

   for ( uint_t l = 0; l <= level; l++ )
   {
      c.interpolate( initialTemperature, l, All );
      q.interpolate( internalHeating, l, All );
   }

   std::function< real_t( const Point3D& ) > normalX = []( const Point3D& x ) { return x[0] / x.norm(); };
   std::function< real_t( const Point3D& ) > normalY = []( const Point3D& x ) { return x[1] / x.norm(); };
   std::function< real_t( const Point3D& ) > normalZ = []( const Point3D& x ) { return x[2] / x.norm(); };

   std::function< void( const Point3D&, Point3D& ) > normalFreeSlip = []( const Point3D& x, Point3D& n ) { n = x / x.norm(); };

   for ( uint_t l = 0; l <= level; l++ )
   {
      outwardNormal.uvw.u.interpolate( normalX, l );
      outwardNormal.uvw.v.interpolate( normalY, l );
      if ( storage->hasGlobalCells() )
      {
         outwardNormal.uvw.w.interpolate( normalZ, l );
      }
   }

   UnsteadyDiffusionOperator diffusionOperator(
       storage, minLevel, level, 1.0, diffusivity, DiffusionTimeIntegrator::ImplicitEuler );
   auto                            A = std::make_shared< StokesOperator >( storage, minLevel, level );
   LaplaceOperator                 L( storage, minLevel, level );
   MassOperatorVelocity            MVelocity( storage, minLevel, level );
   MassOperatorPressure            MPressure( storage, minLevel, level );
   MMOCTransport< ScalarFunction > transport( storage, setupStorage, minLevel, level, TimeSteppingScheme::RK4 );

   std::shared_ptr< Solver< StokesOperator > > stokesSolver;
   std::shared_ptr< Solver< StokesOperator > > stokesSolverBlockPrecMinRes;

   real_t vCycleResidualULast = 0;

   if ( stokesSolverType == StokesSolverType::PETSC_MUMPS )
   {
      stokesSolver = std::make_shared< PETScLUSolver< StokesOperator > >( storage, level );
   }

   else if ( stokesSolverType == StokesSolverType::PETSC_MINRES_JACOBI ||
             stokesSolverType == StokesSolverType::PETSC_MINRES_BOOMER )
   {
      auto velocityPreconditionerType = 1;
      if ( stokesSolverType == StokesSolverType::PETSC_MINRES_BOOMER )
      {
         velocityPreconditionerType = 3;
      }
      auto stokesSolverTmp = std::make_shared< PETScBlockPreconditionedStokesSolver< StokesOperator > >(
          storage, level, stokesAbsoluteResidualTolerance, stokesMaxIter, velocityPreconditionerType, 1 );
      stokesSolverTmp->reassembleMatrix( false );
      stokesSolverTmp->setVerbose( verbose );
      stokesSolver = stokesSolverTmp;
   }

   else if ( stokesSolverType == StokesSolverType::HYTEG_MINRES )
   {
      stokesSolver = solvertemplates::stokesMinResSolver< StokesOperator >(
          storage, level, stokesAbsoluteResidualTolerance, stokesMaxIter, verbose );
   }

   else if ( stokesSolverType == StokesSolverType::HYTEG_MINRES_GMG )
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

      auto stokesSolverTmp =
          std::make_shared< MinResSolver< StokesOperator > >( storage, minLevel, level, stokesMaxIter, 1e-30, preconditioner );
      stokesSolverTmp->setPrintInfo( verbose );
      stokesSolverTmp->setAbsoluteTolerance( stokesAbsoluteResidualTolerance );
      if ( stokesSolverType == StokesSolverType::HYTEG_MINRES_GMG )
      {
         stokesSolver = stokesSolverTmp;
      }
      else
      {
         stokesSolverBlockPrecMinRes = stokesSolverTmp;
      }
   }

   else if ( stokesSolverType == StokesSolverType::HYTEG_UZAWA_V || stokesSolverType == StokesSolverType::HYTEG_UZAWA_FMG )
   {
      auto prolongationOperator = std::make_shared< P2P1StokesToP2P1StokesProlongation >();
      auto restrictionOperator  = std::make_shared< P2P1StokesToP2P1StokesRestriction >( true );

      std::shared_ptr< Solver< typename StokesOperator::VelocityOperator_T > > smoother =
          std::make_shared< WeightedJacobiSmoother< typename StokesOperator::VelocityOperator_T > >(
              storage, minLevel, level, 0.66 );

      auto uzawaVelocityPreconditioner =
          std::make_shared< StokesVelocityBlockBlockDiagonalPreconditioner< StokesOperator > >( storage, smoother );

      auto uzawaSmoother = std::make_shared< UzawaSmoother< StokesOperator > >(
          storage, uzawaVelocityPreconditioner, minLevel, level, 0.5, Inner | NeumannBoundary, 10 );

      std::shared_ptr< Solver< StokesOperator > > coarseGridSolverInternal;

      auto petscSolverInternalTmp = std::make_shared< PETScBlockPreconditionedStokesSolver< StokesOperator > >(
          storage, minLevel, stokesAbsoluteResidualTolerance, stokesMaxIter, 1 );
      petscSolverInternalTmp->setVerbose( verbose );
      auto coarseGridSolver = petscSolverInternalTmp;

      auto multigridSolver = std::make_shared< GeometricMultigridSolver< StokesOperator > >( storage,
                                                                                             uzawaSmoother,
                                                                                             coarseGridSolver,
                                                                                             restrictionOperator,
                                                                                             prolongationOperator,
                                                                                             minLevel,
                                                                                             level,
                                                                                             10,
                                                                                             10,
                                                                                             2,
                                                                                             CycleType::VCYCLE );

      if ( stokesSolverType == StokesSolverType::HYTEG_UZAWA_V )
      {
         auto stopIterationCallback = [&]() {
            real_t r_u;
            real_t r_p;

            calculateStokesResiduals( *A, MVelocity, MPressure, u, f, level, stokesResidual, stokesTmp, r_u, r_p );

            auto reductionRateU = r_u / vCycleResidualULast;

            vCycleResidualULast = r_u;

            if ( verbose )
            {
               WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
                   "[Uzawa] residual u: %10.5e | reduction: %10.5e | residual p: %10.5e", r_u, reductionRateU, r_p ) );
            }

            if ( r_u < stokesAbsoluteResidualTolerance )
            {
               return true;
            }

            if ( reductionRateU > 0.8 )
            {
               return true;
            }

            return false;
         };

         stokesSolver = std::make_shared< SolverLoop< StokesOperator > >( multigridSolver, 10, stopIterationCallback );
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

   if ( diffusionSolverType == DiffusionSolverType::PETSC_MINRES )
   {
      auto internalDiffusionSolver =
          std::make_shared< PETScMinResSolver< UnsteadyDiffusionOperator > >( storage, level, 1e-10, 50000 );
      internalDiffusionSolver->reassembleMatrix( true );
      diffusionLinearSolver = internalDiffusionSolver;
   }
   else if ( diffusionSolverType == DiffusionSolverType::HYTEG_CG )
   {
      auto internalDiffusionSolver = std::make_shared< CGSolver< UnsteadyDiffusionOperator > >( storage, minLevel, level );
      internalDiffusionSolver->setPrintInfo( verbose );
      diffusionLinearSolver = internalDiffusionSolver;
   }

   UnsteadyDiffusion< ScalarFunction, UnsteadyDiffusionOperator, LaplaceOperator, MassOperatorVelocity > diffusionSolver(
       storage, minLevel, level, diffusionLinearSolver );

   real_t timeTotal = 0;
   real_t vMax      = 0;
   real_t vRms      = 0;
   real_t residualU = 0;
   real_t residualP = 0;

   real_t timeStepTotal = 0;
   real_t timeStokes    = 0;
   real_t timeMMOC      = 0;
   real_t timeDiffusion = 0;
   real_t timeVTK       = 0;

   hyteg::VTKOutput vtkOutput( "./vtk", outputBaseName, storage, vtkInterval );

   vtkOutput.add( u );
   vtkOutput.add( f );
   vtkOutput.add( c );
   vtkOutput.add( outwardNormal );
   vtkOutput.add( q );
   vtkOutput.add( stokesResidual );

   timer->stop( "Setup" );

   timer->start( "Simulation" );

   uint_t timeStep = 0;

   for ( uint_t l = 0; l <= level; l++ )
   {
      MVelocity.apply( c, f.uvw.u, l, All );
      MVelocity.apply( c, f.uvw.v, l, All );
      if ( storage->hasGlobalCells() )
      {
         MVelocity.apply( c, f.uvw.w, l, All );
      }

      f.uvw.u.multElementwise( {f.uvw.u, outwardNormal.uvw.u}, l );
      f.uvw.v.multElementwise( {f.uvw.v, outwardNormal.uvw.v}, l );
      f.uvw.w.multElementwise( {f.uvw.w, outwardNormal.uvw.w}, l );
      f.uvw.u.assign( {rayleighNumber}, {f.uvw.u}, l, All );
      f.uvw.v.assign( {rayleighNumber}, {f.uvw.v}, l, All );
      f.uvw.w.assign( {rayleighNumber}, {f.uvw.w}, l, All );
   }

   calculateStokesResiduals( *A, MVelocity, MPressure, u, f, level, stokesResidual, stokesTmp, residualU, residualP );

   vCycleResidualULast = residualU;

   localTimer.start();
   stokesSolver->solve( *A, u, f, level );
   localTimer.end();
   timeStokes = localTimer.last();

   calculateStokesResiduals( *A, MVelocity, MPressure, u, f, level, stokesResidual, stokesTmp, residualU, residualP );

   if ( storage->hasGlobalCells() )
   {
      vMax = velocityMaxMagnitude( u.uvw.u, u.uvw.v, u.uvw.w, uTmp, uTmp2, level, All );
   }
   else
   {
      vMax = velocityMaxMagnitude( u.uvw.u, u.uvw.v, uTmp, uTmp2, level, All );
   }

   localTimer.start();
   if ( vtk )
   {
      vtkOutput.write( level );
   }
   localTimer.end();
   timeVTK = localTimer.last();

   walberla::WcTimer timeStepTimer;

   WALBERLA_LOG_INFO_ON_ROOT(
       " timestep |           dt |   time total | velocity RMS | velocity max magnitude |   residual u |   residual p |  total | Stokes |   MMOC |   diff |    VTK |" )
   WALBERLA_LOG_INFO_ON_ROOT(
       "----------+--------------+--------------+--------------+------------------------+--------------+--------------+--------+--------+--------+--------+--------+" )
   WALBERLA_LOG_INFO_ON_ROOT(
       walberla::format( " %8s | %12s | %12.8f | %12.4f | %22.4f | %12.5e | %12.5e | %6.2f | %6.2f | %6.2f | %6.2f | %6.2f |",
                         "initial",
                         "-",
                         timeTotal,
                         vRms,
                         vMax,
                         residualU,
                         residualP,
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

   db.writeRowOnRoot();

   while ( timeTotal < simulationTime )
   {
      timeStepTimer.start();

      timeStep++;

      // new time step size

      if ( storage->hasGlobalCells() )
      {
         vMax = velocityMaxMagnitude( u.uvw.u, u.uvw.v, u.uvw.w, uTmp, uTmp2, level, All );
      }
      else
      {
         vMax = velocityMaxMagnitude( u.uvw.u, u.uvw.v, uTmp, uTmp2, level, All );
      }

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
      cPr.assign( {1.0}, {c}, level, All );

      // let's just use the current velocity for the prediction
      uLast.assign( {1.0}, {u}, level, All );

      // diffusion

      cPr.interpolate( initialTemperature, level, DirichletBoundary );

      cOld.assign( {1.0}, {cPr}, level, All );

      diffusionOperator.setDt( 0.5 * dt );

      localTimer.start();
      diffusionSolver.step( diffusionOperator, L, MVelocity, cPr, cOld, q, q, level, Inner | NeumannBoundary | FreeslipBoundary );
      localTimer.end();
      timeDiffusion = localTimer.last();

      localTimer.start();
      transport.step( cPr, u.uvw.u, u.uvw.v, u.uvw.w, uLast.uvw.u, uLast.uvw.v, uLast.uvw.w, level, All, dt, 1, true );
      localTimer.end();
      timeMMOC = localTimer.last();

      // diffusion

      cPr.interpolate( initialTemperature, level, DirichletBoundary );

      cOld.assign( {1.0}, {cPr}, level, All );

      localTimer.start();
      diffusionSolver.step( diffusionOperator, L, MVelocity, cPr, cOld, q, q, level, Inner | NeumannBoundary | FreeslipBoundary );
      localTimer.end();
      timeDiffusion += localTimer.last();

      // Stokes

      for ( uint_t l = 0; l <= level; l++ )
      {
         MVelocity.apply( cPr, f.uvw.u, l, All );
         MVelocity.apply( cPr, f.uvw.v, l, All );
         if ( storage->hasGlobalCells() )
         {
            MVelocity.apply( cPr, f.uvw.w, l, All );
         }

         f.uvw.u.multElementwise( {f.uvw.u, outwardNormal.uvw.u}, l );
         f.uvw.v.multElementwise( {f.uvw.v, outwardNormal.uvw.v}, l );
         f.uvw.w.multElementwise( {f.uvw.w, outwardNormal.uvw.w}, l );
         f.uvw.u.assign( {rayleighNumber}, {f.uvw.u}, l, All );
         f.uvw.v.assign( {rayleighNumber}, {f.uvw.v}, l, All );
         f.uvw.w.assign( {rayleighNumber}, {f.uvw.w}, l, All );
      }

      calculateStokesResiduals( *A, MVelocity, MPressure, u, f, level, stokesResidual, stokesTmp, residualU, residualP );

      vCycleResidualULast = residualU;

      localTimer.start();
      stokesSolver->solve( *A, u, f, level );
      localTimer.end();
      timeStokes = localTimer.last();

      calculateStokesResiduals( *A, MVelocity, MPressure, u, f, level, stokesResidual, stokesTmp, residualU, residualP );

      if ( predictorCorrector )
      {
         // energy

         // diffusion

         c.interpolate( initialTemperature, level, DirichletBoundary );

         cOld.assign( {1.0}, {c}, level, All );

         localTimer.start();
         diffusionSolver.step(
             diffusionOperator, L, MVelocity, c, cOld, q, q, level, Inner | NeumannBoundary | FreeslipBoundary );
         localTimer.end();
         timeDiffusion += localTimer.last();

         // advection

         localTimer.start();
         transport.step( c, u.uvw.u, u.uvw.v, u.uvw.w, uLast.uvw.u, uLast.uvw.v, uLast.uvw.w, level, All, dt, 1, true );
         localTimer.end();
         timeMMOC += localTimer.last();

         // diffusion

         c.interpolate( initialTemperature, level, DirichletBoundary );

         cOld.assign( {1.0}, {c}, level, All );

         localTimer.start();
         diffusionSolver.step(
             diffusionOperator, L, MVelocity, c, cOld, q, q, level, Inner | NeumannBoundary | FreeslipBoundary );
         localTimer.end();
         timeDiffusion += localTimer.last();

         // Stokes

         for ( uint_t l = 0; l <= level; l++ )
         {
            MVelocity.apply( c, f.uvw.u, l, All );
            MVelocity.apply( c, f.uvw.v, l, All );
            if ( storage->hasGlobalCells() )
            {
               MVelocity.apply( c, f.uvw.w, l, All );
            }

            f.uvw.u.multElementwise( {f.uvw.u, outwardNormal.uvw.u}, l );
            f.uvw.v.multElementwise( {f.uvw.v, outwardNormal.uvw.v}, l );
            f.uvw.w.multElementwise( {f.uvw.w, outwardNormal.uvw.w}, l );
            f.uvw.u.assign( {rayleighNumber}, {f.uvw.u}, l, All );
            f.uvw.v.assign( {rayleighNumber}, {f.uvw.v}, l, All );
            f.uvw.w.assign( {rayleighNumber}, {f.uvw.w}, l, All );
         }

         calculateStokesResiduals( *A, MVelocity, MPressure, u, f, level, stokesResidual, stokesTmp, residualU, residualP );

         vCycleResidualULast = residualU;

         localTimer.start();
         stokesSolver->solve( *A, u, f, level );
         localTimer.end();
         timeStokes += localTimer.last();

         calculateStokesResiduals( *A, MVelocity, MPressure, u, f, level, stokesResidual, stokesTmp, residualU, residualP );
      }
      else
      {
         // use predicted value
         c.assign( {1.0}, {cPr}, level, All );
      }

      timeTotal += dt;

      vRms = velocityRMS( u, stokesTmp, MVelocity, domainInfo.domainArea(), level );

      localTimer.start();
      if ( vtk )
      {
         vtkOutput.write( level, timeStep );
      }
      localTimer.end();
      timeVTK = localTimer.last();

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

      timeStepTimer.end();
      timeStepTotal = timeStepTimer.last();

      if ( printInterval > 0 && timeStep % printInterval == 0 )
      {
         WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
             " %8d | %12.5e | %12.8f | %12.4f | %22.4f | %12.5e | %12.5e | %6.2f | %6.2f | %6.2f | %6.2f | %6.2f |",
             timeStep,
             dt,
             timeTotal,
             vRms,
             vMax,
             residualU,
             residualP,
             timeStepTotal,
             timeStokes,
             timeMMOC,
             timeDiffusion,
             timeVTK ) )
      }
   }

   timer->stop( "Simulation" );

   timer->stop( "Total" );

   if ( outputTimingJSON )
   {
      writeTimingTreeJSON( *timer, outputBaseName + "_timing.json" );
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

   domainInfo.threeDim = mainConf.getParameter< bool >( "threeDim" );
   domainInfo.rMin     = mainConf.getParameter< hyteg::real_t >( "rMin" );
   domainInfo.rMax     = mainConf.getParameter< hyteg::real_t >( "rMax" );
   domainInfo.nTan     = mainConf.getParameter< hyteg::uint_t >( "nTan" );
   domainInfo.nRad     = mainConf.getParameter< hyteg::uint_t >( "nRad" );

   const hyteg::real_t cflMax         = mainConf.getParameter< hyteg::real_t >( "cflMax" );
   const hyteg::real_t rayleighNumber = mainConf.getParameter< hyteg::real_t >( "rayleighNumber" );
   const bool          fixedTimeStep  = mainConf.getParameter< bool >( "fixedTimeStep" );
   const hyteg::real_t dtConstant     = mainConf.getParameter< hyteg::real_t >( "dtConstant" );
   const uint_t        minLevel       = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t        level          = mainConf.getParameter< uint_t >( "level" );
   const hyteg::real_t simulationTime = mainConf.getParameter< hyteg::real_t >( "simulationTime" );

   const int stokesSolverTypeInt    = mainConf.getParameter< int >( "stokesSolverType" );
   const int diffusionSolverTypeInt = mainConf.getParameter< int >( "diffusionSolverType" );

   auto stokesSolverType    = static_cast< hyteg::StokesSolverType >( stokesSolverTypeInt );
   auto diffusionSolverType = static_cast< hyteg::DiffusionSolverType >( diffusionSolverTypeInt );

   const std::string outputBaseName = mainConf.getParameter< std::string >( "outputBaseName" );
   const bool        vtk            = mainConf.getParameter< bool >( "vtk" );

   const bool verbose = mainConf.getParameter< bool >( "verbose" );

   hyteg::runBenchmark( cflMax,
                        rayleighNumber,
                        fixedTimeStep,
                        dtConstant,
                        minLevel,
                        level,
                        stokesSolverType,
                        diffusionSolverType,
                        domainInfo,
                        true,
                        simulationTime,
                        vtk,
                        1,
                        1,
                        outputBaseName,
                        verbose );
}
