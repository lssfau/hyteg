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
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"
#include "hyteg/solvers/controlflow/SolverLoop.hpp"
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
                               real_t&                     residualV,
                               real_t&                     residualW,
                               real_t&                     residualP )
{
   tmp.interpolate( 0, level, All );
   r.interpolate( 0, level, All );
   A.apply( x, tmp, level, Inner | NeumannBoundary | FreeslipBoundary );
   r.assign( {1.0, -1.0}, {b, tmp}, level, Inner | NeumannBoundary | FreeslipBoundary );
   residualU = normL2( r.uvw.u, tmp.uvw.u, Mu, level, Inner | NeumannBoundary | FreeslipBoundary );
   residualV = normL2( r.uvw.v, tmp.uvw.v, Mu, level, Inner | NeumannBoundary | FreeslipBoundary );
   residualW = real_c( 0 );
   if ( x.getStorage()->hasGlobalCells() )
   {
      residualW = normL2( r.uvw.w, tmp.uvw.w, Mu, level, Inner | NeumannBoundary | FreeslipBoundary );
   }
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

std::shared_ptr< SetupPrimitiveStorage > createSetupStorage( DomainInfo domainInfo, bool freeSlip )
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
   auto innerBoundary = [&]( const Point3D& x ) { return abs( x.norm() - domainInfo.rMin ) <= 1e-8; };
   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   if ( freeSlip )
   {
      setupStorage->setMeshBoundaryFlagsByVertexLocation( 2, innerBoundary );
   }

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

void runBenchmark( real_t      cflMax,
                   real_t      rayleighNumber,
                   bool        fixedTimeStep,
                   real_t      dtConstant,
                   uint_t      level,
                   DomainInfo  domainInfo,
                   bool        freeSlip,
                   bool        predictorCorrector,
                   real_t      simulationTime,
                   bool        vtk,
                   uint_t      printInterval,
                   uint_t      vtkInterval,
                   std::string dbFile )
{
   walberla::WcTimer localTimer;

   const std::string benchmarkName    = "MantleConvection";
   const bool        outputTimingJSON = true;

   auto setupStorage = createSetupStorage( domainInfo, freeSlip );
   auto storage      = std::make_shared< PrimitiveStorage >( *setupStorage );

   if ( vtk )
   {
      writeDomainPartitioningVTK( storage, "vtk/", benchmarkName + "_Domain" );
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

   WALBERLA_LOG_INFO_ON_ROOT( "Benchmark name: " << benchmarkName )
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
   WALBERLA_LOG_INFO_ON_ROOT( "   + free-slip (inner boundary)                    " << ( freeSlip ? "yes" : "no" ) )
   WALBERLA_LOG_INFO_ON_ROOT( " - app settings: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + VTK:                                          " << ( vtk ? "yes" : "no" ) )
   if ( vtk )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "   + VTK interval:                                 " << vtkInterval )
   }
   WALBERLA_LOG_INFO_ON_ROOT( "   + print interval:                               " << printInterval )
   WALBERLA_LOG_INFO_ON_ROOT( "   + database file:                                " << dbFile )
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   FixedSizeSQLDB db( dbFile );

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

   typedef P2P1TaylorHoodFunction< real_t >                                     StokesFunction;
   typedef P2Function< real_t >                                                 ScalarFunction;
   typedef P2P1ElementwiseBlendingStokesOperator                                StokesOperator;
   typedef P2ElementwiseBlendingLaplaceOperator                                 LaplaceOperator;
   typedef P2ElementwiseBlendingMassOperator                                    MassOperatorVelocity;
   typedef P1ElementwiseBlendingMassOperator                                    MassOperatorPressure;
   typedef P2ElementwiseUnsteadyDiffusionOperator                               UnsteadyDiffusionOperator;
   typedef P2ProjectNormalOperator                                              ProjectNormalOperator;
   typedef StrongFreeSlipWrapper< StokesOperator, ProjectNormalOperator, true > StokesOperatorFreeSlip;

   BoundaryCondition bcVelocity;
   BoundaryCondition bcTemperature;

   if ( freeSlip )
   {
      bcVelocity.createDirichletBC( "outerBoundaryVelocity", 1 );
      bcVelocity.createFreeslipBC( "innerBoundaryVelocity", 2 );
      bcTemperature.createDirichletBC( "outerBoundaryTemperature", 1 );
      bcTemperature.createDirichletBC( "innerBoundaryTemperature", 2 );
   }
   else
   {
      bcVelocity.createDirichletBC( "outerBoundaryVelocity", 1 );
      bcVelocity.createDirichletBC( "innerBoundaryVelocity", 1 );
      bcTemperature.createDirichletBC( "outerBoundaryTemperature", 1 );
      bcTemperature.createDirichletBC( "innerBoundaryTemperature", 1 );
   }

   ScalarFunction c( "c", storage, level, level, bcTemperature );
   ScalarFunction cPr( "cPr", storage, level, level, bcTemperature );
   ScalarFunction cOld( "cOld", storage, level, level, bcTemperature );
   ScalarFunction cTmp( "cTmp", storage, level, level, bcTemperature );
   ScalarFunction cTmp2( "cTmp2", storage, level, level, bcTemperature );
   ScalarFunction q( "q", storage, level, level, bcTemperature );

   StokesFunction u( "u", storage, level, level, bcVelocity );
   StokesFunction uLast( "uLast", storage, level, level, bcVelocity );
   StokesFunction f( "f", storage, level, level, bcVelocity );
   StokesFunction outwardNormal( "outwardNormal", storage, level, level, bcVelocity );
   StokesFunction stokesTmp( "stokesTmp", storage, level, level, bcVelocity );
   StokesFunction stokesResidual( "stokesResidual", storage, level, level, bcVelocity );

   ScalarFunction uTmp( "uTmp", storage, level, level, bcVelocity );
   ScalarFunction uTmp2( "uTmp2", storage, level, level, bcVelocity );

   std::function< real_t( const Point3D& ) > initialTemperature = [&]( const Point3D& x ) {
      auto radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
      return std::exp( -initialTemperatureSteepness * ( ( radius - domainInfo.rMin ) / ( domainInfo.rMax - domainInfo.rMin ) ) );
   };

   c.interpolate( initialTemperature, level, All );
   q.interpolate( internalHeating, level, All );

   std::function< real_t( const Point3D& ) > normalX = []( const Point3D& x ) { return x[0] / x.norm(); };
   std::function< real_t( const Point3D& ) > normalY = []( const Point3D& x ) { return x[1] / x.norm(); };
   std::function< real_t( const Point3D& ) > normalZ = []( const Point3D& x ) { return x[2] / x.norm(); };

   std::function< void( const Point3D&, Point3D& ) > normalFreeSlip = []( const Point3D& x, Point3D& n ) { n = x / x.norm(); };

   outwardNormal.uvw.u.interpolate( normalX, level );
   outwardNormal.uvw.v.interpolate( normalY, level );
   if ( storage->hasGlobalCells() )
   {
      outwardNormal.uvw.w.interpolate( normalY, level );
   }

   UnsteadyDiffusionOperator diffusionOperator( storage, level, level, 1.0, diffusivity, DiffusionTimeIntegrator::ImplicitEuler );
   auto projectNormalOperator = std::make_shared< ProjectNormalOperator >( storage, level, level, normalFreeSlip );
   auto A                     = std::make_shared< StokesOperator >( storage, level, level );
   StokesOperatorFreeSlip          AFS( A, projectNormalOperator, FreeslipBoundary );
   LaplaceOperator                 L( storage, level, level );
   MassOperatorVelocity            MVelocity( storage, level, level );
   MassOperatorPressure            MPressure( storage, level, level );
   MMOCTransport< ScalarFunction > transport( storage, setupStorage, level, level, TimeSteppingScheme::RK4 );



   std::shared_ptr< Solver< StokesOperatorFreeSlip > > stokesSolver;
#ifdef HYTEG_BUILD_WITH_PETSC
   if ( freeSlip )
   {
      auto stokesSolverTmp = std::make_shared< PETScBlockPreconditionedStokesSolver< StokesOperatorFreeSlip > >(
          storage, level, 1e-10, 50000, 4, 1 );
      stokesSolverTmp->reassembleMatrix( false );
      stokesSolver = stokesSolverTmp;
   }
   else
   {
//      stokesSolver = std::make_shared< PETScLUSolver< StokesOperatorFreeSlip > >( storage, level );

//      auto stokesSolverTmp = std::make_shared< PETScBlockPreconditionedStokesSolver< StokesOperatorFreeSlip > >(
//          storage, level, 1e-10, 50000, 1, 1 );
//      stokesSolverTmp->reassembleMatrix( false );
//      stokesSolver = stokesSolverTmp;

      stokesSolver = solvertemplates::stokesMinResSolver< StokesOperatorFreeSlip >( storage, level, 1e-10, 50000, true );
   }

   auto internalDiffusionSolver = std::make_shared< PETScMinResSolver< UnsteadyDiffusionOperator > >( storage, level, 1e-10, 50000 );
   internalDiffusionSolver->reassembleMatrix( true );
#else

   auto internalDiffusionSolver =
       std::make_shared< CGSolver< UnsteadyDiffusionOperator > >( storage, level, level, 50000, 1e-10 );

   auto stokesSolver = solvertemplates::stokesMinResSolver< StokesOperatorFreeSlip >( storage, level, 1e-14, 5000 );
#endif

   UnsteadyDiffusion< ScalarFunction, UnsteadyDiffusionOperator, LaplaceOperator, MassOperatorVelocity > diffusionSolver(
       storage, level, level, internalDiffusionSolver );

   real_t timeTotal = 0;
   real_t vMax      = 0;
   real_t nu        = 0;
   real_t vRms      = 0;
   real_t residualU = 0;
   real_t residualV = 0;
   real_t residualW = 0;
   real_t residualP = 0;

   real_t timeStepTotal = 0;
   real_t timeStokes    = 0;
   real_t timeMMOC      = 0;
   real_t timeDiffusion = 0;

   hyteg::VTKOutput vtkOutput( "./vtk", benchmarkName, storage, vtkInterval );

   vtkOutput.add( u );
   vtkOutput.add( f );
   vtkOutput.add( c );
   vtkOutput.add( outwardNormal );
   vtkOutput.add( q );
   vtkOutput.add( stokesResidual );

   WALBERLA_LOG_INFO_ON_ROOT(
       " timestep |           dt | time total | Nusselt number | velocity RMS | velocity max magnitude |   residual u |   residual v |   residual p |  total | Stokes |   MMOC |   diff " )
   WALBERLA_LOG_INFO_ON_ROOT(
       "----------+--------------+------------+----------------+--------------+------------------------+--------------+--------------+--------------+--------+--------+--------+--------" )
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
       " %8s | %12s | %10.5f | %14.4f | %12.4f | %22.4f | %12.5e | %12.5e | %12.5e | %6.2f | %6.2f | %6.2f | %6.2f ",
       "initial",
       "-",
       timeTotal,
       nu,
       vRms,
       vMax,
       residualU,
       residualV,
       residualP,
       timeStepTotal,
       timeStokes,
       timeMMOC,
       timeDiffusion ) )

   timer->stop( "Setup" );

   timer->start( "Simulation" );

   uint_t timeStep = 0;

   MVelocity.apply( c, f.uvw.u, level, All );
   MVelocity.apply( c, f.uvw.v, level, All );
   if ( storage->hasGlobalCells() )
   {
      MVelocity.apply( c, f.uvw.w, level, All );
   }

   f.uvw.u.multElementwise( {f.uvw.u, outwardNormal.uvw.u}, level );
   f.uvw.v.multElementwise( {f.uvw.v, outwardNormal.uvw.v}, level );
   f.uvw.w.multElementwise( {f.uvw.w, outwardNormal.uvw.w}, level );
   f.uvw.u.assign( {rayleighNumber}, {f.uvw.u}, level, All );
   f.uvw.v.assign( {rayleighNumber}, {f.uvw.v}, level, All );
   f.uvw.w.assign( {rayleighNumber}, {f.uvw.w}, level, All );
   projectNormalOperator->apply( f, level, FreeslipBoundary );
   stokesSolver->solve( AFS, u, f, level );

   if ( storage->hasGlobalCells() )
   {
      vMax = velocityMaxMagnitude( u.uvw.u, u.uvw.v, u.uvw.w, uTmp, uTmp2, level, All );
   }
   else
   {
      vMax = velocityMaxMagnitude( u.uvw.u, u.uvw.v, uTmp, uTmp2, level, All );
   }

   if ( vtk )
      vtkOutput.write( level );

   walberla::WcTimer timeStepTimer;

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

      MVelocity.apply( cPr, f.uvw.u, level, All );
      MVelocity.apply( cPr, f.uvw.v, level, All );
      if ( storage->hasGlobalCells() )
      {
         MVelocity.apply( cPr, f.uvw.w, level, All );
      }
      f.uvw.u.multElementwise( {f.uvw.u, outwardNormal.uvw.u}, level );
      f.uvw.v.multElementwise( {f.uvw.v, outwardNormal.uvw.v}, level );
      f.uvw.w.multElementwise( {f.uvw.w, outwardNormal.uvw.w}, level );
      f.uvw.u.assign( {rayleighNumber}, {f.uvw.u}, level, All );
      f.uvw.v.assign( {rayleighNumber}, {f.uvw.v}, level, All );
      f.uvw.w.assign( {rayleighNumber}, {f.uvw.w}, level, All );
      projectNormalOperator->apply( f, level, FreeslipBoundary );

      localTimer.start();
      stokesSolver->solve( AFS, u, f, level );
      localTimer.end();
      timeStokes = localTimer.last();

      calculateStokesResiduals(
          AFS, MVelocity, MPressure, u, f, level, stokesResidual, stokesTmp, residualU, residualV, residualW, residualP );

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

         MVelocity.apply( c, f.uvw.u, level, All );
         MVelocity.apply( c, f.uvw.v, level, All );
         if ( storage->hasGlobalCells() )
         {
            MVelocity.apply( c, f.uvw.w, level, All );
         }
         f.uvw.u.multElementwise( {f.uvw.u, outwardNormal.uvw.u}, level );
         f.uvw.v.multElementwise( {f.uvw.v, outwardNormal.uvw.v}, level );
         f.uvw.w.multElementwise( {f.uvw.w, outwardNormal.uvw.w}, level );
         f.uvw.u.assign( {rayleighNumber}, {f.uvw.u}, level, All );
         f.uvw.v.assign( {rayleighNumber}, {f.uvw.v}, level, All );
         f.uvw.w.assign( {rayleighNumber}, {f.uvw.w}, level, All );
         projectNormalOperator->apply( f, level, FreeslipBoundary );

         localTimer.start();
         stokesSolver->solve( AFS, u, f, level );
         localTimer.end();
         timeStokes += localTimer.last();

         calculateStokesResiduals(
             AFS, MVelocity, MPressure, u, f, level, stokesResidual, stokesTmp, residualU, residualV, residualW, residualP );
      }
      else
      {
         // use predicted value
         c.assign( {1.0}, {cPr}, level, All );
      }

      timeTotal += dt;

      vRms = velocityRMS( u, stokesTmp, MVelocity, domainInfo.domainArea(), level );

      if ( vtk )
      {
         vtkOutput.write( level, timeStep );
      }

      db.setVariableEntry( "ts", uint_c( timeStep ) );
      db.setVariableEntry( "sim_time", timeTotal );
      db.setVariableEntry( "run_time_advection", timeMMOC );
      db.setVariableEntry( "run_time_stokes", timeStokes );
      db.setVariableEntry( "run_time_time_step", timeStepTotal );
      db.setVariableEntry( "v_max", vMax );
      db.setVariableEntry( "v_rms", vRms );
      db.setVariableEntry( "nu", nu );
      db.setVariableEntry( "dt", dt );

      db.writeRowOnRoot();

      timeStepTimer.end();
      timeStepTotal = timeStepTimer.last();

      if ( printInterval > 0 && timeStep % printInterval == 0 )
      {
         WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
             " %8d | %12.5e | %10.5f | %14.4f | %12.4f | %22.4f | %12.5e | %12.5e | %12.5e | %6.2f | %6.2f | %6.2f | %6.2f ",
             timeStep,
             dt,
             timeTotal,
             nu,
             vRms,
             vMax,
             residualU,
             residualV,
             residualP,
             timeStepTotal,
             timeStokes,
             timeMMOC,
             timeDiffusion ) )
      }
   }

   timer->stop( "Simulation" );

   timer->stop( "Total" );

   if ( outputTimingJSON )
   {
      writeTimingTreeJSON( *timer, benchmarkName + "Timing.json" );
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

   const bool freeSlip = mainConf.getParameter< bool >( "freeSlip" );

   const hyteg::real_t cflMax         = mainConf.getParameter< hyteg::real_t >( "cflMax" );
   const hyteg::real_t rayleighNumber = mainConf.getParameter< hyteg::real_t >( "rayleighNumber" );
   const bool          fixedTimeStep  = mainConf.getParameter< bool >( "fixedTimeStep" );
   const hyteg::real_t dtConstant     = mainConf.getParameter< hyteg::real_t >( "dtConstant" );
   const uint_t        level          = mainConf.getParameter< uint_t >( "level" );
   const hyteg::real_t simulationTime = mainConf.getParameter< hyteg::real_t >( "simulationTime" );

   const std::string dbFile = mainConf.getParameter< std::string >( "dbFile" );
   const bool        vtk    = mainConf.getParameter< bool >( "vtk" );

   hyteg::runBenchmark(
       cflMax, rayleighNumber, fixedTimeStep, dtConstant, level, domainInfo, freeSlip, false, simulationTime, vtk, 1, 1, dbFile );
}
