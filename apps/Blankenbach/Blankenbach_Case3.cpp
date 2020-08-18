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
#include "hyteg/composites/MMOCTransport.hpp"
#include "hyteg/composites/UnsteadyDiffusion.hpp"
#include "hyteg/dataexport/TimingOutput.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/numerictools/CFDHelpers.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
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

#include "sqlite/SQLite.h"

namespace hyteg {

using walberla::int_c;
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
                               real_t&                     residualP )
{
   r.interpolate( 0, level, All );
   A.apply( x, tmp, level, Inner | NeumannBoundary );
   r.assign( {1.0, -1.0}, {b, tmp}, level, Inner | NeumannBoundary );
   residualU = normL2( r.u, tmp.u, Mu, level, Inner | NeumannBoundary );
   residualV = normL2( r.v, tmp.v, Mu, level, Inner | NeumannBoundary );
   residualP = normL2( r.p, tmp.p, Mp, level, Inner | NeumannBoundary );
}

template < typename StokesFunction >
real_t velocityRMS( const StokesFunction& u, real_t domainHeight, real_t domainWidth, uint_t level )
{
   const auto norm = std::sqrt(
       ( u.u.dotGlobal( u.u, level, All ) + u.v.dotGlobal( u.v, level, All ) ) /
       real_c( 2 * numberOfGlobalDoFs< typename StokesFunction::VelocityFunction_T::Tag >( *u.u.getStorage(), level ) ) );
   const auto area = domainHeight * domainWidth;
   return norm / std::sqrt( area );
}

std::shared_ptr< SetupPrimitiveStorage > createSetupStorage()
{
   auto meshInfo     = MeshInfo::meshRectangle( Point2D( {0, 0} ), Point2D( {1.5, 1} ), MeshInfo::CRISS, 3, 2 );
   auto setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   auto topBoundary = []( const Point3D& x ) { return std::abs( x[1] - 1.0 ) < 1e-3; };

   auto bottomBoundary = []( const Point3D& x ) { return std::abs( x[1] ) < 1e-3; };

   auto leftBoundary = []( const Point3D& x ) { return std::abs( x[0] ) < 1e-3; };

   auto rightBoundary = []( const Point3D& x ) { return std::abs( x[0] - 1.5 ) < 1e-3; };

   setupStorage->setMeshBoundaryFlagsByVertexLocation( 3, leftBoundary );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( 4, rightBoundary );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( 1, topBoundary );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( 2, bottomBoundary );

   return setupStorage;
}

void runBenchmark( real_t      cflMax,
                   uint_t      level,
                   bool        adjustedAdvection,
                   real_t      simulationTime,
                   bool        vtk,
                   uint_t      printInterval,
                   uint_t      vtkInterval,
                   bool        verbose,
                   std::string dbFile )
{
   walberla::WcTimer localTimer;

   const std::string benchmarkName    = "Blankenbach_Case3";
   const bool        outputTimingJSON = true;

   auto setupStorage = createSetupStorage();
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

   const real_t rayleighNumber  = 216000.0;
   const real_t diffusivity     = 1.0;
   const real_t internalHeating = 1.0;

   WALBERLA_LOG_INFO_ON_ROOT( "Benchmark name: " << benchmarkName )
   WALBERLA_LOG_INFO_ON_ROOT( " - time discretization: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + max CFL:                                      " << cflMax )
   WALBERLA_LOG_INFO_ON_ROOT( "   + simulation time:                              " << simulationTime )
   WALBERLA_LOG_INFO_ON_ROOT( " - space discretization: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + dimensions:                                   " << ( storage->hasGlobalCells() ? "3" : "2" ) )
   WALBERLA_LOG_INFO_ON_ROOT( "   + level:                                        " << level )
   WALBERLA_LOG_INFO_ON_ROOT( "   + unknowns (== particles), including boundary:  " << unknowns )
   WALBERLA_LOG_INFO_ON_ROOT( "   + h_min:                                        " << hMin )
   WALBERLA_LOG_INFO_ON_ROOT( "   + h_max:                                        " << hMax )
   WALBERLA_LOG_INFO_ON_ROOT( " - advection-diffusion settings: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + diffusivity:                                  " << diffusivity )
   WALBERLA_LOG_INFO_ON_ROOT( "   + adjusted advection:                           " << ( adjustedAdvection ? "yes" : "no" ) )
   WALBERLA_LOG_INFO_ON_ROOT( " - app settings: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + VTK:                                          " << ( vtk ? "yes" : "no" ) )
   if ( vtk )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "   + VTK interval:                                 " << vtkInterval )
   }
   WALBERLA_LOG_INFO_ON_ROOT( "   + print interval:                               " << printInterval )
   WALBERLA_LOG_INFO_ON_ROOT( "   + database file:                                " << dbFile )
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   walberla::sqlite::SQLiteDB                 db( dbFile );
   std::map< std::string, walberla::int64_t > sqlIntegerProperties;
   std::map< std::string, double >            sqlRealProperties;
   std::map< std::string, std::string >       sqlStringProperties;

   sqlIntegerProperties["ts"]                   = 0;
   sqlRealProperties["cfl_max"]                 = cflMax;
   sqlRealProperties["simulation_time"]         = simulationTime;
   sqlIntegerProperties["level"]                = int_c( level );
   sqlIntegerProperties["unknowns"]             = int_c( unknowns );
   sqlRealProperties["h_min"]                   = hMin;
   sqlRealProperties["h_max"]                   = hMax;
   sqlIntegerProperties["num_macro_cells"]      = int_c( setupStorage->getNumberOfCells() );
   sqlIntegerProperties["num_macro_faces"]      = int_c( setupStorage->getNumberOfFaces() );
   sqlIntegerProperties["num_macro_edges"]      = int_c( setupStorage->getNumberOfEdges() );
   sqlIntegerProperties["num_macro_vertices"]   = int_c( setupStorage->getNumberOfVertices() );
   sqlIntegerProperties["num_macro_primitives"] = int_c( setupStorage->getNumberOfPrimitives() );
   sqlRealProperties["diffusivity"]             = diffusivity;
   sqlIntegerProperties["adjusted_advection"]   = adjustedAdvection;

   typedef P2P1TaylorHoodFunction< real_t >    StokesFunction;
   typedef P2Function< real_t >                ScalarFunction;
   typedef P2P1TaylorHoodStokesOperator        StokesOperator;
   typedef P2ConstantLaplaceOperator           LaplaceOperator;
   typedef P2ConstantMassOperator              MassOperatorVelocity;
   typedef P1ConstantMassOperator              MassOperatorPressure;
   typedef P2ConstantUnsteadyDiffusionOperator UnsteadyDiffusionOperator;

   //   BoundaryCondition bcVelocity = BoundaryCondition::create012BC();
   //   BoundaryCondition bcTemperature = BoundaryCondition::create012BC();

   BoundaryCondition bcVelocity;
   BoundaryCondition bcTemperature;

   bcVelocity.createDirichletBC( "topAndBottomWalls", {1, 2} );
   bcVelocity.createDirichletBC( "sideWalls", {3, 4} ); // TODO: change to free-slip

   bcTemperature.createDirichletBC( "dirichletWall", 1 );
   bcTemperature.createNeumannBC( "neumannWalls", {2, 3, 4} );

   ScalarFunction c( "c", storage, level, level, bcTemperature );
   ScalarFunction cOld( "cOld", storage, level, level, bcTemperature );
   ScalarFunction cTmp( "cTmp", storage, level, level, bcTemperature );
   ScalarFunction cTmp2( "cTmp2", storage, level, level, bcTemperature );
   ScalarFunction q( "q", storage, level, level, bcTemperature );

   StokesFunction u( "u", storage, level, level, bcVelocity );
   StokesFunction uLast( "uLast", storage, level, level, bcVelocity );
   StokesFunction f( "f", storage, level, level, bcVelocity );
   StokesFunction fLast( "fLast", storage, level, level, bcVelocity );
   StokesFunction upwardNormal( "upwardNormal", storage, level, level, bcVelocity );
   StokesFunction stokesTmp( "stokesTmp", storage, level, level, bcVelocity );
   StokesFunction stokesResidual( "stokesResidual", storage, level, level, bcVelocity );

   ScalarFunction uTmp( "uTmp", storage, level, level, bcVelocity );
   ScalarFunction uTmp2( "uTmp2", storage, level, level, bcVelocity );

   auto initialTemperature = []( const Point3D& x ) {
      return 0.5 * ( 1.0 - x[1] * x[1] ) + 0.01 * std::cos( pi * x[0] / 1.5 ) * std::sin( pi * x[1] / 1.0 );
   };

   c.interpolate( initialTemperature, level, Inner | NeumannBoundary );
   q.interpolate( internalHeating, level );
   upwardNormal.v.interpolate( 1, level );

   UnsteadyDiffusionOperator diffusionOperator( storage, level, level, 1.0, diffusivity, DiffusionTimeIntegrator::ImplicitEuler );
   StokesOperator            A( storage, level, level );
   LaplaceOperator           L( storage, level, level );
   MassOperatorVelocity      MVelocity( storage, level, level );
   MassOperatorPressure      MPressure( storage, level, level );
   MMOCTransport< ScalarFunction > transport( storage, setupStorage, level, level, TimeSteppingScheme::RK4 );

#ifndef HYTEG_BUILD_WITH_PETSC
   PETScManager manager;
   auto internalDiffusionSolver = std::make_shared< PETScMinResSolver< UnsteadyDiffusionOperator > >( storage, level, 1e-06 );
#else
   auto internalDiffusionSolver =
       std::make_shared< CGSolver< UnsteadyDiffusionOperator > >( storage, level, level, 100000, 1e-12 );
#endif

   auto stokesSolver = solvertemplates::stokesMinResSolver< StokesOperator >( storage, level, 1e-12, 5000 );

   UnsteadyDiffusion< ScalarFunction, UnsteadyDiffusionOperator, LaplaceOperator, MassOperatorVelocity > diffusionSolver(
       storage, level, level, internalDiffusionSolver );

   real_t timeTotal = 0;
   real_t vMax      = velocityMaxMagnitude( u.u, u.v, uTmp, uTmp2, level, All );
   real_t nu        = 0;
   real_t vRms      = 0;
   real_t residualU = 0;
   real_t residualV = 0;
   real_t residualP = 0;

   hyteg::VTKOutput vtkOutput( "./vtk", benchmarkName, storage, vtkInterval );

   vtkOutput.add( u );
   vtkOutput.add( f );
   vtkOutput.add( c );
   vtkOutput.add( upwardNormal );
   vtkOutput.add( q );

   if ( vtk )
      vtkOutput.write( level );

   WALBERLA_LOG_INFO_ON_ROOT(
       " timestep |           dt | time total | Nusselt number | velocity RMS | velocity max magnitude |   residual u |   residual v |   residual p " )
   WALBERLA_LOG_INFO_ON_ROOT(
       "----------+--------------+------------+----------------+--------------+------------------------+--------------+--------------+-------------- " )
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %8s | %12s | %10.5f | %14.4f | %12.4f | %22.4f | %12.5e | %12.5e | %12.5e ",
                                                "initial",
                                                "-",
                                                timeTotal,
                                                nu,
                                                vRms,
                                                vMax,
                                                residualU,
                                                residualV,
                                                residualP ) )

   WALBERLA_ROOT_SECTION()
   {
      sqlRealProperties["sim_time"] = timeTotal;
      sqlRealProperties["v_max"]    = vMax;

      db.storeRun( sqlIntegerProperties, sqlStringProperties, sqlRealProperties );
      sqlRealProperties.clear();
      sqlIntegerProperties.clear();
      sqlStringProperties.clear();
   }

   timer->stop( "Setup" );

   timer->start( "Simulation" );

   MVelocity.apply( c, f.u, level, All );
   MVelocity.apply( c, f.v, level, All );
   f.u.multElementwise( {f.u, upwardNormal.u}, level );
   f.v.multElementwise( {f.v, upwardNormal.v}, level );
   f.u.assign( {rayleighNumber}, {f.u}, level, All );
   f.v.assign( {rayleighNumber}, {f.v}, level, All );

   stokesSolver->solve( A, u, f, level );

   uint_t timeStep = 0;

   while ( timeTotal < simulationTime )
   {
      timeStep++;

      if ( verbose )
         WALBERLA_LOG_INFO_ON_ROOT( "timestep " << timeStep )

      // Stokes

      uLast.assign( {1.0}, {u}, level, All );

      MVelocity.apply( c, f.u, level, All );
      MVelocity.apply( c, f.v, level, All );
      f.u.multElementwise( {f.u, upwardNormal.u}, level );
      f.v.multElementwise( {f.v, upwardNormal.v}, level );
      f.u.assign( {rayleighNumber}, {f.u}, level, All );
      f.v.assign( {rayleighNumber}, {f.v}, level, All );

      stokesSolver->solve( A, u, f, level );

      calculateStokesResiduals(
          A, MVelocity, MPressure, u, f, level, stokesResidual, stokesTmp, residualU, residualV, residualP );
      vRms = velocityRMS( u, 1.0, 1.5, level );

      // Energy

      vMax = velocityMaxMagnitude( u.u, u.v, uTmp, uTmp2, level, All );

      const auto dt = ( cflMax / vMax ) * hMin;

      if ( verbose )
         WALBERLA_LOG_INFO_ON_ROOT( "performing transport step" )

      real_t advectionTimeStepRunTime;

      if ( adjustedAdvection )
      {
         const real_t adjustedAdvectionPertubation = 0.1 * ( hMin / vMax );
         localTimer.start();
         transport.step(
             c, u.u, u.v, u.w, uLast.u, uLast.v, uLast.w, level, All, dt, 1, MVelocity, 0.0, adjustedAdvectionPertubation );
         localTimer.end();
         advectionTimeStepRunTime = localTimer.last();
      }
      else
      {
         localTimer.start();
         transport.step( c, u.u, u.v, u.w, uLast.u, uLast.v, uLast.w, level, All, dt, 1, true );
         localTimer.end();
         advectionTimeStepRunTime = localTimer.last();
      }

      cOld.assign( {1.0}, {c}, level, All );

      timeTotal += dt;

      diffusionOperator.setDt( dt );

      if ( verbose )
         WALBERLA_LOG_INFO_ON_ROOT( "performing diffusion step" )
      diffusionSolver.step( diffusionOperator, L, MVelocity, c, cOld, q, q, level, Inner | NeumannBoundary );
      // diffusionSolver.step( diffusionOperator, L, M, c, cOld, level, Inner | NeumannBoundary );

      if ( printInterval > 0 && timeStep % printInterval == 0 )
      {
         WALBERLA_LOG_INFO_ON_ROOT(
             walberla::format( " %8d | %12.5e | %10.5f | %14.4f | %12.4f | %22.4f | %12.5e | %12.5e | %12.5e ",
                               timeStep,
                               dt,
                               timeTotal,
                               nu,
                               vRms,
                               vMax,
                               residualU,
                               residualV,
                               residualP ) )
      }

      if ( vtk )
         vtkOutput.write( level, timeStep );

      WALBERLA_ROOT_SECTION()
      {
         sqlIntegerProperties["ts"]              = int_c( timeStep );
         sqlRealProperties["sim_time"]           = timeTotal;
         sqlRealProperties["run_time_advection"] = advectionTimeStepRunTime;
         sqlRealProperties["v_max"]              = vMax;

         db.storeRun( sqlIntegerProperties, sqlStringProperties, sqlRealProperties );
         sqlRealProperties.clear();
         sqlIntegerProperties.clear();
         sqlStringProperties.clear();
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

   hyteg::runBenchmark( 10.0, 4, false, 100000, true, 1, 1, false, "database.db" );
}