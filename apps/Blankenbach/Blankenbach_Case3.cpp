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
#include "hyteg/dataexport/TimingOutput.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
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
   tmp.interpolate( 0, level, All );
   r.interpolate( 0, level, All );
   A.apply( x, tmp, level, Inner | NeumannBoundary | FreeslipBoundary );
   r.assign( {1.0, -1.0}, {b, tmp}, level, Inner | NeumannBoundary | FreeslipBoundary );
   residualU = normL2( r.u, tmp.u, Mu, level, Inner | NeumannBoundary | FreeslipBoundary );
   residualV = normL2( r.v, tmp.v, Mu, level, Inner | NeumannBoundary | FreeslipBoundary );
   residualP = normL2( r.p, tmp.p, Mp, level, Inner | NeumannBoundary | FreeslipBoundary );
}

template < typename StokesFunction >
real_t velocityRMS( const StokesFunction& u, real_t domainHeight, real_t domainWidth, uint_t level )
{
   const auto norm = std::sqrt(
       ( u.u.dotGlobal( u.u, level, All ) + u.v.dotGlobal( u.v, level, All ) ) /
       real_c( numberOfGlobalDoFs< typename StokesFunction::VelocityFunction_T::Tag >( *u.u.getStorage(), level ) ) );

   const auto area = domainHeight * domainWidth;
   return norm / std::sqrt( area );
}

template < typename FunctionType >
std::vector< real_t > evaluateHorizontalTemperatureSlice( const FunctionType& c,
                                                          uint_t              level,
                                                          real_t              xMin,
                                                          real_t              xMax,
                                                          real_t              y,
                                                          uint_t              numSamples )
{
   std::vector< real_t > samples( numSamples );
   const real_t          dx = ( xMax - xMin ) / real_c( numSamples - 1 );
   for ( uint_t sample = 0; sample < numSamples; sample++ )
   {
      samples[sample] = c.evaluate( Point3D( {xMin + real_c( sample ) * dx, y, 0} ), level );
   }
   return samples;
}

template < typename FunctionType >
std::vector< real_t > verticalGradientSlice( const FunctionType& c,
                                             uint_t              level,
                                             real_t              xMin,
                                             real_t              xMax,
                                             real_t              yTop,
                                             real_t              yBottom,
                                             uint_t              numSamples )
{
   const auto sliceTop    = evaluateHorizontalTemperatureSlice( c, level, xMin, xMax, yTop, numSamples );
   const auto sliceBottom = evaluateHorizontalTemperatureSlice( c, level, xMin, xMax, yBottom, numSamples );

   const real_t          hInv = 1. / ( yTop - yBottom );
   std::vector< real_t > gradient( numSamples );

   for ( uint_t sample = 0; sample < numSamples; sample++ )
   {
      gradient[sample] = ( sliceTop[sample] - sliceBottom[sample] ) * hInv;
   }
   return gradient;
}

real_t integrateSlice( const std::vector< real_t >& slice, real_t xMin, real_t xMax )
{
   real_t result = 0;
   WALBERLA_CHECK_EQUAL( slice.size() % 2, 1 );
   const auto hThird = ( ( xMax - xMin ) / real_c( slice.size() - 1 ) ) / 3.;
   for ( uint_t x0 = 0; x0 < slice.size() - 2; x0 += 2 )
   {
      auto x1 = x0 + 1;
      auto x2 = x0 + 2;
      result += hThird * ( slice[x0] + 4.0 * slice[x1] + slice[x2] );
   }
   return result;
}

template < typename FunctionType >
real_t calculateNusseltNumber( const FunctionType& c, uint_t level, real_t hGradient, real_t epsBoundary, uint_t numSamples )
{
#if 1
   WALBERLA_CHECK_EQUAL( walberla::mpi::MPIManager::instance()->numProcesses(), 1, "Cannot calculate Nusselt number in parallel." );
   const real_t xMin            = epsBoundary;
   const real_t xMax            = 1.5 - epsBoundary;
   const real_t yGradientTop    = 1.0 - epsBoundary;
   const real_t yGradientBottom = yGradientTop - hGradient;
   const real_t yBottom         = epsBoundary;

   auto gradientSlice         = verticalGradientSlice( c, level, xMin, xMax, yGradientTop, yGradientBottom, numSamples );
   auto bottomSlice           = evaluateHorizontalTemperatureSlice( c, level, xMin, xMax, yBottom, numSamples );
   auto gradientSliceIntegral = integrateSlice( gradientSlice, xMin, xMax );
   auto bottomSliceIntegral   = integrateSlice( bottomSlice, xMin, xMax );

   return -gradientSliceIntegral / bottomSliceIntegral;
#else
   return 0;
#endif
}

std::shared_ptr< SetupPrimitiveStorage > createSetupStorage()
{
   auto meshInfo     = MeshInfo::meshRectangle( Point2D( {0, 0} ), Point2D( {1.5, 1} ), MeshInfo::CROSS, 3, 2 );
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
                   real_t      rayleighNumber,
                   bool        fixedTimeStep,
                   real_t      dtConstant,
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

   const real_t diffusivity     = 1.0;
   const real_t internalHeating = 1.0;

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
   WALBERLA_LOG_INFO_ON_ROOT( "   + Rayleigh number                               " << rayleighNumber )
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
   sqlIntegerProperties["fixed_time_step"]      = int_c( fixedTimeStep );
   sqlRealProperties["cfl_max"]                 = cflMax;
   sqlRealProperties["dt_constant"]             = dtConstant;
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

   WALBERLA_ROOT_SECTION()
   {
      db.storeRun( sqlIntegerProperties, sqlStringProperties, sqlRealProperties );
      sqlRealProperties.clear();
      sqlIntegerProperties.clear();
      sqlStringProperties.clear();
   }

   typedef P2P1TaylorHoodFunction< real_t >                                     StokesFunction;
   typedef P2Function< real_t >                                                 ScalarFunction;
   typedef P2P1TaylorHoodStokesOperator                                         StokesOperator;
   typedef P2ConstantLaplaceOperator                                            LaplaceOperator;
   typedef P2ConstantMassOperator                                               MassOperatorVelocity;
   typedef P1ConstantMassOperator                                               MassOperatorPressure;
   typedef P2ConstantUnsteadyDiffusionOperator                                  UnsteadyDiffusionOperator;
   typedef P2ProjectNormalOperator                                              ProjectNormalOperator;
   typedef StrongFreeSlipWrapper< StokesOperator, ProjectNormalOperator, true > StokesOperatorFreeSlip;

   BoundaryCondition bcVelocity;
   BoundaryCondition bcTemperature;

   bcVelocity.createDirichletBC( "topAndBottomWalls", {1, 2} );
   bcVelocity.createFreeslipBC( "sideWalls", {3, 4} );

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
   StokesFunction upwardNormal( "upwardNormal", storage, level, level, bcVelocity );
   StokesFunction stokesTmp( "stokesTmp", storage, level, level, bcVelocity );
   StokesFunction stokesResidual( "stokesResidual", storage, level, level, bcVelocity );

   ScalarFunction uTmp( "uTmp", storage, level, level, bcVelocity );
   ScalarFunction uTmp2( "uTmp2", storage, level, level, bcVelocity );

   auto initialTemperature = []( const Point3D& x ) {
      return 0.5 * ( 1.0 - x[1] * x[1] ) + 0.01 * std::cos( pi * x[0] / 1.5 ) * std::sin( pi * x[1] / 1.0 );
      // return 0.1 * ( 1.0 - x[1] * x[1] ) + 0.1 * ( 1 - x[0] / 1.5 ) * ( 1.0 - x[1] * x[1] );
      // return 0.2 * ( 1.0 - x[1] * x[1] * x[1] ) + 0.2 * ( 1 - x[0] / 1.5 ) * ( 1.0 - x[1] * x[1] * x[1] );
   };

   c.interpolate( initialTemperature, level, Inner | NeumannBoundary | FreeslipBoundary );
   q.interpolate( internalHeating, level, All );
   upwardNormal.v.interpolate( 1, level, All );

   auto surfaceNormalsFreeSlip = []( const Point3D& in, Point3D& out ) {
      if ( in[0] < 0.75 )
      {
         out = Point3D( {-1, 0, 0} );
      }
      else
      {
         out = Point3D( {1, 0, 0} );
      }
   };

   UnsteadyDiffusionOperator diffusionOperator( storage, level, level, 1.0, diffusivity, DiffusionTimeIntegrator::ImplicitEuler );
   auto projectNormalOperator = std::make_shared< ProjectNormalOperator >( storage, level, level, surfaceNormalsFreeSlip );
   auto A                     = std::make_shared< StokesOperator >( storage, level, level );
   StokesOperatorFreeSlip          AFS( A, projectNormalOperator, FreeslipBoundary );
   LaplaceOperator                 L( storage, level, level );
   MassOperatorVelocity            MVelocity( storage, level, level );
   MassOperatorPressure            MPressure( storage, level, level );
   MMOCTransport< ScalarFunction > transport( storage, setupStorage, level, level, TimeSteppingScheme::RK4 );

#ifndef HYTEG_BUILD_WITH_PETSC
   auto internalDiffusionSolver =
       std::make_shared< PETScMinResSolver< UnsteadyDiffusionOperator > >( storage, level, 1e-15, 50000 );
#else
   auto internalDiffusionSolver = std::make_shared< CGSolver< UnsteadyDiffusionOperator > >( storage, level, level, 5000, 1e-14 );
#endif

#ifdef HYTEG_BUILD_WITH_PETSC
   auto stokesSolver =
       std::make_shared< PETScBlockPreconditionedStokesSolver< StokesOperatorFreeSlip > >( storage, level, 1e-08, 5000, 4, 0 );
   stokesSolver->reassembleMatrix( false );
#else
   auto stokesSolver            = solvertemplates::stokesMinResSolver< StokesOperatorFreeSlip >( storage, level, 1e-14, 5000 );
#endif

   UnsteadyDiffusion< ScalarFunction, UnsteadyDiffusionOperator, LaplaceOperator, MassOperatorVelocity > diffusionSolver(
       storage, level, level, internalDiffusionSolver );

   uint_t numSamples = uint_c( 1.5 / (hMin * 0.97) );
   if ( numSamples % 2 == 0 )
      numSamples--;

   real_t timeTotal = 0;
   real_t vMax      = velocityMaxMagnitude( u.u, u.v, uTmp, uTmp2, level, All );
   real_t nu        = calculateNusseltNumber( c, level, 0.5 * hMin, 0.5 * hMin, numSamples );
   real_t vRms      = 0;
   real_t residualU = 0;
   real_t residualV = 0;
   real_t residualP = 0;

   real_t timeStepTotal = 0;
   real_t timeStokes    = 0;
   real_t timeMMOC      = 0;
   real_t timeDiffusion = 0;

   hyteg::VTKOutput vtkOutput( "./vtk", benchmarkName, storage, vtkInterval );

   vtkOutput.add( u );
   vtkOutput.add( f );
   vtkOutput.add( c );
   vtkOutput.add( upwardNormal );
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

   MVelocity.apply( c, f.u, level, All );
   MVelocity.apply( c, f.v, level, All );
   f.u.multElementwise( {f.u, upwardNormal.u}, level );
   f.v.multElementwise( {f.v, upwardNormal.v}, level );
   f.u.assign( {rayleighNumber}, {f.u}, level, All );
   f.v.assign( {rayleighNumber}, {f.v}, level, All );
   projectNormalOperator->apply( f, level, FreeslipBoundary );
   stokesSolver->solve( AFS, u, f, level );

   if ( vtk )
      vtkOutput.write( level );

   walberla::WcTimer timeStepTimer;

   while ( timeTotal < simulationTime )
   {
      timeStepTimer.start();

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
      projectNormalOperator->apply( f, level, FreeslipBoundary );

      localTimer.start();
      stokesSolver->solve( AFS, u, f, level );
      localTimer.end();
      timeStokes = localTimer.last();

      calculateStokesResiduals(
          AFS, MVelocity, MPressure, u, f, level, stokesResidual, stokesTmp, residualU, residualV, residualP );
      vRms = velocityRMS( u, 1.0, 1.5, level );

      // Energy

      vMax = velocityMaxMagnitude( u.u, u.v, uTmp, uTmp2, level, All );

      real_t dt;
      if ( fixedTimeStep )
      {
         dt = dtConstant;
      }
      else
      {
         dt = ( cflMax / vMax ) * hMin;
      }

      if ( verbose )
         WALBERLA_LOG_INFO_ON_ROOT( "performing transport step" )

      if ( adjustedAdvection )
      {
         const real_t adjustedAdvectionPertubation = 0.1 * ( hMin / vMax );
         localTimer.start();
         transport.step(
             c, u.u, u.v, u.w, uLast.u, uLast.v, uLast.w, level, All, dt, 1, MVelocity, 0.0, adjustedAdvectionPertubation );
         localTimer.end();
         timeMMOC = localTimer.last();
      }
      else
      {
         localTimer.start();
         transport.step( c, u.u, u.v, u.w, uLast.u, uLast.v, uLast.w, level, All, dt, 1, true );
         localTimer.end();
         timeMMOC = localTimer.last();
      }

      c.interpolate( 0, level, DirichletBoundary );

      cOld.assign( {1.0}, {c}, level, All );

      timeTotal += dt;

      diffusionOperator.setDt( dt );

      if ( verbose )
         WALBERLA_LOG_INFO_ON_ROOT( "performing diffusion step" )

      localTimer.start();
      diffusionSolver.step( diffusionOperator, L, MVelocity, c, cOld, q, q, level, Inner | NeumannBoundary | FreeslipBoundary );
      localTimer.end();
      timeDiffusion = localTimer.last();

      nu = calculateNusseltNumber( c, level, 0.5 * hMin, 0.5 * hMin, numSamples );

      if ( vtk )
         vtkOutput.write( level, timeStep );

      WALBERLA_ROOT_SECTION()
      {
         sqlIntegerProperties["ts"]              = int_c( timeStep );
         sqlRealProperties["sim_time"]           = timeTotal;
         sqlRealProperties["run_time_advection"] = timeMMOC;
         sqlRealProperties["run_time_stokes"]    = timeStokes;
         sqlRealProperties["run_time_time_step"] = timeStepTotal;
         sqlRealProperties["v_max"]              = vMax;
         sqlRealProperties["v_rms"]              = vRms;
         sqlRealProperties["nu"]                 = nu;
         sqlRealProperties["dt"]                 = dt;

         db.storeRun( sqlIntegerProperties, sqlStringProperties, sqlRealProperties );
         sqlRealProperties.clear();
         sqlIntegerProperties.clear();
         sqlStringProperties.clear();
      }

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
      auto defaultFile = "./Blankenbach_Case3.prm";
      WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   const hyteg::real_t cflMax         = mainConf.getParameter< hyteg::real_t >( "cflMax" );
   const hyteg::real_t rayleighNumber = mainConf.getParameter< hyteg::real_t >( "rayleighNumber" );
   const bool          fixedTimeStep  = mainConf.getParameter< bool >( "fixedTimeStep" );
   const hyteg::real_t dtConstant     = mainConf.getParameter< hyteg::real_t >( "dtConstant" );
   const uint_t        level          = mainConf.getParameter< uint_t >( "level" );
   const hyteg::real_t simulationTime = mainConf.getParameter< hyteg::real_t >( "simulationTime" );

   hyteg::runBenchmark(
       cflMax, rayleighNumber, fixedTimeStep, dtConstant, level, false, simulationTime, true, 1, 1, false, "database.db" );
}