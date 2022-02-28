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

#include "hyteg/Git.hpp"
#include "hyteg/MeshQuality.hpp"
#include "hyteg/composites/StrongFreeSlipWrapper.hpp"
#include "hyteg/composites/UnsteadyDiffusion.hpp"
#include "hyteg/dataexport/SQL.hpp"
#include "hyteg/dataexport/TimingOutput.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
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
   residualU = normL2( r.uvw()[0], tmp.uvw()[0], Mu, level, Inner | NeumannBoundary | FreeslipBoundary );
   residualV = normL2( r.uvw()[1], tmp.uvw()[1], Mu, level, Inner | NeumannBoundary | FreeslipBoundary );
   residualP = normL2( r.p(), tmp.p(), Mp, level, Inner | NeumannBoundary | FreeslipBoundary );
}

template < typename StokesFunction, typename VelocityMass >
real_t velocityRMS( const StokesFunction& u,
                    const StokesFunction& tmp,
                    const VelocityMass&   M,
                    real_t                domainHeight,
                    real_t                domainWidth,
                    uint_t                level )
{
   auto norm = std::pow( normL2( u.uvw()[0], tmp.uvw()[0], M, level, All ), 2.0 ) +
               std::pow( normL2( u.uvw()[1], tmp.uvw()[1], M, level, All ), 2.0 );
   const auto area = domainHeight * domainWidth;
   return std::sqrt( norm / area );
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
   std::vector< bool >   sampleLocallyAvailable( numSamples, false );
   const real_t          dx = ( xMax - xMin ) / real_c( numSamples - 1 );
   for ( uint_t sample = 0; sample < numSamples; sample++ )
   {
      Point3D pos( {xMin + real_c( sample ) * dx, y, 0} );
      sampleLocallyAvailable[sample] = c.evaluate( pos, level, samples[sample], 1e-5 );
   }

   walberla::mpi::SendBuffer sendbuffer;
   walberla::mpi::RecvBuffer recvbuffer;

   sendbuffer << samples;
   sendbuffer << sampleLocallyAvailable;

   walberla::mpi::allGathervBuffer( sendbuffer, recvbuffer );

   std::vector< real_t > samplesGlobal( numSamples );
   std::vector< bool >   sampleGloballyAvailable( numSamples, false );

   for ( uint_t rank = 0; rank < uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ); rank++ )
   {
      std::vector< real_t > recvSamples( numSamples );
      std::vector< bool >   recvSamplesAvailable( numSamples, false );
      recvbuffer >> recvSamples;
      recvbuffer >> recvSamplesAvailable;

      for ( uint_t sample = 0; sample < numSamples; sample++ )
      {
         if ( recvSamplesAvailable[sample] )
         {
            sampleGloballyAvailable[sample] = true;
            samplesGlobal[sample]           = recvSamples[sample];
         }
      }
   }

   for ( uint_t sample = 0; sample < numSamples; sample++ )
   {
      WALBERLA_CHECK( sampleGloballyAvailable[sample] );
   }

   return samplesGlobal;
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
}

std::shared_ptr< SetupPrimitiveStorage > createSetupStorage( uint_t nx )
{
   WALBERLA_CHECK_EQUAL( nx % 3, 0 );

   auto meshInfo     = MeshInfo::meshRectangle( Point2D( {0, 0} ), Point2D( {1.5, 1} ), MeshInfo::CROSS, nx, ( nx / 3 ) * 2 );
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
                   uint_t      nx,
                   bool        predictorCorrector,
                   real_t      simulationTime,
                   bool        vtk,
                   uint_t      printInterval,
                   uint_t      vtkInterval,
                   bool        strangSplitting,
                   std::string dbFile )
{
   walberla::WcTimer localTimer;

   auto gitHash = gitSHA1();
   printGitInfo();

   const std::string benchmarkName    = "Blankenbach_Case3";
   const bool        outputTimingJSON = true;

   auto setupStorage = createSetupStorage( nx );
   auto storage      = std::make_shared< PrimitiveStorage >( *setupStorage, 1 );

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
   WALBERLA_LOG_INFO_ON_ROOT( "   + Strang splitting:                             " << ( strangSplitting ? "yes" : "no" ) )
   WALBERLA_LOG_INFO_ON_ROOT( "   + simulation time:                              " << simulationTime )
   WALBERLA_LOG_INFO_ON_ROOT( " - space discretization: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + dimensions:                                   " << ( storage->hasGlobalCells() ? "3" : "2" ) )
   WALBERLA_LOG_INFO_ON_ROOT( "   + level:                                        " << level )
   WALBERLA_LOG_INFO_ON_ROOT( "   + unknowns (== particles), including boundary:  " << unknowns )
   WALBERLA_LOG_INFO_ON_ROOT( "   + h_min:                                        " << hMin )
   WALBERLA_LOG_INFO_ON_ROOT( "   + h_max:                                        " << hMax )
   WALBERLA_LOG_INFO_ON_ROOT( " - benchmark settings: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + Rayleigh number                               " << rayleighNumber )
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

   db.setConstantEntry( "git_sha1", gitHash );
   db.setConstantEntry( "fixed_time_step", fixedTimeStep );
   db.setConstantEntry( "cfl_max", cflMax );
   db.setConstantEntry( "dt_constant", dtConstant );
   db.setConstantEntry( "simulation_time", simulationTime );
   db.setConstantEntry( "level", level );
   db.setConstantEntry( "unknowns", unknowns );
   db.setConstantEntry( "h_min", hMin );
   db.setConstantEntry( "h_max", hMax );
   db.setConstantEntry( "num_macro_cells", setupStorage->getNumberOfCells() );
   db.setConstantEntry( "num_macro_faces", setupStorage->getNumberOfFaces() );
   db.setConstantEntry( "num_macro_edges", setupStorage->getNumberOfEdges() );
   db.setConstantEntry( "num_macro_vertices", setupStorage->getNumberOfVertices() );
   db.setConstantEntry( "num_macro_primitives", setupStorage->getNumberOfPrimitives() );
   db.setConstantEntry( "diffusivity", diffusivity );
   db.setConstantEntry( "strang_splitting", strangSplitting );

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
   ScalarFunction cPr( "cPr", storage, level, level, bcTemperature );
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
   };

   c.interpolate( initialTemperature, level, Inner | NeumannBoundary | FreeslipBoundary );
   q.interpolate( internalHeating, level, All );
   upwardNormal.uvw()[1].interpolate( 1, level, All );

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

   UnsteadyDiffusionOperator diffusionOperator( storage, level, level, 1.0, diffusivity, DiffusionTimeIntegrator::CrankNicolson );
   auto projectNormalOperator = std::make_shared< ProjectNormalOperator >( storage, level, level, surfaceNormalsFreeSlip );
   auto A                     = std::make_shared< StokesOperator >( storage, level, level );
   StokesOperatorFreeSlip          AFS( A, projectNormalOperator, FreeslipBoundary );
   LaplaceOperator                 L( storage, level, level );
   MassOperatorVelocity            MVelocity( storage, level, level );
   MassOperatorPressure            MPressure( storage, level, level );
   MMOCTransport< ScalarFunction > transport( storage, level, level, TimeSteppingScheme::RK4 );

   auto internalDiffusionSolver = std::make_shared< CGSolver< UnsteadyDiffusionOperator > >( storage, level, level, 5000, 1e-14 );

#ifdef HYTEG_BUILD_WITH_PETSC
   auto stokesSolver =
       std::make_shared< PETScBlockPreconditionedStokesSolver< StokesOperatorFreeSlip > >( storage, level, 1e-08, 5000, 4, 1 );
   stokesSolver->reassembleMatrix( false );
#else
   auto stokesSolver = solvertemplates::stokesMinResSolver< StokesOperatorFreeSlip >( storage, level, 1e-14, 5000 );
#endif

   UnsteadyDiffusion< ScalarFunction, UnsteadyDiffusionOperator, LaplaceOperator, MassOperatorVelocity > diffusionSolver(
       storage, level, level, internalDiffusionSolver );

   uint_t numSamples = uint_c( 1.5 / ( hMin * 0.97 ) );
   if ( numSamples % 2 == 0 )
      numSamples--;

   const real_t nusseltEpsBoundary = 1e-12;
   const real_t nusseltDiffH       = 2.5e-1 * hMin;

   real_t timeTotal = 0;
   real_t vMax      = velocityMaxMagnitude( u.uvw(), uTmp, uTmp2, level, All );
   real_t nu        = 0;
   real_t vRms      = 0;
   real_t residualU = 0;
   real_t residualV = 0;
   real_t residualP = 0;

   real_t timeStepTotal = 0;
   real_t timeStokes    = 0;
   real_t timeMMOC      = 0;
   real_t timeDiffusion = 0;

   VTKOutput vtkOutput( "./vtk", benchmarkName, storage, vtkInterval );

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

   MVelocity.apply( c, f.uvw()[0], level, All );
   MVelocity.apply( c, f.uvw()[1], level, All );
   f.uvw()[0].multElementwise( {f.uvw()[0], upwardNormal.uvw()[0]}, level );
   f.uvw()[1].multElementwise( {f.uvw()[1], upwardNormal.uvw()[1]}, level );
   f.uvw()[0].assign( {rayleighNumber}, {f.uvw()[0]}, level, All );
   f.uvw()[1].assign( {rayleighNumber}, {f.uvw()[1]}, level, All );
   projectNormalOperator->project( f, level, FreeslipBoundary );

   localTimer.start();
   stokesSolver->solve( AFS, u, f, level );
   localTimer.end();
   timeStokes = localTimer.last();

   vMax = velocityMaxMagnitude( u.uvw(), uTmp, uTmp2, level, All );

   if ( vtk )
   {
      vtkOutput.write( level );
   }

   walberla::WcTimer timeStepTimer;

   vRms = velocityRMS( u, stokesTmp, MVelocity, 1.0, 1.5, level );
   nu   = calculateNusseltNumber( c, level, nusseltDiffH, nusseltEpsBoundary, numSamples );

   db.setVariableEntry( "ts", uint_c( 0 ) );
   db.setVariableEntry( "sim_time", timeTotal );
   db.setVariableEntry( "run_time_advection", timeMMOC );
   db.setVariableEntry( "run_time_stokes", timeStokes );
   db.setVariableEntry( "run_time_time_step", timeStepTotal );
   db.setVariableEntry( "v_max", vMax );
   db.setVariableEntry( "v_rms", vRms );
   db.setVariableEntry( "nu", nu );
   db.setVariableEntry( "dt", real_c( 0 ) );

   db.writeRowOnRoot();

   while ( timeTotal < simulationTime )
   {
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
      cPr.assign( {1.0}, {c}, level, All );

      // let's just use the current velocity for the prediction
      uLast.assign( {1.0}, {u}, level, All );

      // diffusion

      cPr.interpolate( 0, level, DirichletBoundary );

      cOld.assign( {1.0}, {cPr}, level, All );

      if ( strangSplitting )
      {
         diffusionOperator.setDt( 0.5 * dt );
      }
      else
      {
         diffusionOperator.setDt( dt );
      }

      timeDiffusion = 0;
      if ( strangSplitting )
      {
         localTimer.start();
         diffusionSolver.step(
             diffusionOperator, L, MVelocity, cPr, cOld, q, q, level, Inner | NeumannBoundary | FreeslipBoundary );
         localTimer.end();
         timeDiffusion += localTimer.last();
      }

      localTimer.start();
      transport.step( cPr, u.uvw(), uLast.uvw(), level, All, dt, 1, true );
      localTimer.end();
      timeMMOC = localTimer.last();

      // diffusion

      cPr.interpolate( 0, level, DirichletBoundary );

      cOld.assign( {1.0}, {cPr}, level, All );

      localTimer.start();
      diffusionSolver.step( diffusionOperator, L, MVelocity, cPr, cOld, q, q, level, Inner | NeumannBoundary | FreeslipBoundary );
      localTimer.end();
      timeDiffusion += localTimer.last();

      // Stokes

      MVelocity.apply( cPr, f.uvw()[0], level, All );
      MVelocity.apply( cPr, f.uvw()[1], level, All );
      f.uvw()[0].multElementwise( {f.uvw()[0], upwardNormal.uvw()[0]}, level );
      f.uvw()[1].multElementwise( {f.uvw()[1], upwardNormal.uvw()[1]}, level );
      f.uvw()[0].assign( {rayleighNumber}, {f.uvw()[0]}, level, All );
      f.uvw()[1].assign( {rayleighNumber}, {f.uvw()[1]}, level, All );
      projectNormalOperator->project( f, level, FreeslipBoundary );

      localTimer.start();
      stokesSolver->solve( AFS, u, f, level );
      localTimer.end();
      timeStokes = localTimer.last();

      calculateStokesResiduals(
          AFS, MVelocity, MPressure, u, f, level, stokesResidual, stokesTmp, residualU, residualV, residualP );

      if ( predictorCorrector )
      {
         // energy

         // diffusion

         c.interpolate( 0, level, DirichletBoundary );

         cOld.assign( {1.0}, {c}, level, All );

         if ( strangSplitting )
         {
            localTimer.start();
            diffusionSolver.step(
                diffusionOperator, L, MVelocity, c, cOld, q, q, level, Inner | NeumannBoundary | FreeslipBoundary );
            localTimer.end();
            timeDiffusion += localTimer.last();
         }

         // advection

         localTimer.start();
         transport.step( c, u.uvw(), uLast.uvw(), level, All, dt, 1, true );
         localTimer.end();
         timeMMOC += localTimer.last();

         // diffusion

         c.interpolate( 0, level, DirichletBoundary );

         cOld.assign( {1.0}, {c}, level, All );

         localTimer.start();
         diffusionSolver.step(
             diffusionOperator, L, MVelocity, c, cOld, q, q, level, Inner | NeumannBoundary | FreeslipBoundary );
         localTimer.end();
         timeDiffusion += localTimer.last();

         // Stokes

         MVelocity.apply( c, f.uvw()[0], level, All );
         MVelocity.apply( c, f.uvw()[1], level, All );
         f.uvw()[0].multElementwise( {f.uvw()[0], upwardNormal.uvw()[0]}, level );
         f.uvw()[1].multElementwise( {f.uvw()[1], upwardNormal.uvw()[1]}, level );
         f.uvw()[0].assign( {rayleighNumber}, {f.uvw()[0]}, level, All );
         f.uvw()[1].assign( {rayleighNumber}, {f.uvw()[1]}, level, All );
         projectNormalOperator->project( f, level, FreeslipBoundary );

         localTimer.start();
         stokesSolver->solve( AFS, u, f, level );
         localTimer.end();
         timeStokes += localTimer.last();

         calculateStokesResiduals(
             AFS, MVelocity, MPressure, u, f, level, stokesResidual, stokesTmp, residualU, residualV, residualP );
      }
      else
      {
         // use predicted value
         c.assign( {1.0}, {cPr}, level, All );
      }

      timeTotal += dt;

      vRms = velocityRMS( u, stokesTmp, MVelocity, 1.0, 1.5, level );
      nu   = calculateNusseltNumber( c, level, nusseltDiffH, nusseltEpsBoundary, numSamples );

      if ( vtk )
      {
         vtkOutput.write( level, timeStep );
      }

      db.setVariableEntry( "ts", timeStep );
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
      auto defaultFile = "./Blankenbach_Case3.prm";
      WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   const hyteg::real_t cflMax          = mainConf.getParameter< hyteg::real_t >( "cflMax" );
   const hyteg::real_t rayleighNumber  = mainConf.getParameter< hyteg::real_t >( "rayleighNumber" );
   const bool          fixedTimeStep   = mainConf.getParameter< bool >( "fixedTimeStep" );
   const hyteg::real_t dtConstant      = mainConf.getParameter< hyteg::real_t >( "dtConstant" );
   const uint_t        level           = mainConf.getParameter< uint_t >( "level" );
   const hyteg::real_t simulationTime  = mainConf.getParameter< hyteg::real_t >( "simulationTime" );
   const hyteg::uint_t nx              = mainConf.getParameter< hyteg::uint_t >( "nx" );
   const std::string   dbFile          = mainConf.getParameter< std::string >( "dbFile" );
   const bool          vtk             = mainConf.getParameter< bool >( "vtk" );
   const bool          strangSplitting = mainConf.getParameter< bool >( "strangSplitting", true );

   hyteg::runBenchmark(
       cflMax, rayleighNumber, fixedTimeStep, dtConstant, level, nx, true, simulationTime, vtk, 1, 1, strangSplitting, dbFile );
}
