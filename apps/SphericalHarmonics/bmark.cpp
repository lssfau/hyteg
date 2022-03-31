/*
 * Copyright (c) 2020 Marcus Mohr
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

/* =====================================================================
This app implements the Stokes benchmark described in:

@article{Horbach:2020:GEM,
  author = {Andr{\'e} Horbach and Marcus Mohr and Hans-Peter Bunge},
  journal = {GEM - International Journal of Geomathematics},
  number = {1},
  title = {{A Semi-Analytic Accuracy Benchmark for Stokes Flow in
           3-D Spherical Mantle Convection Codes}},
  volume = {11},
  year = {2020},
  doi = {10.1007/s13137-019-0137-3},
}
===================================================================== */

#include <cmath>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/HytegDefinitions.hpp"
#include "hyteg/composites/StrongFreeSlipWrapper.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/numerictools/SphericalHarmonicsTool.hpp"
#include "hyteg/p1functionspace/P1ProjectNormalOperator.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/petsc/PETScBlockPreconditionedStokesSolver.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"

using walberla::real_c;
using walberla::real_t;
using namespace hyteg;

#include "bmark_tools.hpp"

// radii of spherical shell
const real_t innerRadius = 1.0;
const real_t outerRadius = 2.0;

namespace terraneo {

// =========
//  Meshing
// =========
std::shared_ptr< hyteg::PrimitiveStorage > generateMesh( uint_t nRad, uint_t nTan, bool reportPrimitives )
{
   hyteg::MeshInfo              meshInfo = hyteg::MeshInfo::meshSphericalShell( nTan, nRad, innerRadius, outerRadius );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                              walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   hyteg::loadbalancing::roundRobin( setupStorage );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   IcosahedralShellMap::setMap( setupStorage );

   auto surface = []( const Point3D& p ) { return std::abs( p.norm() - outerRadius ) < 1e-10; };
   auto cmb     = []( const Point3D& p ) { return std::abs( p.norm() - innerRadius ) < 1e-10; };
   setupStorage.setMeshBoundaryFlagsByVertexLocation( 1, surface );
   setupStorage.setMeshBoundaryFlagsByVertexLocation( 2, cmb );

   if ( reportPrimitives )
      WALBERLA_LOG_INFO_ON_ROOT( "" << setupStorage );

   std::shared_ptr< walberla::WcTimingTree >  timingTree( new walberla::WcTimingTree() );
   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, timingTree );

   return storage;
}

// ===================
//  runBenchmarkTests
// ===================
template < template < class > class feFuncType, typename massOpType, typename stokesOpType, typename projOpType >
void runBenchmarkTests( std::shared_ptr< walberla::config::Config > cfg,
                        std::shared_ptr< hyteg::PrimitiveStorage >  storage,
                        uint_t                                      targetLevel )
{
   // parameter handling
   const walberla::Config::BlockHandle problemCfg = cfg->getBlock( "Problem" );

   // as long as we solve with PETSc we just need one level
   const uint_t minLevel = targetLevel;
   const uint_t maxLevel = targetLevel;

   std::string       bc = problemCfg.getParameter< std::string >( "bcType" );
   bcType            bmarkBC;
   DoFType           bmarkFsDoFs;
   BoundaryCondition bcVelocity;

   if ( bc == "no-slip" )
   {
      bmarkBC     = BC_NOSLIP;
      bmarkFsDoFs = None;
      bcVelocity.createDirichletBC( "all boundaries", {1, 2} );
   }
   else if ( bc == "mixed" )
   {
      bmarkBC     = BC_MIXED;
      bmarkFsDoFs = FreeslipBoundary;
      bcVelocity.createDirichletBC( "surface", {1} );
      bcVelocity.createFreeslipBC( "CMB", {2} );
   }
   else if ( bc == "free-slip" )
   {
      bmarkBC     = BC_FREESLIP;
      bmarkFsDoFs = FreeslipBoundary;
      bcVelocity.createFreeslipBC( "all boundaries", {1, 2} );
      WALBERLA_ABORT( "Unsupported bcType = '" << bc << "' detected!" );
   }
   else
   {
      WALBERLA_ABORT( "Unsupported bcType = '" << bc << "' detected!" );
   }

   // prepare spherical harmonics computations
   const uint_t                              degree  = problemCfg.getParameter< uint_t >( "degree" );
   const int                                 order   = problemCfg.getParameter< int >( "order" );
   const uint_t                              lmax    = degree;
   std::shared_ptr< SphericalHarmonicsTool > sphTool = std::make_shared< SphericalHarmonicsTool >( lmax );

   // some short-hands
   typedef typename feFuncType< real_t >::VelocityFunction_T vf_t;
   typedef typename feFuncType< real_t >::PressureFunction_T pf_t;

   // determine number of degrees of freedom
   if ( problemCfg.getParameter< bool >( "reportDoFNum" ) )
   {
      uint_t pressureDoFs = numberOfGlobalDoFs< typename pf_t::Tag >( *storage, targetLevel );
      WALBERLA_LOG_INFO_ON_ROOT( "Number of pressure DoFs on target level = " << pressureDoFs );
      uint_t velocityDoFs = numberOfGlobalDoFs< typename vf_t::Tag >( *storage, targetLevel );
      WALBERLA_LOG_INFO_ON_ROOT( "Number of velocity DoFs on target level = " << velocityDoFs );
   }

   // ------------------------
   //  (Weak) right-hand side
   // ------------------------

   // set the strong right-hand side
   WALBERLA_LOG_PROGRESS_ON_ROOT( "Computing forcing field ..." );
   vf_t forcing( "forcing", storage, minLevel, maxLevel );
   setForceField( forcing, maxLevel, degree, order, sphTool, bmarkBC );

   // derive from this the weak right-hand side
   WALBERLA_LOG_PROGRESS_ON_ROOT( "Computing weak right-hand side ..." );
   massOpType           mass( storage, maxLevel, maxLevel );
   feFuncType< real_t > rhs( "right-hand side", storage, minLevel, maxLevel );
   mass.apply( forcing[0], rhs.uvw()[0], maxLevel, All );
   mass.apply( forcing[1], rhs.uvw()[1], maxLevel, All );
   mass.apply( forcing[2], rhs.uvw()[2], maxLevel, All );

   // ---------------------
   //  Analytical Solution
   // ---------------------
   WALBERLA_LOG_PROGRESS_ON_ROOT( "Interpolating analytical solution to FEspace ..." );
   vf_t u_exact( "Analytical Solution (mapped to FE space)", storage, maxLevel, maxLevel );
   setAnalyticSolution( u_exact, maxLevel, degree, order, sphTool, bmarkBC );

   // -------------------
   //  Setup FE Solution
   // -------------------
   WALBERLA_LOG_PROGRESS_ON_ROOT( "Preparing FE function for discrete solution ... " );
   feFuncType< real_t > feSol( "Discrete Solution", storage, minLevel, maxLevel, bcVelocity );

   // set everything to zero (including boundaries)
   feSol.uvw()[0].setToZero( maxLevel );
   feSol.uvw()[1].setToZero( maxLevel );
   feSol.uvw()[2].setToZero( maxLevel );
   feSol.p().setToZero( maxLevel );

   // ---------------------
   //  Stuff for Free-Slip
   // ---------------------
   auto normalFunc = []( const Point3D& p, Point3D& n ) -> void {
      real_t norm = p.norm();
      real_t sign = ( norm > 0.95 * outerRadius ) ? 1.0 : -1.0;
      n           = sign / norm * p;
   };
   auto projectNormalsOp        = std::make_shared< projOpType >( storage, maxLevel, maxLevel, normalFunc );
   using StokesOperatorFreeSlip = StrongFreeSlipWrapper< stokesOpType, projOpType >;

   // -----------------------
   //  Setup Stokes Operator
   // -----------------------
   WALBERLA_LOG_PROGRESS_ON_ROOT( "Preparing Stokes operator ... " );
   std::shared_ptr< stokesOpType > stokesOp = std::make_shared< stokesOpType >( storage, maxLevel, maxLevel );
   StokesOperatorFreeSlip          stokesOpFS( stokesOp, projectNormalsOp, bmarkFsDoFs );
   WALBERLA_LOG_PROGRESS_ON_ROOT( "Finished with Stokes operator ... " );

   // ---------------------
   //  Solve linear system
   // ---------------------
   std::string solverType = problemCfg.getParameter< std::string >( "solverType" );

   switch ( bmarkBC )
   {
   case BC_NOSLIP:
   {
      if ( solverType.compare( "PETScLU" ) == 0 )
      {
#ifdef HYTEG_BUILD_WITH_PETSC
         WALBERLA_LOG_PROGRESS_ON_ROOT( "Sparse direct solution with PETSc ... " );
         PETScManager                  petscManager;
         PETScLUSolver< stokesOpType > PETScLU( storage, maxLevel );
         PETScLU.solve( *stokesOp, feSol, rhs, maxLevel );
#else
         WALBERLA_ABORT( "Recompile with PETSc support to use PETSc solvers!" );
#endif
      }
      else if ( solverType.compare( "PETScBPSS" ) == 0 )
      {
#ifdef HYTEG_BUILD_WITH_PETSC
         WALBERLA_LOG_PROGRESS_ON_ROOT( "Solving iteratively with PETScBPSS ... " );
         PETScManager                                         petscManager;
         PETScBlockPreconditionedStokesSolver< stokesOpType > stokesSolver( storage, maxLevel, 1e-08, 5000, 1, 1 );
         stokesSolver.setVerbose( true );
         stokesSolver.solve( *stokesOp, feSol, rhs, maxLevel );
#else
         WALBERLA_ABORT( "Recompile with PETSc support to use PETSc solvers!" );
#endif
      }
      else if ( solverType.compare( "Minres" ) == 0 )
      {
         auto stokesSolver = solvertemplates::stokesMinResSolver< StokesOperatorFreeSlip >( storage, maxLevel, 1e-08, 200, true );
         WALBERLA_LOG_PROGRESS_ON_ROOT( "Iterative solution with MINRES ... " );
         stokesSolver->solve( stokesOpFS, feSol, rhs, maxLevel );
      }
      else
      {
         WALBERLA_ABORT( "Solver for no-slip '" << solverType << "' unsupported!" );
      }
      break;
   }
   case BC_MIXED:
   {
      // PETScMinResSolver< StokesOperatorFreeSlip > stokesSolver( storage, maxLevel, 1e-30, 1e-08, 200 );
      // PETScBlockPreconditionedStokesSolver< StokesOperatorFreeSlip > stokesSolver( storage, maxLevel, 1e-08, 5000, 4, 1 );
      auto stokesSolver = solvertemplates::stokesMinResSolver< StokesOperatorFreeSlip >( storage, maxLevel, 1e-08, 200, true );
      WALBERLA_LOG_PROGRESS_ON_ROOT( "Iterative solution with MINRES ... " );
      stokesSolver->solve( stokesOpFS, feSol, rhs, maxLevel );
      break;
   }
   case BC_FREESLIP:
      WALBERLA_ABORT( "No Solver for BC_FREESLIP, yes!" );
   }

   // -----------------
   //  Determine error
   // -----------------
   vf_t error( "error", storage, maxLevel, maxLevel );
   error.assign( {1.0, -1.0}, {u_exact, feSol.uvw()}, maxLevel );

   if ( problemCfg.getParameter< bool >( "compute_L2_error" ) )
   {
      real_t                             L2error = 0.0;
      typename vf_t::VectorComponentType aux( "aux", storage, maxLevel, maxLevel );

      for ( uint_t idx = 0; idx < 3; ++idx )
      {
         mass.apply( error[idx], aux, maxLevel, All );
         L2error += error[idx].dotGlobal( aux, maxLevel );
      }
      L2error = std::sqrt( L2error );
      WALBERLA_LOG_INFO_ON_ROOT( "L2 error = " << std::scientific << L2error );
   }

   // We estimate the L_infty error on the next finer mesh
   if ( problemCfg.getParameter< bool >( "compute_max_error" ) )
   {
      real_t maxError = 0.0;
      uint_t fLevel   = maxLevel + 1;
      vf_t   u_fine( "u_fine", storage, fLevel, fLevel );
      setAnalyticSolution( u_fine, fLevel, degree, order, sphTool, bmarkBC );

      typename vf_t::VectorComponentType aux( "aux", storage, maxLevel, maxLevel + 1 );
      P2toP2QuadraticProlongation        embedder;

      for ( uint_t idx = 0; idx < 3; ++idx )
      {
         aux.assign( {1.0}, {feSol.uvw()[idx]}, maxLevel );
         embedder.prolongate( aux, maxLevel, All );
         aux.assign( {1.0, -1.0}, {u_fine[idx], aux}, fLevel );
         real_t val = aux.getMaxMagnitude( fLevel );
         maxError   = maxError < val ? val : maxError;
      }

      WALBERLA_LOG_INFO_ON_ROOT( "Maximum error = " << std::scientific << maxError );
   }

   // -----------------------------
   //  Output stuff for inspection
   // -----------------------------
   if ( problemCfg.getParameter< bool >( "VTKOutput" ) )
   {
      hyteg::VTKOutput vtkOutput( "./output", "bmark", storage );
      vtkOutput.add( forcing );
      vtkOutput.add( rhs.uvw() );
      vtkOutput.add( u_exact );
      vtkOutput.add( feSol.uvw() );
      vtkOutput.add( error );
      vtkOutput.write( maxLevel, 0 );
   }
}
} // namespace terraneo

// Hunting the NaN
#include <cfenv>

int main( int argc, char* argv[] )
{
#ifndef __APPLE__
   #ifndef _MSC_VER
      feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );
   #endif
#endif

   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );

   // ============
   //  Parameters
   // ============

   // check if a config was given on command line or load default file otherwise
   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./bmark.prm";
      WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle problemCfg = cfg->getBlock( "Problem" );
   if ( walberla::MPIManager::instance()->worldRank() == 0 )
   {
      problemCfg.listParameters();
   }

   // =========
   //  Meshing
   // =========
   const uint_t nRad    = problemCfg.getParameter< uint_t >( "nRad" );
   const uint_t nTan    = problemCfg.getParameter< uint_t >( "nTan" );
   auto         storage = terraneo::generateMesh( nRad, nTan, problemCfg.getParameter< bool >( "reportPrimitives" ) );

   // =================
   //  Benchmark Stuff
   // =================
   uint_t bmInitLevel = problemCfg.getParameter< uint_t >( "bmInitLevel" );
   uint_t bmStopLevel = problemCfg.getParameter< uint_t >( "bmStopLevel" );

   // delegate actual work to templated function to allow different FE spaces
   std::string feSpace = problemCfg.getParameter< std::string >( "feSpace" );

   for ( uint_t bmCurrentLevel = bmInitLevel; bmCurrentLevel <= bmStopLevel; bmCurrentLevel++ )
   {
      if ( feSpace == "P2P1TaylorHood" )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "=================================================" );
         WALBERLA_LOG_INFO_ON_ROOT( " REFINEMENT LEVEL = " << bmCurrentLevel );
         WALBERLA_LOG_INFO_ON_ROOT( "=================================================" );
         terraneo::runBenchmarkTests< P2P1TaylorHoodFunction,
                                      P2ElementwiseBlendingMassOperator,
                                      P2P1ElementwiseBlendingStokesOperator,
                                      P2ProjectNormalOperator >( cfg, storage, bmCurrentLevel );
      }
   }

   return EXIT_SUCCESS;
}
