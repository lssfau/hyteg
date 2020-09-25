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

#include "hyteg/FunctionProperties.hpp"
#include "hyteg/HytegDefinitions.hpp"
#include "hyteg/composites/StrongFreeSlipWrapper.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/numerictools/SphericalHarmonicsTool.hpp"
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

// radii of spherical shell
const real_t innerRadius = 1.0;
const real_t outerRadius = 2.0;

namespace terraneo {

// type to distinguish different scenarios w.r.t. boundary conditions
typedef enum
{
   BC_NOSLIP,
   BC_FREESLIP,
   BC_MIXED
} bcType;

// ===============
//  SetForceField
// ===============
template < typename vecFuncType >
void setForceField( vecFuncType&                              f,
                    uint_t                                    level,
                    uint_t                                    degree,
                    int                                       order,
                    std::shared_ptr< SphericalHarmonicsTool > sphTool,
                    bcType                                    bmCase )
{
   // safety check
   if ( (uint_t) std::abs( order ) > degree )
   {
      WALBERLA_ABORT( "Spherical harmonics order must be smaller equal to degree!" );
   }

   // first we compute the spherical harmonics contribution which is the same
   // for all three component functions
   std::function< real_t( const Point3D& ) > sphFunc = [sphTool, degree, order]( const Point3D& x ) {
      return sphTool->shconvert_eval( degree, order, x[0], x[1], x[2] );
   };

   typename vecFuncType::template VectorComponentType sph( "spherical harmonic", f[0].getStorage(), level, level );
   // feFuncType sph( "spherical harmonic", f[0].getStorage(), level, level );
   sph.interpolate( sphFunc, level, All );

   // set polynomial coefficients for radial component
   real_t                  elfac = real_c( degree * ( degree + 1 ) );
   std::array< real_t, 5 > pc;

   pc[0] = ( 12.0 - elfac ) * ( 2.0 - elfac );
   pc[1] = elfac * ( 6.0 - elfac );
   pc[2] = elfac * ( elfac - 2.0 );
   pc[3] = elfac * ( 2.0 - elfac );
   pc[4] = elfac * ( elfac - 6.0 );

   switch ( bmCase )
   {
   case BC_NOSLIP:
      pc[0] *= 0.25;
      pc[1] *= 1.50;
      pc[2] *= 3.25;
      pc[3] *= 3.00;
      pc[4] *= 1.00;
      break;

   case BC_FREESLIP:
      pc[0] *= 7.0 / 24.0;
      pc[1] *= 1.875;
      pc[2] *= 49.0 / 12.0;
      pc[3] *= 3.5;
      pc[4] *= 1.0;
      break;

   case BC_MIXED:
      pc[0] *= 0.375;
      pc[1] *= 2.125;
      pc[2] *= 4.250;
      pc[3] *= 3.500;
      pc[4] *= 1.000;
      break;
   }

   // x-component
   std::function< real_t( const Point3D& ) > funcX = [pc]( const Point3D& x ) {
      real_t rad = sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
      WALBERLA_ASSERT( rad < 2.0 + 1e-12 );
      WALBERLA_ASSERT( rad > 1.0 - 1e-12 );
      real_t iRad  = 1.0 / rad;
      real_t value = ( pc[0] + ( pc[1] + ( pc[2] + ( pc[3] + pc[4] * iRad ) * iRad ) * iRad ) * iRad );

      return value * x[0] * iRad;
   };

   f[0].interpolate( funcX, level, All );
   f[0].multElementwise( {f[0], sph}, level, All );

   // y-component
   std::function< real_t( const Point3D& ) > funcY = [pc]( const Point3D& x ) {
      real_t rad = sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
      WALBERLA_ASSERT( rad < 2.0 + 1e-12 );
      WALBERLA_ASSERT( rad > 1.0 - 1e-12 );
      real_t iRad  = 1.0 / rad;
      real_t value = ( pc[0] + ( pc[1] + ( pc[2] + ( pc[3] + pc[4] * iRad ) * iRad ) * iRad ) * iRad );

      return value * x[1] * iRad;
   };

   f[1].interpolate( funcY, level, All );
   f[1].multElementwise( {f[1], sph}, level, All );

   // z-component
   std::function< real_t( const Point3D& ) > funcZ = [pc]( const Point3D& x ) {
      real_t rad = sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
      WALBERLA_ASSERT( rad < 2.0 + 1e-12 );
      WALBERLA_ASSERT( rad > 1.0 - 1e-12 );
      real_t iRad  = 1.0 / rad;
      real_t value = ( pc[0] + ( pc[1] + ( pc[2] + ( pc[3] + pc[4] * iRad ) * iRad ) * iRad ) * iRad );

      return value * x[2] * iRad;
   };

   f[2].interpolate( funcZ, level, All );
   f[2].multElementwise( {f[2], sph}, level, All );
}

// =======================
//  SetAnalyticalSolution
// =======================
template < typename vecFuncType >
void setAnalyticSolution( vecFuncType&                              u,
                          uint_t                                    level,
                          uint_t                                    degree,
                          int                                       order,
                          std::shared_ptr< SphericalHarmonicsTool > sphTool,
                          bcType                                    bmCase )
{
   // safety check
   if ( (uint_t) std::abs( order ) > degree )
   {
      WALBERLA_ABORT( "Spherical harmonics order must be smaller equal to degree!" );
   }

   // set polynomial coefficients for radial component
   real_t                  elfac = real_c( degree * ( degree + 1 ) );
   std::array< real_t, 5 > pCoeff;

   switch ( bmCase )
   {
   case BC_NOSLIP:
      pCoeff[0] = 1.00;
      pCoeff[1] = -3.00;
      pCoeff[2] = 3.25;
      pCoeff[3] = -1.50;
      pCoeff[4] = 0.25;
      break;

   case BC_FREESLIP:
      pCoeff[0] = 1.00;
      pCoeff[1] = -7.0 / 2.0;
      pCoeff[2] = 49.0 / 12.0;
      pCoeff[3] = -15.0 / 8.0;
      pCoeff[4] = 7.0 / 24.0;
      break;

   case BC_MIXED:
      pCoeff[0] = 1.00;
      pCoeff[1] = -7.0 / 2.0;
      pCoeff[2] = 17.0 / 4.0;
      pCoeff[3] = -17.0 / 8.0;
      pCoeff[4] = 3.0 / 8.0;
      break;
   }

   // define function for evaluating the different components
   uint_t                                    component;
   std::function< real_t( const Point3D& ) > vshFunc = [sphTool, degree, order, &component, pCoeff, elfac]( const Point3D& x ) {
      real_t value = real_c( 0 );
      real_t rad   = sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );

      // check validity of radial value
      WALBERLA_ASSERT( rad < 2.0 + 1e-12 );
      WALBERLA_ASSERT( rad > 1.0 - 1e-12 );

      real_t poly = ( ( ( pCoeff[4] * rad + pCoeff[3] ) * rad + pCoeff[2] ) * rad + pCoeff[1] ) * rad + pCoeff[0];

      real_t polyDeriv = ( ( 4.0 * pCoeff[4] * rad + 3.0 * pCoeff[3] ) * rad + 2.0 * pCoeff[2] ) * rad + pCoeff[1];

      // Add Y_(l,m)^0 contribution
      value = sphTool->evalVSH( degree, order, x[0], x[1], x[2], 0, component ) * poly * elfac / ( rad * rad );

      // Add Y_(l,m)^1 contribution
      value += sphTool->evalVSH( degree, order, x[0], x[1], x[2], 1, component ) * polyDeriv * sqrt( elfac ) / rad;

      return value;
   };

   // interpolate component functions
   for( component = 0; component < 3; ++component ) u[component].interpolate( vshFunc, level, All );

}

// ===================
//  runBenchmarkTests
// ===================
template < template < class > class feFuncType, typename massOpType, typename stokesOpType >
void runBenchmarkTests( std::shared_ptr< walberla::config::Config > cfg, std::shared_ptr< hyteg::PrimitiveStorage > storage )
{
   // parameter handling
   const walberla::Config::BlockHandle problemCfg = cfg->getBlock( "Problem" );

   const uint_t targetLevel = problemCfg.getParameter< uint_t >( "targetLevel" );

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
   mass.apply( forcing[0], rhs.uvw[0], maxLevel, All );
   mass.apply( forcing[1], rhs.uvw[1], maxLevel, All );
   mass.apply( forcing[2], rhs.uvw[2], maxLevel, All );

   // ---------------------
   //  Analytical Solution
   // ---------------------
   WALBERLA_LOG_PROGRESS_ON_ROOT( "Interpolating analytical solution to FEspace ..." );
   vf_t u_exact( "Analytical Solution (mapped to FE space)", storage, maxLevel, maxLevel );
   setAnalyticSolution( u_exact, maxLevel, degree, order, sphTool, bmarkBC );

   // -------------------
   //  Setup FE Solution
   // -------------------
   WALBERLA_LOG_PROGRESS_ON_ROOT( "Preparing FE function for discete solution ... " );
   feFuncType< real_t > feSol( "Discrete Solution", storage, minLevel, maxLevel, bcVelocity );

   // set everything to zero (including boundaries)
   feSol.uvw[0].setToZero( maxLevel );
   feSol.uvw[1].setToZero( maxLevel );
   feSol.uvw[2].setToZero( maxLevel );
   feSol.p.setToZero( maxLevel );

   // ---------------------
   //  Stuff for Free-Slip
   // ---------------------
   auto normalFunc = []( const Point3D& p, Point3D& n ) -> void {
      real_t norm = p.norm();
      real_t sign = ( norm > 0.95 * outerRadius ) ? 1.0 : -1.0;
      n           = sign / norm * p;
   };
   auto projectNormalsOp        = std::make_shared< P2ProjectNormalOperator >( storage, maxLevel, maxLevel, normalFunc );
   using StokesOperatorFreeSlip = StrongFreeSlipWrapper< stokesOpType, P2ProjectNormalOperator >;

   // -----------------------
   //  Setup Stokes Operator
   // -----------------------
   WALBERLA_LOG_PROGRESS_ON_ROOT( "Preparing Stokes operator ... " );
   std::shared_ptr< stokesOpType > stokesOp = std::make_shared< stokesOpType >( storage, maxLevel, maxLevel );
   // StrongFreeSlipWrapper< stokesOpType, P2ProjectNormalOperator > stokesOpFS( stokesOp, projectNormalsOp, bmarkFsDoFs );
   StokesOperatorFreeSlip stokesOpFS( stokesOp, projectNormalsOp, bmarkFsDoFs );
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
         WALBERLA_LOG_PROGRESS_ON_ROOT( "Iterative with PETScBPSS ... " );
         WALBERLA_ABORT( "Need to define the block preconditioner for the P2P1ElementwiseBlendingStokesOperator for this one!" );
         // PETScBlockPreconditionedStokesSolver< stokesOpType > stokesSolver( storage, maxLevel, 1e-08, 500, 1 );
         // stokesSolver.solve( stokesOpFS, feSol, rhs, maxLevel );
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
      // PETScMinResSolver< StokesOperatorFreeSlip > stokesSolver( storage, maxLevel, 1e-08, 200 );
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
   error.assign( {1.0, -1.0}, {u_exact, feSol.uvw}, maxLevel );

   if ( problemCfg.getParameter< bool >( "compute_L2_error" ) )
   {
      real_t L2error = 0.0;
      typename vf_t::template VectorComponentType aux( "aux", storage, maxLevel, maxLevel );

      for( uint_t idx = 0; idx < 3; ++idx )
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

      typename vf_t::template VectorComponentType aux( "aux", storage, maxLevel, maxLevel+1 );
      P2toP2QuadraticProlongation embedder;

      // TEST
      // feSol.prolongate( maxLevel, All ); there is no P2Function::prolongate !!!

      for( uint_t idx = 0; idx < 3; ++idx )
      {
        aux.assign( {1.0}, {feSol.uvw[idx]}, maxLevel );
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
      for( uint_t idx = 0; idx < 3; ++idx )
      {
        vtkOutput.add( forcing[idx] );
        vtkOutput.add( rhs.uvw[idx] );
        vtkOutput.add( u_exact[idx] );
        vtkOutput.add( feSol.uvw[idx] );
        vtkOutput.add( error[idx] );
      }
      vtkOutput.write( maxLevel, 0 );
   }
}
} // namespace terraneo

// Hunting the NaN
#include <cfenv>

int main( int argc, char* argv[] )
{
#ifndef __APPLE__
   feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );
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
   const uint_t nRad = problemCfg.getParameter< uint_t >( "nRad" );
   const uint_t nTan = problemCfg.getParameter< uint_t >( "nTan" );

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

   std::shared_ptr< walberla::WcTimingTree >  timingTree( new walberla::WcTimingTree() );
   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, timingTree );

   if ( problemCfg.getParameter< bool >( "reportPrimitives" ) )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "" << setupStorage );
   }

   // =================
   //  Benchmark Stuff
   // =================

   // delegate actual work to templated function to allow different FE spaces
   std::string feSpace = problemCfg.getParameter< std::string >( "feSpace" );

   if ( feSpace == "P2P1TaylorHood" )
   {
      terraneo::runBenchmarkTests< P2P1TaylorHoodFunction,
                                   P2ElementwiseBlendingMassOperator,
                                   P2P1ElementwiseBlendingStokesOperator >( cfg, storage );
   }

   return EXIT_SUCCESS;
}
