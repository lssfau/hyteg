/*
 * Copyright (c) 2017-2019 Benjamin Mann
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
#define FP_FAST_FMA
#define FP_FAST_FMAF
#define FP_FAST_FMAL

#include <core/Environment.h>
#include <core/config/Create.h>
#include <core/timing/Timer.h>
#include <hyteg/composites/P2P1BlendingTaylorHoodStokesOperator.hpp>
#include <hyteg/composites/P2P1SurrogateTaylorHoodStokesOperator.hpp>
#include <hyteg/composites/P2P1TaylorHoodStokesOperator.hpp>
#include <hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp>
#include <hyteg/mesh/MeshInfo.hpp>
#include <hyteg/p1functionspace/P1VariableOperator.hpp>
#include <hyteg/p2functionspace/P2ConstantOperator.hpp>
#include <hyteg/p2functionspace/P2VariableOperator.hpp>

#include "core/Format.hpp"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/CircularMap.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"

using walberla::real_t;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

using c_function = std::function< real_t( const hyteg::Point3D& ) >;

enum StencilType
{
   NONE     = -1,
   CONST    = 0,
   VARIABLE = 1,
   LSQP     = 2,
   ELWISE   = 3
};

// cartesian to polar coordinates
c_function radius = []( const hyteg::Point3D& x ) { return std::hypot( x[0], x[1] ); };
c_function angle  = []( const hyteg::Point3D& x ) { return std::atan2( x[1], x[0] ); };

template < StencilType ST, class StokesOperator >
real_t solve( std::shared_ptr< StokesOperator >   L,
              std::shared_ptr< PrimitiveStorage > storage,
              const uint_t                        minLevel,
              const uint_t                        maxLevel,
              const uint_t                        max_outer_iter,
              const uint_t                        max_cg_iter,
              const real_t                        mg_tolerance,
              const real_t                        coarse_tolerance,
              const bool                          vtk,
              c_function&                         u_exact,
              c_function&                         v_exact,
              c_function&                         p_exact,
              c_function&                         u_boundary,
              c_function&                         v_boundary,
              c_function&                         rhs_x,
              c_function&                         rhs_y,
              c_function&                         T_field )
{
   // functions and operators
   P2BlendingMassOperator M2( storage, maxLevel, maxLevel );
   P1BlendingMassOperator M1( storage, maxLevel, maxLevel );

   P2P1TaylorHoodFunction< real_t > r( "r", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > f( "f", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > up( "up", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > Lu( "Lu", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > up_exact( "up_exact", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > err( "err", storage, minLevel, maxLevel );
   P2Function< real_t >             T( "T", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > tmp( "", storage, minLevel, maxLevel );

   Lu.uvw()[0].setToZero( maxLevel );
   Lu.uvw()[1].setToZero( maxLevel );
   Lu.p().setToZero( maxLevel );

   // init boundary
   up.uvw().interpolate( { u_boundary, v_boundary }, maxLevel, DirichletBoundary );

   // init rhs
   tmp.uvw().interpolate( { rhs_x, rhs_y }, maxLevel, All );
   M2.apply( tmp.uvw()[0], f.uvw()[0], maxLevel, All );
   M2.apply( tmp.uvw()[1], f.uvw()[1], maxLevel, All );
   f.p().setToZero( maxLevel );

   // init exact solution
   T.interpolate( T_field, maxLevel, All );
   up_exact.uvw().interpolate( { u_exact, v_exact }, maxLevel, All );
   up_exact.p().interpolate( p_exact, maxLevel, All );
   vertexdof::projectMean( up_exact.p(), maxLevel );

   communication::syncP2FunctionBetweenPrimitives( up_exact.uvw()[0], maxLevel );
   communication::syncP2FunctionBetweenPrimitives( up_exact.uvw()[1], maxLevel );
   communication::syncFunctionBetweenPrimitives( up_exact.p(), maxLevel );

   // solver
   const real_t uzawaOmega = 0.37;

   auto coarseGridSolver =
       solvertemplates::stokesMinResSolver< StokesOperator >( storage, minLevel, coarse_tolerance, max_cg_iter );
   auto fineGridSolver = solvertemplates::stokesMinResSolver< StokesOperator >( storage, maxLevel, mg_tolerance, max_cg_iter );

   auto restriction  = std::make_shared< P2P1StokesToP2P1StokesRestriction >( true );
   auto prolongation = std::make_shared< P2P1StokesToP2P1StokesProlongation >();
   auto smoother     = std::make_shared< GaussSeidelSmoother< typename StokesOperator::VelocityOperator_T > >();
   auto preconditioner =
       std::make_shared< StokesVelocityBlockBlockDiagonalPreconditioner< StokesOperator > >( storage, smoother );
   auto uzawaSmoother = std::make_shared< UzawaSmoother< StokesOperator > >(
       storage, preconditioner, minLevel, maxLevel, uzawaOmega, Inner | NeumannBoundary, 4 );
   auto gmgSolver = std::make_shared< GeometricMultigridSolver< StokesOperator > >(
       storage, uzawaSmoother, coarseGridSolver, restriction, prolongation, minLevel, maxLevel, 10, 10, 2 );

   // solve iteratively
   uint_t       iter       = 0;
   real_t       res        = 0, res_old, discr_l2_err_uv, discr_l2_err_u, discr_l2_err_v, discr_l2_err_p;
   real_t       vCycleTime = 0, solveTime = 0;
   real_t       averageConvergenceRate = 0;
   const uint_t convergenceStartIter   = 3;

   WALBERLA_LOG_INFO_ON_ROOT( "Starting V cycles" );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%6s|%10s|%10s|%10s|%10s|%10s|%10s|%10s",
                                                "iter",
                                                "res",
                                                "conv",
                                                "L2-error |u|",
                                                "L2-error u",
                                                "L2-error v",
                                                "L2-error p",
                                                "Cycle-Time" ) );

   while ( 1 )
   {
      // compute error
      vertexdof::projectMean( up.p(), maxLevel );
      err.assign( { 1.0, -1.0 }, { up, up_exact }, maxLevel );
      M2.apply( err.uvw()[0], tmp.uvw()[0], maxLevel, Inner | NeumannBoundary, Replace );
      M2.apply( err.uvw()[1], tmp.uvw()[1], maxLevel, Inner | NeumannBoundary, Replace );
      M1.apply( err.p(), tmp.p(), maxLevel, Inner | NeumannBoundary, Replace );
      discr_l2_err_u  = std::sqrt( err.uvw()[0].dotGlobal( tmp.uvw()[0], maxLevel, Inner | NeumannBoundary ) );
      discr_l2_err_v  = std::sqrt( err.uvw()[1].dotGlobal( tmp.uvw()[1], maxLevel, Inner | NeumannBoundary ) );
      discr_l2_err_p  = std::sqrt( err.p().dotGlobal( tmp.p(), maxLevel, Inner | NeumannBoundary ) );
      discr_l2_err_uv = std::sqrt( discr_l2_err_u * discr_l2_err_u + discr_l2_err_v * discr_l2_err_v );

      // compute residual
      res_old = res;
      L->apply( up, Lu, maxLevel, hyteg::Inner | NeumannBoundary );
      r.assign( { 1.0, -1.0 }, { f, Lu }, maxLevel, hyteg::Inner | NeumannBoundary );
      M2.apply( r.uvw()[0], tmp.uvw()[0], maxLevel, hyteg::All, Replace );
      M2.apply( r.uvw()[1], tmp.uvw()[1], maxLevel, hyteg::All, Replace );
      M1.apply( r.p(), tmp.p(), maxLevel, hyteg::All, Replace );
      res = std::sqrt( r.dotGlobal( tmp, maxLevel, All ) );

      // compute convergence rate
      real_t convRate = res / res_old;

      if ( iter >= convergenceStartIter )
      {
         averageConvergenceRate += convRate;
      }

      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%6d|%10.3e|%10.3e|%10.3e|%10.3e|%10.3e|%10.3e|%10.3e",
                                                   iter,
                                                   res,
                                                   convRate,
                                                   discr_l2_err_uv,
                                                   discr_l2_err_u,
                                                   discr_l2_err_v,
                                                   discr_l2_err_p,
                                                   vCycleTime ) );

      // stopping criterion
      if ( ++iter > max_outer_iter || ( iter > convergenceStartIter && convRate > 0.999 ) )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Ending multigrid without reaching desired tolerance!" );
         break;
      }

      if ( res < mg_tolerance )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Multigrid converged!" );
         break;
      }

      // solve
      auto start = walberla::timing::getWcTime();
      gmgSolver->solve( *L, up, f, maxLevel );
      auto end   = walberla::timing::getWcTime();
      vCycleTime = end - start;
      solveTime += vCycleTime;
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Solve time " << std::defaultfloat << solveTime );
   WALBERLA_LOG_INFO_ON_ROOT( "Avg. convergence rate: " << std::scientific
                                                        << averageConvergenceRate / real_c( iter - convergenceStartIter ) );
   WALBERLA_LOG_INFO_ON_ROOT( "L^2 error U: " << std::scientific << discr_l2_err_uv );
   WALBERLA_LOG_INFO_ON_ROOT( "L^2 error u: " << std::scientific << discr_l2_err_u );
   WALBERLA_LOG_INFO_ON_ROOT( "L^2 error v: " << std::scientific << discr_l2_err_v );
   WALBERLA_LOG_INFO_ON_ROOT( "L^2 error p: " << std::scientific << discr_l2_err_p );

   if ( vtk )
   {
      std::string name = "Stokes_Surrogates_";

      switch ( ST )
      {
      case StencilType::CONST:
         name += "const";
         break;

      case StencilType::VARIABLE:
         name += "var";
         break;

      case StencilType::LSQP:
         name += "surrogate";
         break;
      }

      hyteg::VTKOutput vtkOutput( "../../output", name, storage );
      vtkOutput.add( up );
      vtkOutput.add( up_exact );
      vtkOutput.add( err );
      vtkOutput.add( T );
      vtkOutput.write( maxLevel, 0 );
   }

   return solveTime;
}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   walberla::shared_ptr< walberla::config::Config > cfg;

   if ( argc == 1 )
   {
      walberla::shared_ptr< walberla::config::Config > cfg_( new walberla::config::Config );
      cfg_->readParameterFile( "../../data/param/Stokes2d_Surrogates.prm" );
      cfg = cfg_;
   }
   else
   {
      cfg = walberla::config::create( argc, argv );
   }

   // read parameter file
   WALBERLA_LOG_INFO_ON_ROOT( "config = " << *cfg );
   walberla::Config::BlockHandle parameters = cfg->getOneBlock( "Parameters" );

   const uint_t minLevel = parameters.getParameter< uint_t >( "level_h_coarse" );
   const uint_t maxLevel = parameters.getParameter< uint_t >( "level_h_fine" );

   const int         opTypeInt = parameters.getParameter< int >( "operatorType" );
   const StencilType opType    = StencilType( opTypeInt );

   const bool   blending   = parameters.getParameter< bool >( "blending" );
   const bool   annulus    = parameters.getParameter< bool >( "annulus" );
   uint_t       nX         = parameters.getParameter< uint_t >( "nX" );
   uint_t       nY         = parameters.getParameter< uint_t >( "nY" );
   const uint_t macroLevel = parameters.getParameter< uint_t >( "macro_refinement" );
   nX                      = nX << macroLevel;
   nY                      = nY << macroLevel;

   const uint_t polyDegree            = parameters.getParameter< uint_t >( "polyDegree" );
   const uint_t maxInterpolationLevel = parameters.getParameter< uint_t >( "interpolationLevel" );
   const uint_t interpolationLevel    = std::min( maxLevel, maxInterpolationLevel );

   const uint_t max_outer_iter   = parameters.getParameter< uint_t >( "max_outer_iter" );
   const uint_t max_cg_iter      = parameters.getParameter< uint_t >( "max_cg_iter" );
   const real_t mg_tolerance     = parameters.getParameter< real_t >( "mg_tolerance" );
   const real_t coarse_tolerance = parameters.getParameter< real_t >( "coarse_tolerance" );

   const bool vtk = parameters.getParameter< bool >( "vtkOutput" );

   if ( annulus )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Geometry: Annulus" );

      if ( blending )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Geometry blending enabled" );
      }
      else
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Geometry blending disabled" );
      }
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Geometry: Rectangle" );
   }

   // define functions and domain

   /// case rectangle
   c_function u_exact    = []( const hyteg::Point3D& x ) { return sin( 2 * pi * x[0] ) * cos( pi * x[1] ); };
   c_function v_exact    = []( const hyteg::Point3D& x ) { return -2 * cos( 2 * pi * x[0] ) * sin( pi * x[1] ); };
   c_function p_exact    = []( const hyteg::Point3D& x ) { return 2.5 * pi * cos( 2 * pi * x[0] ) * cos( pi * x[1] ); };
   c_function u_boundary = u_exact;
   c_function v_boundary = v_exact;
   c_function T_field    = []( const hyteg::Point3D& x ) { return -12.5 * pi * pi * cos( 2 * pi * x[0] ) * sin( pi * x[1] ); };
   c_function rhs_x      = []( const hyteg::Point3D& ) { return 0; };
   c_function rhs_y      = [&]( const hyteg::Point3D& x ) { return T_field( x ); };

   MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( { 0.0, 0.0 } ), Point2D( { 1.0, 1.0 } ), MeshInfo::CRISS, nX, nY );

   /// case annulus
   Point3D                                               circleCenter{ { 0, 0, 0 } };
   const real_t                                          rMin    = pi;
   const real_t                                          rMax    = 2 * pi;
   const real_t                                          middle  = ( rMax + rMin ) / 2.0;
   const real_t                                          k       = 8.0;
   std::function< real_t( const real_t, const real_t ) > u_polar = [&]( const real_t r, const real_t phi ) {
      const auto sinPhi  = std::sin( phi );
      const auto cosPhi  = std::cos( phi );
      const auto sinKPhi = std::sin( k * phi );
      const auto cosKPhi = std::cos( k * phi );
      const auto sinR    = std::sin( r );
      const auto cosR    = std::cos( r );
      return cosPhi * ( k / r ) * sinR * cosKPhi + sinPhi * cosR * sinKPhi;
   };
   std::function< real_t( const real_t, const real_t ) > v_polar = [&]( const real_t r, const real_t phi ) {
      const auto sinPhi  = std::sin( phi );
      const auto cosPhi  = std::cos( phi );
      const auto sinKPhi = std::sin( k * phi );
      const auto cosKPhi = std::cos( k * phi );
      const auto sinR    = std::sin( r );
      const auto cosR    = std::cos( r );
      return sinPhi * ( k / r ) * sinR * cosKPhi - cosPhi * cosR * sinKPhi;
   };

   if ( annulus )
   {
      u_exact = [&]( const hyteg::Point3D& x ) { return u_polar( radius( x ), angle( x ) ); };
      v_exact = [&]( const hyteg::Point3D& x ) { return v_polar( radius( x ), angle( x ) ); };
      p_exact = [&]( const hyteg::Point3D& x ) {
         const auto r       = radius( x );
         const auto phi     = angle( x );
         const auto cosKPhi = std::cos( k * phi );
         const auto sinR    = std::sin( r );
         const auto cosR    = std::cos( r );
         return ( 2.0 / ( r * r ) ) * k * sinR * cosKPhi - ( 1.0 / k ) * sinR * cosKPhi -
                ( ( k * k + r * r + 1.0 ) / ( k * r ) ) * cosR * cosKPhi;
      };
      u_boundary = [&]( const hyteg::Point3D& x ) {
         const real_t r = ( radius( x ) < middle ) ? rMin : rMax;
         return u_polar( r, angle( x ) );
      };
      v_boundary = [&]( const hyteg::Point3D& x ) {
         const real_t r = ( radius( x ) < middle ) ? rMin : rMax;
         return v_polar( r, angle( x ) );
      };
      T_field = [&]( const hyteg::Point3D& x ) {
         const auto r       = radius( x );
         const auto phi     = angle( x );
         const auto cosKPhi = std::cos( k * phi );
         const auto sinR    = std::sin( r );
         const auto cosR    = std::cos( r );
         const auto velU    = ( k / r ) * sinR * cosKPhi;
         const auto r_2     = r * r;
         const auto k_2     = k * k;
         return ( k_2 + r_2 - 4.0 + ( r_2 / k_2 ) * ( k_2 + r_2 + 1.0 ) ) * velU / r_2 +
                ( 2.0 * k_2 - 2.0 * r_2 + 1.0 ) / k * cosR * cosKPhi / r_2;
      };
      rhs_x = [&]( const hyteg::Point3D& x ) {
         const auto phi    = angle( x );
         const auto cosPhi = std::cos( phi );
         return cosPhi * T_field( x );
      };
      rhs_y = [&]( const hyteg::Point3D& x ) {
         const auto phi    = angle( x );
         const auto sinPhi = std::sin( phi );
         return sinPhi * T_field( x );
      };

      meshInfo = MeshInfo::meshAnnulus( rMin, rMax, MeshInfo::CRISS, nX, nY );
   }

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   if ( annulus && blending )
      AnnulusMap::setMap( setupStorage );

   hyteg::loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // WALBERLA_LOG_INFO_ON_ROOT(storage->getGlobalInfo());
   WALBERLA_LOG_INFO_ON_ROOT( "Refinement levels: " << minLevel << "->" << maxLevel );

   real_t setupTime = 0, startSetup = 0, endSetup = 0, solveTime;

   auto L1 = std::make_shared< P2P1TaylorHoodStokesOperator >( storage, minLevel, maxLevel );
   auto L2 = std::make_shared< P2P1BlendingTaylorHoodStokesOperator >( storage, minLevel, maxLevel );
   auto L3 = std::make_shared< P2P1ElementwiseBlendingStokesOperator >( storage, minLevel, maxLevel );
   auto L4 = std::make_shared< P2P1SurrogateTaylorHoodStokesOperator >( storage, minLevel, maxLevel, interpolationLevel );

   switch ( opType )
   {
   case CONST:
      WALBERLA_LOG_INFO_ON_ROOT( "Operatortype: Constant Stencil" );
      solveTime = solve< CONST, P2P1TaylorHoodStokesOperator >( L1,
                                                                storage,
                                                                minLevel,
                                                                maxLevel,
                                                                max_outer_iter,
                                                                max_cg_iter,
                                                                mg_tolerance,
                                                                coarse_tolerance,
                                                                vtk,
                                                                u_exact,
                                                                v_exact,
                                                                p_exact,
                                                                u_boundary,
                                                                v_boundary,
                                                                rhs_x,
                                                                rhs_y,
                                                                T_field );
      break;

   case VARIABLE:
      WALBERLA_LOG_INFO_ON_ROOT( "Operatortype: Variable Stencil" );
      solveTime = solve< VARIABLE, P2P1BlendingTaylorHoodStokesOperator >( L2,
                                                                           storage,
                                                                           minLevel,
                                                                           maxLevel,
                                                                           max_outer_iter,
                                                                           max_cg_iter,
                                                                           mg_tolerance,
                                                                           coarse_tolerance,
                                                                           vtk,
                                                                           u_exact,
                                                                           v_exact,
                                                                           p_exact,
                                                                           u_boundary,
                                                                           v_boundary,
                                                                           rhs_x,
                                                                           rhs_y,
                                                                           T_field );
      break;

   case LSQP:
      WALBERLA_LOG_INFO_ON_ROOT( "Operatortype: Surrogate Polynomial Stencil" );
      WALBERLA_LOG_INFO_ON_ROOT( "Interpolation level: " << interpolationLevel );
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Apply LSQ-fit with q = %d", polyDegree ) );
      startSetup = walberla::timing::getWcTime();
      L4->interpolateStencils( polyDegree );
      endSetup = walberla::timing::getWcTime();
      L4->useDegree( polyDegree );
      setupTime = ( endSetup - startSetup );
      solveTime = solve< LSQP, P2P1SurrogateTaylorHoodStokesOperator >( L4,
                                                                        storage,
                                                                        minLevel,
                                                                        maxLevel,
                                                                        max_outer_iter,
                                                                        max_cg_iter,
                                                                        mg_tolerance,
                                                                        coarse_tolerance,
                                                                        vtk,
                                                                        u_exact,
                                                                        v_exact,
                                                                        p_exact,
                                                                        u_boundary,
                                                                        v_boundary,
                                                                        rhs_x,
                                                                        rhs_y,
                                                                        T_field );
      break;

   default:
      WALBERLA_ABORT( "The desired Operator Type is not supported!" );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Setup time: " << std::defaultfloat << setupTime );
   WALBERLA_LOG_INFO_ON_ROOT( "Time to solution: " << std::defaultfloat << setupTime + solveTime );
   WALBERLA_LOG_INFO_ON_ROOT( "=======================================================" );
}
