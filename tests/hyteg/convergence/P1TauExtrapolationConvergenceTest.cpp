/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
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

#include "core/Environment.h"
#include "core/logging/Logging.h"

#include "hyteg/gridtransferoperators/P1toP1InjectionRestriction.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/gridtransferoperators/P2toP1Conversion.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/FAS.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"

#include "constant_stencil_operator/P1ConstantOperator.hpp"
#include "constant_stencil_operator/P2ConstantOperator.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   const uint_t      minLevel                  = 2;
   const uint_t      maxLevel                  = 5;
   const std::string meshFile                  = hyteg::prependHyTeGMeshDir( "2D/quad_8el.msh" );
   const real_t      coarseGridSolverTolerance = real_c( 1e-16 );
   const uint_t      maxCoarseGridSolverIter   = 10000;
   const uint_t      numVCycles                = 10;

   auto meshInfo = MeshInfo::fromGmshFile( meshFile );
   auto setupStorage =
       std::make_shared< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage );

   hyteg::P1Function< real_t > r( "r", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > f( "f", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > u( "u", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > Au( "Au", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > AIu( "AIu", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > IAu( "IAu", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > u_exact( "u_exact", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > err( "err", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > npoints_helper( "npoints_helper", storage, minLevel, maxLevel );

   hyteg::P2Function< real_t > quadraticTmp( "tmp_P2", storage, minLevel, maxLevel - 1 );
   hyteg::P2Function< real_t > quadraticRhs( "f_P2", storage, minLevel, maxLevel - 1 );

   hyteg::P1ConstantMassOperator    M( storage, minLevel, maxLevel );
   hyteg::P2ConstantMassOperator    M_P2( storage, minLevel, maxLevel );
   hyteg::P1ConstantLaplaceOperator L( storage, minLevel, maxLevel );

#if 1
   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& x ) {
      return ( 1.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
   };
   std::function< real_t( const hyteg::Point3D& ) > rhs = []( const hyteg::Point3D& x ) {
      return ( 3.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
   };
#else
   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& x ) { return sin( x[0] ) * sinh( x[1] ); };
   std::function< real_t( const hyteg::Point3D& ) > rhs   = []( const hyteg::Point3D& ) { return real_c( 0 ); };
#endif

   std::function< real_t( const hyteg::Point3D& ) > ones = []( const hyteg::Point3D& ) { return 1.0; };

   u.interpolate( exact, maxLevel, hyteg::DirichletBoundary );
   u_exact.interpolate( exact, maxLevel );
   u_exact.interpolate( exact, maxLevel - 1 );
   npoints_helper.interpolate( rhs, maxLevel );
   M.apply( npoints_helper, f, maxLevel, hyteg::All );

   auto smoother         = std::make_shared< hyteg::GaussSeidelSmoother< hyteg::P1ConstantLaplaceOperator > >();
   auto coarseGridSolver = std::make_shared< hyteg::CGSolver< hyteg::P1ConstantLaplaceOperator > >(
       storage, minLevel, minLevel, maxCoarseGridSolverIter, coarseGridSolverTolerance );
   auto restrictionOperator         = std::make_shared< hyteg::P1toP1LinearRestriction<> >();
   auto solutionRestrictionOperator = std::make_shared< hyteg::P1toP1InjectionRestriction >();
   auto prolongationOperator        = std::make_shared< hyteg::P1toP1LinearProlongation<> >();

   auto gmgSolver = hyteg::FASSolver< hyteg::P1ConstantLaplaceOperator >( storage,
                                                                          smoother,
                                                                          coarseGridSolver,
                                                                          restrictionOperator,
                                                                          solutionRestrictionOperator,
                                                                          prolongationOperator,
                                                                          minLevel,
                                                                          maxLevel,
                                                                          3,
                                                                          3 );

   auto gmgSolverTau = hyteg::GeometricMultigridSolver< hyteg::P1ConstantLaplaceOperator >(
       storage, smoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel - 1, 1, 1 );

   npoints_helper.interpolate( ones, maxLevel );
   const real_t npoints = npoints_helper.dotGlobal( npoints_helper, maxLevel );

   npoints_helper.interpolate( ones, maxLevel - 1 );
   const real_t npoints_tau = npoints_helper.dotGlobal( npoints_helper, maxLevel - 1 );

   // init residual once
   L.apply( u, Au, maxLevel, hyteg::Inner );
   r.assign( { 1.0, -1.0 }, { f, Au }, maxLevel, hyteg::Inner );

   real_t discr_l2_err;
   real_t discr_l2_res           = std::sqrt( r.dotGlobal( r, maxLevel, DoFType::Inner ) / npoints );
   real_t discr_l2_res_last_step = discr_l2_res;

   for ( uint_t vCycleCount = 0; vCycleCount < numVCycles; vCycleCount++ )
   {
      gmgSolver.solve( L, u, f, maxLevel );

      err.assign( { 1.0, -1.0 }, { u, u_exact }, maxLevel );

      L.apply( u, Au, maxLevel, hyteg::Inner );
      r.assign( { 1.0, -1.0 }, { f, Au }, maxLevel, hyteg::Inner );

      discr_l2_err           = std::sqrt( err.dotGlobal( err, maxLevel, DoFType::All ) / npoints );
      discr_l2_res_last_step = discr_l2_res;
      discr_l2_res           = std::sqrt( r.dotGlobal( r, maxLevel, DoFType::Inner ) / npoints );

      const real_t convRate = discr_l2_res / discr_l2_res_last_step;

      WALBERLA_LOG_INFO_ON_ROOT( "After " << vCycleCount
                                          << " V-cycles: "
                                             "discrete L2 error = "
                                          << discr_l2_err << ", discrete L2 residual = " << discr_l2_res
                                          << ", conv rate = " << convRate )
   }

   WALBERLA_CHECK_LESS( discr_l2_err, 2.9e-06 );

   // perform tau-extrapolation

   // I * A * u
   L.apply( u, IAu, maxLevel, hyteg::Inner );
   restrictionOperator->restrict( IAu, maxLevel, hyteg::All );

   // A * I * u (injection)
   solutionRestrictionOperator->restrict( u, maxLevel, hyteg::All );
   L.apply( u, AIu, maxLevel - 1, hyteg::All );

   // build RHS (quadratic, then restrict with weighting)
   quadraticTmp.interpolate( rhs, maxLevel - 1, hyteg::All );
   M_P2.apply( quadraticTmp, quadraticRhs, maxLevel - 1, hyteg::All );
   // f.assign( quadraticRhs, maxLevel, hyteg::All );
   P2toP1Conversion( quadraticRhs, f, maxLevel, hyteg::All );

   restrictionOperator->restrict( f, maxLevel, All );
   f.assign( { 1.0, real_c( -4.0 / 3.0 ), real_c( 4.0 / 3.0 ) }, { f, IAu, AIu }, maxLevel - 1 );

   const uint_t tauMaxLevel = maxLevel - 1;

   for ( uint_t vCycleCount = 0; vCycleCount < numVCycles; vCycleCount++ )
   {
      gmgSolverTau.solve( L, u, f, tauMaxLevel );

      err.assign( { 1.0, -1.0 }, { u, u_exact }, tauMaxLevel );

      L.apply( u, Au, tauMaxLevel, hyteg::Inner );
      r.assign( { 1.0, -1.0 }, { f, Au }, tauMaxLevel, hyteg::Inner );

      discr_l2_err           = std::sqrt( err.dotGlobal( err, tauMaxLevel, DoFType::All ) / npoints_tau );
      discr_l2_res_last_step = discr_l2_res;
      discr_l2_res           = std::sqrt( r.dotGlobal( r, tauMaxLevel, DoFType::Inner ) / npoints_tau );

      const real_t convRate = discr_l2_res / discr_l2_res_last_step;

      WALBERLA_LOG_INFO_ON_ROOT( "After " << vCycleCount
                                          << " V-cycles: "
                                             "discrete L2 error = "
                                          << discr_l2_err << ", discrete L2 residual = " << discr_l2_res
                                          << ", conv rate = " << convRate )
   }

   auto dp = std::is_same< real_t, double >();
   WALBERLA_CHECK_LESS( discr_l2_err, dp ? 3.3e-09 : 2e-7 );

   return 0;
}
