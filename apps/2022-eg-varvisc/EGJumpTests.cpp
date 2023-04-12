/*
* Copyright (c) 2023, Fabian Böhm.
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

#include "core/DataTypes.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/MeshQuality.hpp"
#include "hyteg/composites/P1DGEP0StokesOperator.hpp"
#include "hyteg/egfunctionspace/EGConvTestUtils.hpp"
#include "hyteg/egfunctionspace/EGOperators.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.cpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/MinresSolver.hpp"

using walberla::real_t;
using walberla::uint_t;

using hyteg::Point3D;
using hyteg::dg::eg::EGMassOperator;
using hyteg::dg::eg::EGP0EpsilonStokesOperator;
using hyteg::dg::eg::EGP0StokesOperator;
using hyteg::dg::eg::EGSIPGLaplaceOperator;

namespace hyteg {
namespace dg {
namespace eg {

void StraightJump3D( const uint_t minLevel, const uint_t maxLevel );

void IncreasingJump3D( const uint_t minLevel, const uint_t maxLevel );
} // namespace eg
} // namespace dg
} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::PETScManager petscManager( &argc, &argv );

   /* commandline arguments for petsc solver:
   -ksp_monitor -ksp_rtol 1e-7 -ksp_type minres  -pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type diag  -fieldsplit_0_ksp_type cg -fieldsplit_1_ksp_type cg -pc_fieldsplit_detect_saddle_point -fieldsplit_1_ksp_constant_null_space
   */

   uint_t minLevel = 3;

   if constexpr ( false )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Testing Straight jumps Epsilon 3D ###" )
      hyteg::dg::eg::StraightJump3D( minLevel, 5 );
   }

   if constexpr ( true )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Testing Increasing Jump Epsilon 3D ###" )
      hyteg::dg::eg::IncreasingJump3D( minLevel, 4 );
   }
   return 0;
}

namespace hyteg {
namespace dg {
namespace eg {

void IncreasingJump3D( const uint_t minLevel, const uint_t maxLevel )
{
   {
      real_t alpha = 1;
      auto   mu    = [&alpha]( const hyteg::Point3D& p ) {
         const real_t x = p[0];
         const real_t y = p[1];
         const real_t z = p[2];
         return std::tanh( alpha * ( x + y + z ) - 1.5 ) + 2;
      };

      // cube_6el, inhom. solution, EGP0
      if ( true )
      {
         auto resNormsP2P1 = { 1e-5, 1e-7, 1e-7, 1e-5, 1e-6, 1e-7, 1e-5, 1e-6 };

         auto meshInfo = hyteg::MeshInfo::meshCuboid( Point3D( { -1, -1, -1 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );

         hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                                    walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
         setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
         auto storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, 1 );

         // cube, hom. solution
         WALBERLA_LOG_INFO_ON_ROOT( "### EGP0, cube_6el, IncreasingJump3D, alpha = 1  ###" );

         StokesConvergenceOrderTest< EGP0EpsilonOperatorStokesNitscheBC >(
             "EGP0EpsilonOperatorStokesNitscheBC_IncreasingJump3D_alpha5",
             std::make_tuple(
                 []( const hyteg::Point3D& p ) {
                    const real_t x = p[0];
                    const real_t y = p[1];
                    const real_t z = p[2];
                    return y + std::sin( M_PI * ( x + y + z ) );
                 },
                 []( const hyteg::Point3D& p ) {
                    const real_t x = p[0];
                    const real_t y = p[1];
                    const real_t z = p[2];

                    return z - 2 * std::sin( M_PI * ( x + y + z ) );
                 },
                 []( const hyteg::Point3D& p ) {
                    const real_t x = p[0];
                    const real_t y = p[1];
                    const real_t z = p[2];

                    return x + std::sin( M_PI * ( x + y + z ) );
                 },
                 []( const hyteg::Point3D& p ) {
                    const real_t x = p[0];
                    const real_t y = p[1];
                    const real_t z = p[2];
                    return 2 * x - 4 * y + 10 * z;
                 } ),
             std::make_tuple(
                 [&alpha]( const hyteg::Point3D& p ) {
                    const real_t x  = p[0];
                    const real_t y  = p[1];
                    const real_t z  = p[2];
                    const real_t x0 = x + y + z;
                    const real_t x1 = M_PI * x0;
                    const real_t x2 = std::tanh( alpha * x0 - 1.5 );
                    const real_t x3 = M_PI * std::cos( x1 );
                    const real_t x4 = alpha * ( 1 - std::pow( x2, 2 ) );
                    const real_t x5 = 2 * x4;
                    return -2.0 * x3 * x4 - x5 * ( 0.5 - 0.5 * x3 ) - x5 * ( 1.0 * x3 + 0.5 ) +
                           3.0 * std::pow( M_PI, 2 ) * ( x2 + 2 ) * std::sin( x1 ) + 2;
                 },
                 [&alpha]( const hyteg::Point3D& p ) {
                    const real_t x  = p[0];
                    const real_t y  = p[1];
                    const real_t z  = p[2];
                    const real_t x0 = x + y + z;
                    const real_t x1 = M_PI * x0;
                    const real_t x2 = std::tanh( alpha * x0 - 1.5 );
                    const real_t x3 = std::cos( x1 );
                    const real_t x4 = 1 - std::pow( x2, 2 );
                    return 4.0 * M_PI * alpha * x3 * x4 - 4 * alpha * x4 * ( -0.5 * M_PI * x3 + 0.5 ) -
                           6.0 * std::pow( M_PI, 2 ) * ( x2 + 2 ) * std::sin( x1 ) - 4;
                 },
                 [&alpha]( const hyteg::Point3D& p ) {
                    const real_t x  = p[0];
                    const real_t y  = p[1];
                    const real_t z  = p[2];
                    const real_t x0 = x + y + z;
                    const real_t x1 = M_PI * x0;
                    const real_t x2 = std::tanh( alpha * x0 - 1.5 );
                    const real_t x3 = M_PI * std::cos( x1 );
                    const real_t x4 = alpha * ( 1 - std::pow( x2, 2 ) );
                    const real_t x5 = 2 * x4;
                    return -2.0 * x3 * x4 - x5 * ( 0.5 - 0.5 * x3 ) - x5 * ( 1.0 * x3 + 0.5 ) +
                           3.0 * std::pow( M_PI, 2 ) * ( x2 + 2 ) * std::sin( x1 ) + 10;
                 },
                 []( const hyteg::Point3D& ) { return 0; } ),
             std::make_shared< EGP0EpsilonOperatorStokesNitscheBC >( storage, minLevel, maxLevel, mu ),
             storage,
             minLevel,
             maxLevel,
             2,
             false,
             false,
             NULL,
             NULL,
             NULL,
             1 );
      }

      // cube_6el, inhom. solution, P2P1
      if ( true )
      {
         auto resNormsP2P1 = { 1e-5, 1e-7, 1e-7, 1e-5, 1e-6, 1e-7, 1e-5, 1e-6 };

         auto meshInfo = hyteg::MeshInfo::meshCuboid( Point3D( { -1, -1, -1 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );

         hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                                    walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
         setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
         auto storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, 1 );

         // cube, hom. solution
         WALBERLA_LOG_INFO_ON_ROOT( "### P2P1, cube_6el, IncreasingJump3D, alpha = 1 ###" );

         StokesConvergenceOrderTest< P2P1ElementwiseAffineEpsilonStokesOperator >(
             "P2P1EpsilonStokesOp3D_IncreasingJump3D_alpha5",
             std::make_tuple(
                 []( const hyteg::Point3D& p ) {
                    const real_t x = p[0];
                    const real_t y = p[1];
                    const real_t z = p[2];
                    return y + std::sin( M_PI * ( x + y + z ) );
                 },
                 []( const hyteg::Point3D& p ) {
                    const real_t x = p[0];
                    const real_t y = p[1];
                    const real_t z = p[2];

                    return z - 2 * std::sin( M_PI * ( x + y + z ) );
                 },
                 []( const hyteg::Point3D& p ) {
                    const real_t x = p[0];
                    const real_t y = p[1];
                    const real_t z = p[2];

                    return x + std::sin( M_PI * ( x + y + z ) );
                 },
                 []( const hyteg::Point3D& p ) {
                    const real_t x = p[0];
                    const real_t y = p[1];
                    const real_t z = p[2];
                    return 2 * x - 4 * y + 10 * z;
                 } ),
             std::make_tuple(
                 [&alpha]( const hyteg::Point3D& p ) {
                    const real_t x  = p[0];
                    const real_t y  = p[1];
                    const real_t z  = p[2];
                    const real_t x0 = x + y + z;
                    const real_t x1 = M_PI * x0;
                    const real_t x2 = std::tanh( alpha * x0 - 1.5 );
                    const real_t x3 = M_PI * std::cos( x1 );
                    const real_t x4 = alpha * ( 1 - std::pow( x2, 2 ) );
                    const real_t x5 = 2 * x4;
                    return -2.0 * x3 * x4 - x5 * ( 0.5 - 0.5 * x3 ) - x5 * ( 1.0 * x3 + 0.5 ) +
                           3.0 * std::pow( M_PI, 2 ) * ( x2 + 2 ) * std::sin( x1 ) + 2;
                 },
                 [&alpha]( const hyteg::Point3D& p ) {
                    const real_t x  = p[0];
                    const real_t y  = p[1];
                    const real_t z  = p[2];
                    const real_t x0 = x + y + z;
                    const real_t x1 = M_PI * x0;
                    const real_t x2 = std::tanh( alpha * x0 - 1.5 );
                    const real_t x3 = std::cos( x1 );
                    const real_t x4 = 1 - std::pow( x2, 2 );
                    return 4.0 * M_PI * alpha * x3 * x4 - 4 * alpha * x4 * ( -0.5 * M_PI * x3 + 0.5 ) -
                           6.0 * std::pow( M_PI, 2 ) * ( x2 + 2 ) * std::sin( x1 ) - 4;
                 },
                 [&alpha]( const hyteg::Point3D& p ) {
                    const real_t x  = p[0];
                    const real_t y  = p[1];
                    const real_t z  = p[2];
                    const real_t x0 = x + y + z;
                    const real_t x1 = M_PI * x0;
                    const real_t x2 = std::tanh( alpha * x0 - 1.5 );
                    const real_t x3 = M_PI * std::cos( x1 );
                    const real_t x4 = alpha * ( 1 - std::pow( x2, 2 ) );
                    const real_t x5 = 2 * x4;
                    return -2.0 * x3 * x4 - x5 * ( 0.5 - 0.5 * x3 ) - x5 * ( 1.0 * x3 + 0.5 ) +
                           3.0 * std::pow( M_PI, 2 ) * ( x2 + 2 ) * std::sin( x1 ) + 10;
                 },
                 []( const hyteg::Point3D& ) { return 0; } ),
             std::make_shared< P2P1ElementwiseAffineEpsilonStokesOperator >( storage, minLevel, maxLevel, mu ),
             storage,
             minLevel,
             maxLevel,
             2,
             false,
             false,
             NULL,
             NULL,
             NULL,
             1 );
      }
   }
}

void StraightJump3D( const uint_t minLevel, const uint_t maxLevel )
{
   auto mu = []( const hyteg::Point3D& p ) {
      const real_t z = p[2];
      return ( ( z < 0.5 ) ? ( 1 ) : ( 100 ) );
   };

   // cube_6el, inhom. solution, P2P1
   if ( false )
   {
      auto resNormsP2P1 = { 1e-5, 1e-7, 1e-7, 1e-5, 1e-6 };

      auto meshInfo = hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_6el.msh" );

      hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                                 walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
      auto storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, 1 );

      // cube, hom. solution
      WALBERLA_LOG_INFO_ON_ROOT( "### cube_6el, StraightJump3D, P2P1 ###" );

      hyteg::dg::eg::StokesConvergenceOrderTest< P2P1ElementwiseAffineEpsilonStokesOperator >(
          "P2P1EpsilonStokesOp3D_StraightJump3D_inhom",
          std::make_tuple(
              []( const hyteg::Point3D& p ) {
                 const real_t x = p[0];
                 const real_t y = p[1];
                 const real_t z = p[2];
                 return z + 2 * std::sin( M_PI * ( x + y + z ) ) + 4;
              },
              []( const hyteg::Point3D& p ) {
                 const real_t x = p[0];
                 const real_t y = p[1];
                 const real_t z = p[2];
                 return x - std::sin( M_PI * ( x + y + z ) ) + 4;
              },
              []( const hyteg::Point3D& p ) {
                 const real_t x = p[0];
                 const real_t y = p[1];
                 const real_t z = p[2];
                 return y - std::sin( M_PI * ( x + y + z ) ) + 4;
              },
              []( const hyteg::Point3D& xx ) { return std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] ); } ),
          std::make_tuple(
              []( const hyteg::Point3D& p ) {
                 const real_t x = p[0];
                 const real_t y = p[1];
                 const real_t z = p[2];
                 return 6.0 * std::pow( M_PI, 2 ) * ( ( z < 0.5 ) ? ( 1 ) : ( 100 ) ) * std::sin( M_PI * ( x + y + z ) ) +
                        4 * std::sin( 8 * y ) * std::sin( 2 * z ) * std::cos( 4 * x );
              },
              []( const hyteg::Point3D& p ) {
                 const real_t x = p[0];
                 const real_t y = p[1];
                 const real_t z = p[2];
                 return -3.0 * std::pow( M_PI, 2 ) * ( ( z < 0.5 ) ? ( 1 ) : ( 100 ) ) * std::sin( M_PI * ( x + y + z ) ) +
                        8 * std::sin( 4 * x ) * std::sin( 2 * z ) * std::cos( 8 * y );
              },
              []( const hyteg::Point3D& p ) {
                 const real_t x = p[0];
                 const real_t y = p[1];
                 const real_t z = p[2];
                 return -3.0 * std::pow( M_PI, 2 ) * ( ( z < 0.5 ) ? ( 1 ) : ( 100 ) ) * std::sin( M_PI * ( x + y + z ) ) +
                        2 * std::sin( 4 * x ) * std::sin( 8 * y ) * std::cos( 2 * z );
              },
              []( const hyteg::Point3D& ) { return 0; } ),
          std::make_shared< P2P1ElementwiseAffineEpsilonStokesOperator >( storage, minLevel, maxLevel, mu ),
          storage,
          minLevel,
          maxLevel,
          1,
          true,
          false,
          NULL,
          std::make_shared< std::vector< real_t > >( resNormsP2P1 ) );
   }

   / cube_6el, inhom.solution, EGP0 if ( true )
   {
      auto meshInfo = hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_6el.msh" );

      hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                                 walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
      auto storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, 1 );

      WALBERLA_LOG_INFO_ON_ROOT( "### cube_6el, StraightJump3D, Nitsche Bc ###" );
      auto resNormsEGP0 = { 1e-5, 1e-7, 1e-7, 1e-6 };
      hyteg::dg::eg::StokesConvergenceOrderTest< hyteg::dg::eg::EGP0EpsilonOperatorStokesNitscheBC >(
          "EGP0EpsilonStokesOp3DNitscheBC_StraightJump3D_inhom",
          std::make_tuple(
              []( const hyteg::Point3D& p ) {
                 const real_t x = p[0];
                 const real_t y = p[1];
                 const real_t z = p[2];
                 return z + 2 * std::sin( M_PI * ( x + y + z ) ) + 4;
              },
              []( const hyteg::Point3D& p ) {
                 const real_t x = p[0];
                 const real_t y = p[1];
                 const real_t z = p[2];
                 return x - std::sin( M_PI * ( x + y + z ) ) + 4;
              },
              []( const hyteg::Point3D& p ) {
                 const real_t x = p[0];
                 const real_t y = p[1];
                 const real_t z = p[2];
                 return y - std::sin( M_PI * ( x + y + z ) ) + 4;
              },
              []( const hyteg::Point3D& xx ) { return std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] ); } ),
          std::make_tuple(
              []( const hyteg::Point3D& p ) {
                 const real_t x = p[0];
                 const real_t y = p[1];
                 const real_t z = p[2];
                 return 6.0 * std::pow( M_PI, 2 ) * ( ( z < 0.5 ) ? ( 1 ) : ( 100 ) ) * std::sin( M_PI * ( x + y + z ) ) +
                        4 * std::sin( 8 * y ) * std::sin( 2 * z ) * std::cos( 4 * x );
              },
              []( const hyteg::Point3D& p ) {
                 const real_t x = p[0];
                 const real_t y = p[1];
                 const real_t z = p[2];
                 return -3.0 * std::pow( M_PI, 2 ) * ( ( z < 0.5 ) ? ( 1 ) : ( 100 ) ) * std::sin( M_PI * ( x + y + z ) ) +
                        8 * std::sin( 4 * x ) * std::sin( 2 * z ) * std::cos( 8 * y );
              },
              []( const hyteg::Point3D& p ) {
                 const real_t x = p[0];
                 const real_t y = p[1];
                 const real_t z = p[2];
                 return -3.0 * std::pow( M_PI, 2 ) * ( ( z < 0.5 ) ? ( 1 ) : ( 100 ) ) * std::sin( M_PI * ( x + y + z ) ) +
                        2 * std::sin( 4 * x ) * std::sin( 8 * y ) * std::cos( 2 * z );
              },
              []( const hyteg::Point3D& ) { return 0; } ),
          std::make_shared< hyteg::dg::eg::EGP0EpsilonOperatorStokesNitscheBC >( storage, minLevel, maxLevel, mu ),
          storage,
          minLevel,
          maxLevel,
          1,
          false,
          false,
          NULL,
          std::make_shared< std::vector< real_t > >( resNormsEGP0 ) );
   }
}
} // namespace eg
} // namespace dg
} // namespace hyteg
