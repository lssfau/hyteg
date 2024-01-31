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
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/MinresSolver.hpp"

#include "mixedOperator/EGConvTestUtils.hpp"
#include "mixedOperator/EGOperators.hpp"
#include "constantStencilOperator/P1ConstantOperator.cpp"

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

void MatfreeLaplaceConvTest2D( const uint_t                               minLevel,
                               const uint_t                               maxLevel,
                               const std::shared_ptr< PrimitiveStorage >& storage,
                               const uint_t                               solver );
void MatfreeEpsilonConvTest2D( const uint_t                               minLevel,
                               const uint_t                               maxLevel,
                               const std::shared_ptr< PrimitiveStorage >& storage,
                               const uint_t                               solver );

} // namespace eg
} // namespace dg
} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::PETScManager petscManager( &argc, &argv );

   if constexpr ( true )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Testing matfree Laplace solver convergence 2D ###" )
      auto meshInfo = hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" );

      hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                                 walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
      auto storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, 1 );

      hyteg::dg::eg::MatfreeEpsilonConvTest2D( 3, 3, storage, 6 );
      hyteg::dg::eg::MatfreeEpsilonConvTest2D( 3, 3, storage, 0 );

      hyteg::dg::eg::MatfreeLaplaceConvTest2D( 3, 3, storage, 0 );
      hyteg::dg::eg::MatfreeLaplaceConvTest2D( 3, 3, storage, 6 );
   }

   return 0;
}

namespace hyteg {
namespace dg {
namespace eg {

void MatfreeLaplaceConvTest2D( const uint_t                               minLevel,
                               const uint_t                               maxLevel,
                               const std::shared_ptr< PrimitiveStorage >& storage,
                               const uint_t                               solver )
{
   auto dummyLambda = []( const Point3D& ) -> real_t { return 0; };

   // nitsche BCs
   WALBERLA_LOG_INFO_ON_ROOT( "### Nitsche BCs, quad_4el, inhom. solution ###" );
   hyteg::dg::eg::StokesConvergenceOrderTest< hyteg::dg::eg::EGP0StokesOperatorNitscheBC >(
       "EGP0StokesOpNitscheBC2D_matfree",
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return y + 2 * std::sin( M_PI * ( x + y ) ) + 4;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return -x - 2 * std::sin( M_PI * ( x + y ) ) + 3;
           },
           dummyLambda,
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return 2 * x - y + 1;
           } ),
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return 4 * std::pow( M_PI, 2 ) * std::sin( M_PI * ( x + y ) ) + 2;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return -4 * std::pow( M_PI, 2 ) * std::sin( M_PI * ( x + y ) ) - 1;
           },
           dummyLambda,
           []( const Point3D& ) -> real_t { return 0; } ),
       std::make_shared< hyteg::dg::eg::EGP0StokesOperatorNitscheBC >( storage, minLevel, maxLevel ),
       storage,
       minLevel,
       maxLevel,
       solver,
       1e-6,
       false,
       std::make_pair( true, 2.0e-02 ) );
}

void MatfreeEpsilonConvTest2D( const uint_t                               minLevel,
                               const uint_t                               maxLevel,
                               const std::shared_ptr< PrimitiveStorage >& storage,
                               const uint_t                               solver )
{
   auto dummyLambda = []( const Point3D& ) -> real_t { return 0; };

   // nitsche BCs
   WALBERLA_LOG_INFO_ON_ROOT( "### Checking matfree convergence for Nitsche BCs, quad_4el, inhom. solution ###" );
   hyteg::dg::eg::StokesConvergenceOrderTest< hyteg::dg::eg::EGP0EpsilonOperatorStokesNitscheBC >(
       "EGP0EpsilonOpNitscheBC2D_divFree_smoothVisc",
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return y + 2 * std::sin( M_PI * ( x + y ) ) + 4;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return -x - 2 * std::sin( M_PI * ( x + y ) ) + 3;
           },
           dummyLambda,
           []( const Point3D& p ) -> real_t {
              const real_t x = p[0];
              const real_t y = p[1];
              return 2 * x - y + 1;
           } ),
       std::make_tuple(
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = M_PI * ( x + y );
              const real_t x1 = M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 );
              const real_t x2 = std::exp( x ) * std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) );
              const real_t x3 = x2 * std::sin( x1 );
              return 4.0 * std::pow( M_PI, 2 ) * ( x3 + 1 ) * std::sin( x0 ) -
                     4.0 * M_PI * ( ( 1.0 / 2.0 ) * M_PI * x2 * std::cos( x1 ) + x3 ) * std::cos( x0 ) + 2;
           },
           []( const Point3D& p ) -> real_t {
              const real_t x  = p[0];
              const real_t y  = p[1];
              const real_t x0 = std::exp( x );
              const real_t x1 = std::pow( M_PI, 2 );
              const real_t x2 = M_PI * ( x + y );
              const real_t x3 = std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) );
              const real_t x4 = M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 );
              return 2.0 * x0 * x1 * x3 * std::cos( x2 ) * std::cos( x4 ) -
                     4.0 * x1 * ( x0 * x3 * std::sin( x4 ) + 1 ) * std::sin( x2 ) - 1;
           },
           dummyLambda,
           dummyLambda ),
       std::make_shared< hyteg::dg::eg::EGP0EpsilonOperatorStokesNitscheBC >(
           storage,
           minLevel,
           maxLevel,
           []( const hyteg::Point3D& p ) {
              const real_t x = p[0];
              const real_t y = p[1];
              return std::exp( x ) * std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) *
                         std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) ) +
                     1;
           } ),
       storage,
       minLevel,
       maxLevel,
       solver,
       1e-6,
       false,
       std::make_pair( true, 2.0e-02 ) );
}
} // namespace eg
} // namespace dg
} // namespace hyteg
