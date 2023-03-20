/*
* Copyright (c) 2017-2022 Nils Kohl.
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
#include "hyteg/composites/P1P0StokesOperator.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
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
void Stokes2D( const uint_t minLevel, const uint_t maxLevel );

void HytegSolverCheck2D( const uint_t minLevel, const uint_t maxLevel, const std::shared_ptr< PrimitiveStorage >& storage );

void Stokes3D( const uint_t minLevel, const uint_t maxLevel );
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

   uint_t minLevel = 2;

   if constexpr ( true )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Testing 2D ###" )

      hyteg::dg::eg::Stokes2D( minLevel, 7 );
   }

   if constexpr ( false )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Testing 3D ###" )

      hyteg::dg::eg::Stokes3D( minLevel, 5 );
   }

   return 0;
}

namespace hyteg {
namespace dg {
namespace eg {

void Stokes2D( const uint_t minLevel, const uint_t maxLevel )
{
   auto dummyLambda = []( const Point3D& ) -> real_t { return 0; };
    if constexpr ( true )
    {
        auto                         meshInfo = hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
        hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                                   walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
        setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
        auto storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, 1 );


        auto discrErrorsEGP0Nitsche = { 0.0187979, 0.00489507, 0.00124212, 0.00031431, 7.83271e-05 };

        auto resNormsEGP0 = { 1e-5, 1e-5, 1e-5, 1e-5, 1e-6, 1e-6 };



        if(true) {
            StokesConvergenceOrderTest<EGP0StokesOperatorNitscheBC>(
                    "EGP0StokesOpNitscheBC2D_quad_4el_inhom",
                    std::make_tuple(
                            [](const Point3D &p) -> real_t {
                                const real_t x = p[0];
                                const real_t y = p[1];
                                return y + 2 * std::sin(M_PI * (x + y)) + 4;
                            },
                            [](const Point3D &p) -> real_t {
                                const real_t x = p[0];
                                const real_t y = p[1];
                                return -x - 2 * std::sin(M_PI * (x + y)) + 3;
                            },
                            dummyLambda,
                            [](const Point3D &p) -> real_t {
                                const real_t x = p[0];
                                const real_t y = p[1];
                                return 2 * x - y + 1;
                            }),
                    std::make_tuple(
                            [](const Point3D &p) -> real_t {
                                const real_t x = p[0];
                                const real_t y = p[1];
                                return 4 * std::pow(M_PI, 2) * std::sin(M_PI * (x + y)) + 2;
                            },
                            [](const Point3D &p) -> real_t {
                                const real_t x = p[0];
                                const real_t y = p[1];
                                return -4 * std::pow(M_PI, 2) * std::sin(M_PI * (x + y)) - 1;
                            },
                            dummyLambda,
                            [](const Point3D &) -> real_t { return 0; }),
                    std::make_shared<EGP0StokesOperatorNitscheBC>(storage, minLevel, maxLevel),
                    storage,
                    minLevel,
                    maxLevel,
                    2,
                    false, false, std::make_shared<std::vector<real_t>>(discrErrorsEGP0Nitsche),
                    std::make_shared<std::vector<real_t>>(resNormsEGP0));
        }
    }
   // quad_4el, inhom
   if constexpr ( false )
   {
      auto                         meshInfo = hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" );
      hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                                 walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
      auto storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, 1 );


      auto discrErrorsP2P1       = { 0.000454598, 5.78939e-05, 7.31194e-06, 9.18941e-07, 1.15185e-07 };
    auto discrErrorsEGP0Nitsche = { 0.0187979, 0.00489507, 0.00124212, 0.00031431, 7.83271e-05 };

      auto resNormsEGP0 = { 1e-5, 1e-5, 1e-5, 1e-5, 1e-6, 1e-6 };

      auto resNormsP2P1 = { 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-6 };



        if(false) {
            StokesConvergenceOrderTest<P2P1TaylorHoodStokesOperator>(
                    "P2P1StokesOp2D_quad_4el_inhom",
                    std::make_tuple(
                            [](const Point3D &p) -> real_t {
                                const real_t x = p[0];
                                const real_t y = p[1];
                                return y + 2 * std::sin(M_PI * (x + y)) + 4;
                            },
                            [](const Point3D &p) -> real_t {
                                const real_t x = p[0];
                                const real_t y = p[1];
                                return -x - 2 * std::sin(M_PI * (x + y)) + 3;
                            },
                            dummyLambda,
                            [](const Point3D &p) -> real_t {
                                const real_t x = p[0];
                                const real_t y = p[1];
                                return 2 * x - y + 1;
                            }),
                    std::make_tuple(
                            [](const Point3D &p) -> real_t {
                                const real_t x = p[0];
                                const real_t y = p[1];
                                return 4 * std::pow(M_PI, 2) * std::sin(M_PI * (x + y)) + 2;
                            },
                            [](const Point3D &p) -> real_t {
                                const real_t x = p[0];
                                const real_t y = p[1];
                                return -4 * std::pow(M_PI, 2) * std::sin(M_PI * (x + y)) - 1;
                            },
                            dummyLambda,
                            [](const Point3D &) -> real_t { return 0; }),
                    std::make_shared<P2P1TaylorHoodStokesOperator>(storage, minLevel, maxLevel),
                    storage,
                    minLevel,
                    maxLevel,
                    2,
                    false,
                    false,
                    std::make_shared<std::vector<real_t> >(discrErrorsP2P1),
                    std::make_shared<std::vector<real_t>>(resNormsP2P1));
        }

            if(true) {
                StokesConvergenceOrderTest<EGP0StokesOperatorNitscheBC>(
                        "EGP0StokesOpNitscheBC2D_quad_4el_inhom",
                        std::make_tuple(
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    return y + 2 * std::sin(M_PI * (x + y)) + 4;
                                },
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    return -x - 2 * std::sin(M_PI * (x + y)) + 3;
                                },
                                dummyLambda,
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    return 2 * x - y + 1;
                                }),
                        std::make_tuple(
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    return 4 * std::pow(M_PI, 2) * std::sin(M_PI * (x + y)) + 2;
                                },
                                [](const Point3D &p) -> real_t {
                                    const real_t x = p[0];
                                    const real_t y = p[1];
                                    return -4 * std::pow(M_PI, 2) * std::sin(M_PI * (x + y)) - 1;
                                },
                                dummyLambda,
                                [](const Point3D &) -> real_t { return 0; }),
                        std::make_shared<EGP0StokesOperatorNitscheBC>(storage, minLevel, maxLevel),
                        storage,
                        minLevel,
                        maxLevel,
                        2,
                        false, false, std::make_shared<std::vector<real_t>>(discrErrorsEGP0Nitsche),
                        std::make_shared<std::vector<real_t>>(resNormsEGP0));
            }
   }

   // IIPG: quad_4el, inhom
   if constexpr ( false )
   {
      auto                         meshInfo = hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" );
      hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                                 walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
      auto storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, 1 );

      EGP0IIPGStokesOperator EGP0StokesOp( storage, minLevel, maxLevel );
      StokesConvergenceOrderTest< EGP0IIPGStokesOperator >(
          "EGP0IIPGStokesOp2D_quad_4el_inhom",
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
          std::make_shared< EGP0IIPGStokesOperator >( storage, minLevel, maxLevel ),
          storage,
          minLevel,
          maxLevel );
   }

   // quad_16el, inhom
   if constexpr ( false )
   {
      auto                         meshInfo = hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_16el.msh" );
      hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                                 walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
      auto storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, 1 );

      StokesConvergenceOrderTest< EGP0StokesOperator >( "EGP0StokesOp2D_quad_16el_inhom",
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
                                                        std::make_shared< EGP0StokesOperator >( storage, minLevel, maxLevel ),
                                                        storage,
                                                        minLevel,
                                                        maxLevel );
   }

   // IIPG: quad_16el, inhom
   if constexpr ( false )
   {
      auto                         meshInfo = hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_16el.msh" );
      hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                                 walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
      auto storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, 1 );

      EGP0IIPGStokesOperator EGP0StokesOp( storage, minLevel, maxLevel );
      StokesConvergenceOrderTest< EGP0IIPGStokesOperator >(
          "EGP0IIPGStokesOp2D_quad_16el_inhom",
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
          std::make_shared< EGP0IIPGStokesOperator >( storage, minLevel, maxLevel ),
          storage,
          minLevel,
          maxLevel );
   }
}

void Stokes3D( const uint_t minLevel, const uint_t maxLevel )
{
   // tet
   {
      auto                         meshInfo = hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" );
      hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                                 walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
      auto storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, 1 );

      hyteg::P2P1TaylorHoodStokesOperator P2P1StokesOp( storage, minLevel, maxLevel + 2 );

      // tet_1el, simple
      if ( false )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "### tet_1el, simple ###" );
         StokesConvergenceOrderTest< EGP0StokesOperator >(
             "EGP0StokesOp3D_tet_1el_simple",
             std::make_tuple( []( const Point3D& xx ) -> real_t { return 2 * xx[0]; },
                              []( const Point3D& xx ) -> real_t { return -3 * xx[1]; },
                              []( const Point3D& xx ) -> real_t { return 1 * xx[2]; },
                              []( const Point3D& ) -> real_t { return 1; } ),
             std::make_tuple( []( const Point3D& ) -> real_t { return 0; },
                              []( const Point3D& ) -> real_t { return 0; },
                              []( const Point3D& ) -> real_t { return 0; },
                              []( const Point3D& ) -> real_t { return 0; } ),
             std::make_shared< EGP0StokesOperator >( storage, minLevel, maxLevel ),
             storage,
             minLevel,
             maxLevel,
             5,
             true );
      }

      // tet_1el, hom.
      if ( false )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "### tet_1el, hom. ###" );
         auto u = []( const Point3D& xx ) -> real_t {
            const real_t x  = xx[0];
            const real_t y  = xx[1];
            const real_t z  = xx[2];
            const real_t x0 = 2 * x;
            const real_t x1 = 2 * y;
            const real_t x2 = 2 * z;
            return std::sin( M_PI * x0 ) * std::sin( M_PI * x1 ) * std::sin( M_PI * x2 ) *
                   std::sin( M_PI * ( x0 + x1 + x2 - 2 ) );
         };

         StokesConvergenceOrderTest< EGP0StokesOperator >(
             "EGP0StokesOp3D_tet_1el_hom",
             std::make_tuple( u,
                              u,
                              u,
                              []( const Point3D& xx ) -> real_t {
                                 const real_t x = xx[0];
                                 const real_t y = xx[1];
                                 const real_t z = xx[2];
                                 return std::sin( 4 * x ) * std::sin( 8 * y ) * std::sin( 2 * z );
                              } ),
             std::make_tuple(
                 []( const Point3D& xx ) -> real_t {
                    const real_t x   = xx[0];
                    const real_t y   = xx[1];
                    const real_t z   = xx[2];
                    const real_t x0  = 2 * z;
                    const real_t x1  = std::pow( M_PI, 2 );
                    const real_t x2  = 2 * x;
                    const real_t x3  = M_PI * x2;
                    const real_t x4  = std::sin( x3 );
                    const real_t x5  = 2 * y;
                    const real_t x6  = M_PI * x5;
                    const real_t x7  = std::sin( x6 );
                    const real_t x8  = M_PI * x0;
                    const real_t x9  = std::sin( x8 );
                    const real_t x10 = M_PI * ( x0 + x2 + x5 - 2 );
                    const real_t x11 = 8 * std::cos( x10 );
                    const real_t x12 = x1 * x11 * x4;
                    return -x1 * x11 * x7 * x9 * std::cos( x3 ) + 24 * x1 * x4 * x7 * x9 * std::sin( x10 ) -
                           x12 * x7 * std::cos( x8 ) - x12 * x9 * std::cos( x6 ) +
                           4 * std::sin( x0 ) * std::sin( 8 * y ) * std::cos( 4 * x );
                 },
                 []( const Point3D& xx ) -> real_t {
                    const real_t x   = xx[0];
                    const real_t y   = xx[1];
                    const real_t z   = xx[2];
                    const real_t x0  = 2 * z;
                    const real_t x1  = std::pow( M_PI, 2 );
                    const real_t x2  = 2 * x;
                    const real_t x3  = M_PI * x2;
                    const real_t x4  = std::sin( x3 );
                    const real_t x5  = 2 * y;
                    const real_t x6  = M_PI * x5;
                    const real_t x7  = std::sin( x6 );
                    const real_t x8  = M_PI * x0;
                    const real_t x9  = std::sin( x8 );
                    const real_t x10 = M_PI * ( x0 + x2 + x5 - 2 );
                    const real_t x11 = 8 * std::cos( x10 );
                    const real_t x12 = x1 * x11 * x4;
                    return -x1 * x11 * x7 * x9 * std::cos( x3 ) + 24 * x1 * x4 * x7 * x9 * std::sin( x10 ) -
                           x12 * x7 * std::cos( x8 ) - x12 * x9 * std::cos( x6 ) +
                           8 * std::sin( 4 * x ) * std::sin( x0 ) * std::cos( 8 * y );
                 },
                 []( const Point3D& xx ) -> real_t {
                    const real_t x   = xx[0];
                    const real_t y   = xx[1];
                    const real_t z   = xx[2];
                    const real_t x0  = 2 * z;
                    const real_t x1  = std::pow( M_PI, 2 );
                    const real_t x2  = 2 * x;
                    const real_t x3  = M_PI * x2;
                    const real_t x4  = std::sin( x3 );
                    const real_t x5  = 2 * y;
                    const real_t x6  = M_PI * x5;
                    const real_t x7  = std::sin( x6 );
                    const real_t x8  = M_PI * x0;
                    const real_t x9  = std::sin( x8 );
                    const real_t x10 = M_PI * ( x0 + x2 + x5 - 2 );
                    const real_t x11 = 8 * std::cos( x10 );
                    const real_t x12 = x1 * x11 * x4;
                    return -x1 * x11 * x7 * x9 * std::cos( x3 ) + 24 * x1 * x4 * x7 * x9 * std::sin( x10 ) -
                           x12 * x7 * std::cos( x8 ) - x12 * x9 * std::cos( x6 ) +
                           2 * std::sin( 4 * x ) * std::sin( 8 * y ) * std::cos( x0 );
                 },
                 []( const Point3D& xx ) -> real_t {
                    const real_t x   = xx[0];
                    const real_t y   = xx[1];
                    const real_t z   = xx[2];
                    const real_t x0  = 2 * x;
                    const real_t x1  = 2 * y;
                    const real_t x2  = 2 * z;
                    const real_t x3  = M_PI * ( x0 + x1 + x2 - 2 );
                    const real_t x4  = M_PI * x0;
                    const real_t x5  = std::sin( x4 );
                    const real_t x6  = M_PI * x1;
                    const real_t x7  = std::sin( x6 );
                    const real_t x8  = M_PI * x2;
                    const real_t x9  = std::sin( x8 );
                    const real_t x10 = M_PI * x7 * x9;
                    const real_t x11 = 2 * std::sin( x3 );
                    const real_t x12 = M_PI * x11 * x5;
                    return -x10 * x11 * std::cos( x4 ) - 6 * x10 * x5 * std::cos( x3 ) - x12 * x7 * std::cos( x8 ) -
                           x12 * x9 * std::cos( x6 );
                 } ),

             std::make_shared< EGP0StokesOperator >( storage, minLevel, maxLevel ),
             storage,
             minLevel,
             maxLevel + 1,
             5,
             true );
      }

      // tet_1el, hom. exp-enhanced
      if ( false )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "### tet_1el, hom. exp-enhanced ###" );
         auto u = []( const Point3D& xx ) -> real_t {
            const real_t x  = xx[0];
            const real_t y  = xx[1];
            const real_t z  = xx[2];
            const real_t x0 = 2 * x;
            const real_t x1 = 2 * y;
            const real_t x2 = 2 * z;
            return std::exp( 5 * z ) * std::sin( M_PI * x0 ) * std::sin( M_PI * x1 ) * std::sin( M_PI * x2 ) *
                   std::sin( M_PI * ( x0 + x1 + x2 - 2 ) );
         };

         auto v = []( const Point3D& xx ) -> real_t {
            const real_t x  = xx[0];
            const real_t y  = xx[1];
            const real_t z  = xx[2];
            const real_t x0 = 2 * x;
            const real_t x1 = 2 * y;
            const real_t x2 = 2 * z;
            return std::exp( 10 * x ) * std::sin( M_PI * x0 ) * std::sin( M_PI * x1 ) * std::sin( M_PI * x2 ) *
                   std::sin( M_PI * ( x0 + x1 + x2 - 2 ) );
         };

         auto w = []( const Point3D& xx ) -> real_t {
            const real_t x  = xx[0];
            const real_t y  = xx[1];
            const real_t z  = xx[2];
            const real_t x0 = 2 * x;
            const real_t x1 = 2 * y;
            const real_t x2 = 2 * z;
            return std::exp( 15 * y ) * std::sin( M_PI * x0 ) * std::sin( M_PI * x1 ) * std::sin( M_PI * x2 ) *
                   std::sin( M_PI * ( x0 + x1 + x2 - 2 ) );
         };

         StokesConvergenceOrderTest< EGP0StokesOperator >(
             "EGP0StokesOp3D_tet_1el_hom_exp",
             std::make_tuple( u,
                              v,
                              w,
                              []( const Point3D& xx ) -> real_t {
                                 const real_t x = xx[0];
                                 const real_t y = xx[1];
                                 const real_t z = xx[2];
                                 return std::sin( 4 * x ) * std::sin( 8 * y ) * std::sin( 2 * z );
                              } ),
             std::make_tuple(
                 []( const Point3D& xx ) -> real_t {
                    const real_t x   = xx[0];
                    const real_t y   = xx[1];
                    const real_t z   = xx[2];
                    const real_t x0  = 2 * z;
                    const real_t x1  = 2 * x;
                    const real_t x2  = 2 * y;
                    const real_t x3  = M_PI * ( x0 + x1 + x2 - 2 );
                    const real_t x4  = std::sin( x3 );
                    const real_t x5  = M_PI * x1;
                    const real_t x6  = std::sin( x5 );
                    const real_t x7  = std::exp( 5 * z );
                    const real_t x8  = M_PI * x2;
                    const real_t x9  = std::sin( x8 );
                    const real_t x10 = M_PI * x0;
                    const real_t x11 = std::sin( x10 );
                    const real_t x12 = x11 * x7 * x9;
                    const real_t x13 = x12 * x6;
                    const real_t x14 = std::cos( x3 );
                    const real_t x15 = 20 * M_PI;
                    const real_t x16 = x6 * x7;
                    const real_t x17 = x16 * x9 * std::cos( x10 );
                    const real_t x18 = std::pow( M_PI, 2 );
                    const real_t x19 = 8 * x14 * x18;
                    return -x11 * x16 * x19 * std::cos( x8 ) + 24 * x11 * x18 * x4 * x6 * x7 * x9 - x12 * x19 * std::cos( x5 ) -
                           x13 * x14 * x15 - 25 * x13 * x4 - x15 * x17 * x4 - x17 * x19 +
                           4 * std::sin( x0 ) * std::sin( 8 * y ) * std::cos( 4 * x );
                 },
                 []( const Point3D& xx ) -> real_t {
                    const real_t x   = xx[0];
                    const real_t y   = xx[1];
                    const real_t z   = xx[2];
                    const real_t x0  = 2 * z;
                    const real_t x1  = 2 * x;
                    const real_t x2  = 2 * y;
                    const real_t x3  = M_PI * ( x0 + x1 + x2 - 2 );
                    const real_t x4  = std::sin( x3 );
                    const real_t x5  = M_PI * x1;
                    const real_t x6  = std::sin( x5 );
                    const real_t x7  = std::exp( 10 * x );
                    const real_t x8  = M_PI * x2;
                    const real_t x9  = std::sin( x8 );
                    const real_t x10 = M_PI * x0;
                    const real_t x11 = std::sin( x10 );
                    const real_t x12 = x11 * x7 * x9;
                    const real_t x13 = x12 * x6;
                    const real_t x14 = std::cos( x3 );
                    const real_t x15 = 40 * M_PI;
                    const real_t x16 = x12 * std::cos( x5 );
                    const real_t x17 = std::pow( M_PI, 2 );
                    const real_t x18 = 8 * x14 * x17;
                    const real_t x19 = x18 * x6 * x7;
                    return 24 * x11 * x17 * x4 * x6 * x7 * x9 - x11 * x19 * std::cos( x8 ) - x13 * x14 * x15 - 100 * x13 * x4 -
                           x15 * x16 * x4 - x16 * x18 - x19 * x9 * std::cos( x10 ) +
                           8 * std::sin( 4 * x ) * std::sin( x0 ) * std::cos( 8 * y );
                 },
                 []( const Point3D& xx ) -> real_t {
                    const real_t x   = xx[0];
                    const real_t y   = xx[1];
                    const real_t z   = xx[2];
                    const real_t x0  = 2 * z;
                    const real_t x1  = 2 * x;
                    const real_t x2  = 2 * y;
                    const real_t x3  = M_PI * ( x0 + x1 + x2 - 2 );
                    const real_t x4  = std::sin( x3 );
                    const real_t x5  = M_PI * x1;
                    const real_t x6  = std::sin( x5 );
                    const real_t x7  = std::exp( 15 * y );
                    const real_t x8  = M_PI * x2;
                    const real_t x9  = std::sin( x8 );
                    const real_t x10 = M_PI * x0;
                    const real_t x11 = std::sin( x10 );
                    const real_t x12 = x11 * x7 * x9;
                    const real_t x13 = x12 * x6;
                    const real_t x14 = std::cos( x3 );
                    const real_t x15 = 60 * M_PI;
                    const real_t x16 = x6 * x7;
                    const real_t x17 = x11 * x16 * std::cos( x8 );
                    const real_t x18 = std::pow( M_PI, 2 );
                    const real_t x19 = 8 * x14 * x18;
                    return 24 * x11 * x18 * x4 * x6 * x7 * x9 - x12 * x19 * std::cos( x5 ) - x13 * x14 * x15 - 225 * x13 * x4 -
                           x15 * x17 * x4 - x16 * x19 * x9 * std::cos( x10 ) - x17 * x19 +
                           2 * std::sin( 4 * x ) * std::sin( 8 * y ) * std::cos( x0 );
                 },
                 []( const Point3D& xx ) -> real_t {
                    const real_t x   = xx[0];
                    const real_t y   = xx[1];
                    const real_t z   = xx[2];
                    const real_t x0  = std::exp( 10 * x );
                    const real_t x1  = 2 * y;
                    const real_t x2  = M_PI * x1;
                    const real_t x3  = std::sin( x2 );
                    const real_t x4  = 2 * x;
                    const real_t x5  = 2 * z;
                    const real_t x6  = M_PI * ( x1 + x4 + x5 - 2 );
                    const real_t x7  = M_PI * x4;
                    const real_t x8  = std::sin( x7 );
                    const real_t x9  = M_PI * x5;
                    const real_t x10 = std::sin( x9 );
                    const real_t x11 = 2 * M_PI * x10 * x8;
                    const real_t x12 = x11 * std::cos( x6 );
                    const real_t x13 = x12 * x3;
                    const real_t x14 = std::sin( x6 );
                    const real_t x15 = std::exp( 15 * y );
                    const real_t x16 = 2 * M_PI * x14;
                    const real_t x17 = x3 * std::exp( 5 * z );
                    return -x0 * x11 * x14 * std::cos( x2 ) - x0 * x13 - x10 * x16 * x17 * std::cos( x7 ) - x12 * x17 -
                           x13 * x15 - x15 * x16 * x3 * x8 * std::cos( x9 );
                 } ),

             std::make_shared< EGP0StokesOperator >( storage, minLevel, maxLevel ),
             storage,
             minLevel,
             maxLevel,
             5,
             true );
      }

      // tet_1el, inhom.
      if ( false )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "### tet_1el, inhom. ###" );

         auto u = []( const Point3D& xx ) -> real_t {
            const real_t x  = xx[0];
            const real_t y  = xx[1];
            const real_t z  = xx[2];
            const real_t x0 = 2 * x;
            const real_t x1 = 2 * y;
            const real_t x2 = 2 * z;
            return std::sin( M_PI * ( x0 + 1.0 ) ) * std::sin( M_PI * ( x1 + 1.0 ) ) * std::sin( M_PI * ( x2 + 1.0 ) ) *
                   std::sin( M_PI * ( x0 + x1 + x2 - 2 ) );
         };

         StokesConvergenceOrderTest< EGP0StokesOperator >(
             "EGP0StokesOp3D_tet_1el_inhom",
             std::make_tuple( u,
                              u,
                              u,
                              []( const Point3D& xx ) -> real_t {
                                 const real_t x = xx[0];
                                 const real_t y = xx[1];
                                 const real_t z = xx[2];
                                 return std::sin( 4 * x ) * std::sin( 8 * y ) * std::sin( 2 * z );
                              } ),
             std::make_tuple(
                 []( const Point3D& xx ) -> real_t {
                    const real_t x   = xx[0];
                    const real_t y   = xx[1];
                    const real_t z   = xx[2];
                    const real_t x0  = 2 * z;
                    const real_t x1  = std::pow( M_PI, 2 );
                    const real_t x2  = 2 * x;
                    const real_t x3  = M_PI * ( x2 + 1.0 );
                    const real_t x4  = std::sin( x3 );
                    const real_t x5  = 2 * y;
                    const real_t x6  = M_PI * ( x5 + 1.0 );
                    const real_t x7  = std::sin( x6 );
                    const real_t x8  = M_PI * ( x0 + 1.0 );
                    const real_t x9  = std::sin( x8 );
                    const real_t x10 = M_PI * ( x0 + x2 + x5 - 2 );
                    const real_t x11 = 8 * std::cos( x10 );
                    const real_t x12 = x1 * x11 * x4;
                    return -x1 * x11 * x7 * x9 * std::cos( x3 ) + 24 * x1 * x4 * x7 * x9 * std::sin( x10 ) -
                           x12 * x7 * std::cos( x8 ) - x12 * x9 * std::cos( x6 ) +
                           4 * std::sin( x0 ) * std::sin( 8 * y ) * std::cos( 4 * x );
                 },
                 []( const Point3D& xx ) -> real_t {
                    const real_t x   = xx[0];
                    const real_t y   = xx[1];
                    const real_t z   = xx[2];
                    const real_t x0  = 2 * z;
                    const real_t x1  = std::pow( M_PI, 2 );
                    const real_t x2  = 2 * x;
                    const real_t x3  = M_PI * ( x2 + 1.0 );
                    const real_t x4  = std::sin( x3 );
                    const real_t x5  = 2 * y;
                    const real_t x6  = M_PI * ( x5 + 1.0 );
                    const real_t x7  = std::sin( x6 );
                    const real_t x8  = M_PI * ( x0 + 1.0 );
                    const real_t x9  = std::sin( x8 );
                    const real_t x10 = M_PI * ( x0 + x2 + x5 - 2 );
                    const real_t x11 = 8 * std::cos( x10 );
                    const real_t x12 = x1 * x11 * x4;
                    return -x1 * x11 * x7 * x9 * std::cos( x3 ) + 24 * x1 * x4 * x7 * x9 * std::sin( x10 ) -
                           x12 * x7 * std::cos( x8 ) - x12 * x9 * std::cos( x6 ) +
                           8 * std::sin( 4 * x ) * std::sin( x0 ) * std::cos( 8 * y );
                 },
                 []( const Point3D& xx ) -> real_t {
                    const real_t x   = xx[0];
                    const real_t y   = xx[1];
                    const real_t z   = xx[2];
                    const real_t x0  = 2 * z;
                    const real_t x1  = std::pow( M_PI, 2 );
                    const real_t x2  = 2 * x;
                    const real_t x3  = M_PI * ( x2 + 1.0 );
                    const real_t x4  = std::sin( x3 );
                    const real_t x5  = 2 * y;
                    const real_t x6  = M_PI * ( x5 + 1.0 );
                    const real_t x7  = std::sin( x6 );
                    const real_t x8  = M_PI * ( x0 + 1.0 );
                    const real_t x9  = std::sin( x8 );
                    const real_t x10 = M_PI * ( x0 + x2 + x5 - 2 );
                    const real_t x11 = 8 * std::cos( x10 );
                    const real_t x12 = x1 * x11 * x4;
                    return -x1 * x11 * x7 * x9 * std::cos( x3 ) + 24 * x1 * x4 * x7 * x9 * std::sin( x10 ) -
                           x12 * x7 * std::cos( x8 ) - x12 * x9 * std::cos( x6 ) +
                           2 * std::sin( 4 * x ) * std::sin( 8 * y ) * std::cos( x0 );
                 },
                 []( const Point3D& xx ) -> real_t {
                    const real_t x   = xx[0];
                    const real_t y   = xx[1];
                    const real_t z   = xx[2];
                    const real_t x0  = 2 * x;
                    const real_t x1  = 2 * y;
                    const real_t x2  = 2 * z;
                    const real_t x3  = M_PI * ( x0 + x1 + x2 - 2 );
                    const real_t x4  = M_PI * ( x0 + 1.0 );
                    const real_t x5  = std::sin( x4 );
                    const real_t x6  = M_PI * ( x1 + 1.0 );
                    const real_t x7  = std::sin( x6 );
                    const real_t x8  = M_PI * ( x2 + 1.0 );
                    const real_t x9  = std::sin( x8 );
                    const real_t x10 = M_PI * x7 * x9;
                    const real_t x11 = 2 * std::sin( x3 );
                    const real_t x12 = M_PI * x11 * x5;
                    return -x10 * x11 * std::cos( x4 ) - 6 * x10 * x5 * std::cos( x3 ) - x12 * x7 * std::cos( x8 ) -
                           x12 * x9 * std::cos( x6 );
                 } ),

             std::make_shared< EGP0StokesOperator >( storage, minLevel, maxLevel ),
             storage,
             minLevel,
             maxLevel,
             5,
             true );
      }

      // tet_1el, inhom. split solution
      if ( false )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "### tet_1el, inhom. split solution ###" );
         StokesConvergenceOrderTest< EGP0StokesOperator >(
             "EGP0StokesOp3D_tet_1el_inhom_split",
             std::make_tuple( []( const Point3D& xx ) -> real_t { return -real_c( 4 ) * std::cos( real_c( 4 ) * xx[2] ); },
                              []( const Point3D& xx ) -> real_t { return real_c( 8 ) * std::cos( real_c( 8 ) * xx[0] ); },
                              []( const Point3D& xx ) -> real_t { return -real_c( 2 ) * std::cos( real_c( 2 ) * xx[1] ); },
                              []( const Point3D& xx ) -> real_t {
                                 return std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] );
                              } ),
             std::make_tuple(
                 []( const Point3D& xx ) -> real_t {
                    return 4 * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] ) * std::cos( 4 * xx[0] ) - 64 * std::cos( 4 * xx[2] );
                 },
                 []( const Point3D& xx ) -> real_t {
                    return 8 * std::sin( 4 * xx[0] ) * std::sin( 2 * xx[2] ) * std::cos( 8 * xx[1] ) +
                           512 * std::cos( 8 * xx[0] );
                 },
                 []( const Point3D& xx ) -> real_t {
                    return 2 * std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::cos( 2 * xx[2] ) - 8 * std::cos( 2 * xx[1] );
                 },
                 []( const Point3D& ) -> real_t { return 0; } ),

             std::make_shared< EGP0StokesOperator >( storage, minLevel, maxLevel ),
             storage,
             minLevel,
             maxLevel + 1,
             2,
             true );
      }
   }

   // pyramid
   if ( false )
   {
      auto                         meshInfo = hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_2el.msh" );
      hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                                 walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
      auto storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, 1 );

      // pyramid_2el, inhom. split solution
      WALBERLA_LOG_INFO_ON_ROOT( "### pyramid_2el, inhom. split solution ###" );
      StokesConvergenceOrderTest< EGP0StokesOperator >(
          "EGP0StokesOp3D_pyramid_2el_inhom_split",
          std::make_tuple( []( const Point3D& xx ) -> real_t { return -real_c( 4 ) * std::cos( real_c( 4 ) * xx[2] ); },
                           []( const Point3D& xx ) -> real_t { return real_c( 8 ) * std::cos( real_c( 8 ) * xx[0] ); },
                           []( const Point3D& xx ) -> real_t { return -real_c( 2 ) * std::cos( real_c( 2 ) * xx[1] ); },
                           []( const Point3D& xx ) -> real_t {
                              return std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] );
                           } ),
          std::make_tuple(
              []( const Point3D& xx ) -> real_t {
                 return 4 * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] ) * std::cos( 4 * xx[0] ) - 64 * std::cos( 4 * xx[2] );
              },
              []( const Point3D& xx ) -> real_t {
                 return 8 * std::sin( 4 * xx[0] ) * std::sin( 2 * xx[2] ) * std::cos( 8 * xx[1] ) + 512 * std::cos( 8 * xx[0] );
              },
              []( const Point3D& xx ) -> real_t {
                 return 2 * std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::cos( 2 * xx[2] ) - 8 * std::cos( 2 * xx[1] );
              },
              []( const Point3D& ) -> real_t { return 0; } ),

          std::make_shared< EGP0StokesOperator >( storage, minLevel, maxLevel ),
          storage,
          minLevel,
          maxLevel,
          2,
          false,
          false );

      StokesConvergenceOrderTest< hyteg::P2P1TaylorHoodStokesOperator >(
          "P2P1StokesOperator3D_pyramid_2el_inhom_split",
          std::make_tuple(
              []( const hyteg::Point3D& xx ) { return -real_c( 4 ) * std::cos( real_c( 4 ) * xx[2] ); },
              []( const hyteg::Point3D& xx ) { return real_c( 8 ) * std::cos( real_c( 8 ) * xx[0] ); },
              []( const hyteg::Point3D& xx ) { return -real_c( 2 ) * std::cos( real_c( 2 ) * xx[1] ); },
              []( const hyteg::Point3D& xx ) { return std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] ); } ),
          std::make_tuple(
              []( const hyteg::Point3D& xx ) {
                 return 4 * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] ) * std::cos( 4 * xx[0] ) - 64 * std::cos( 4 * xx[2] );
              },
              []( const hyteg::Point3D& xx ) {
                 return 8 * std::sin( 4 * xx[0] ) * std::sin( 2 * xx[2] ) * std::cos( 8 * xx[1] ) + 512 * std::cos( 8 * xx[0] );
              },
              []( const hyteg::Point3D& xx ) {
                 return 2 * std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::cos( 2 * xx[2] ) - 8 * std::cos( 2 * xx[1] );
              },
              []( const hyteg::Point3D& ) { return 0; } ),
          std::make_shared< hyteg::P2P1TaylorHoodStokesOperator >( storage, minLevel, maxLevel ),
          storage,
          minLevel,
          maxLevel,
          5,
          false,
          false );
   }

   // cube_6el, simple
   if ( false )
   {
      MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid( Point3D( { -1, -1, -1 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );
      //auto meshInfo = hyteg::MeshInfo::fromGmshFile("../../data/meshes/3D/cube_6el.msh");

      hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                                 walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
      auto storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, 1 );

      WALBERLA_LOG_INFO_ON_ROOT( "### cube_6el, simple ###" );
      StokesConvergenceOrderTest< EGP0StokesOperator >( "EGP0StokesOp3D_cube_6el_simple",
                                                        std::make_tuple( []( const Point3D& xx ) -> real_t { return 2 * xx[0]; },
                                                                         []( const Point3D& xx ) -> real_t { return -3 * xx[1]; },
                                                                         []( const Point3D& xx ) -> real_t { return 1 * xx[2]; },
                                                                         []( const Point3D& ) -> real_t { return 1; } ),
                                                        std::make_tuple( []( const Point3D& ) -> real_t { return 0; },
                                                                         []( const Point3D& ) -> real_t { return 0; },
                                                                         []( const Point3D& ) -> real_t { return 0; },
                                                                         []( const Point3D& ) -> real_t { return 0; } ),

                                                        std::make_shared< EGP0StokesOperator >( storage, minLevel, maxLevel ),
                                                        storage,
                                                        minLevel,
                                                        maxLevel,
                                                        5,
                                                        true );
   }

   // cube_6el, hom. solution
   if ( false )
   {
      // MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid( Point3D( { 0, 0, 0 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );
      auto meshInfo = hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_6el.msh" );

      hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                                 walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
      auto storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, 1 );

      // cube, hom. solution
      WALBERLA_LOG_INFO_ON_ROOT( "### cube_6el, hom. solution ###" );
      StokesConvergenceOrderTest< EGP0StokesOperator >(
          "EGP0StokesOp3D_cube_6el_hom",
          std::make_tuple(
              []( const Point3D& xx ) -> real_t {
                 const real_t x  = xx[0];
                 const real_t y  = xx[1];
                 const real_t z  = xx[2];
                 const real_t x0 = 2 * M_PI;
                 return std::sin( x * x0 ) * std::sin( x0 * y ) * std::sin( x0 * z );
              },
              []( const Point3D& xx ) -> real_t {
                 const real_t x  = xx[0];
                 const real_t y  = xx[1];
                 const real_t z  = xx[2];
                 const real_t x0 = 2 * M_PI;
                 return std::sin( x * x0 ) * std::sin( x0 * y ) * std::sin( x0 * z );
              },
              []( const Point3D& xx ) -> real_t {
                 const real_t x  = xx[0];
                 const real_t y  = xx[1];
                 const real_t z  = xx[2];
                 const real_t x0 = 2 * M_PI;
                 return std::sin( x * x0 ) * std::sin( x0 * y ) * std::sin( x0 * z );
              },
              []( const Point3D& xx ) -> real_t {
                 return std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] );
              } ),
          std::make_tuple(
              []( const Point3D& xx ) -> real_t {
                 const real_t x  = xx[0];
                 const real_t y  = xx[1];
                 const real_t z  = xx[2];
                 const real_t x0 = 2 * z;
                 const real_t x1 = 2 * M_PI;
                 return 4 * std::sin( x0 ) * std::sin( 8 * y ) * std::cos( 4 * x ) +
                        12 * std::pow( M_PI, 2 ) * std::sin( M_PI * x0 ) * std::sin( x * x1 ) * std::sin( x1 * y );
              },
              []( const Point3D& xx ) -> real_t {
                 const real_t x  = xx[0];
                 const real_t y  = xx[1];
                 const real_t z  = xx[2];
                 const real_t x0 = 2 * z;
                 const real_t x1 = 2 * M_PI;
                 return 8 * std::sin( 4 * x ) * std::sin( x0 ) * std::cos( 8 * y ) +
                        12 * std::pow( M_PI, 2 ) * std::sin( M_PI * x0 ) * std::sin( x * x1 ) * std::sin( x1 * y );
              },
              []( const Point3D& xx ) -> real_t {
                 const real_t x  = xx[0];
                 const real_t y  = xx[1];
                 const real_t z  = xx[2];
                 const real_t x0 = 2 * z;
                 const real_t x1 = 2 * M_PI;
                 return 2 * std::sin( 4 * x ) * std::sin( 8 * y ) * std::cos( x0 ) +
                        12 * std::pow( M_PI, 2 ) * std::sin( M_PI * x0 ) * std::sin( x * x1 ) * std::sin( x1 * y );
              },
              []( const Point3D& xx ) -> real_t {
                 const real_t x  = xx[0];
                 const real_t y  = xx[1];
                 const real_t z  = xx[2];
                 const real_t x0 = 2 * M_PI;
                 const real_t x1 = x0 * y;
                 const real_t x2 = std::sin( x1 );
                 const real_t x3 = x * x0;
                 const real_t x4 = std::sin( x3 );
                 const real_t x5 = x0 * z;
                 const real_t x6 = x0 * std::sin( x5 );
                 return -x0 * x2 * x4 * std::cos( x5 ) - x2 * x6 * std::cos( x3 ) - x4 * x6 * std::cos( x1 );
              } ),

          std::make_shared< EGP0StokesOperator >( storage, minLevel, maxLevel ),
          storage,
          minLevel,
          maxLevel,
          2,
          true );
   }

   // cube_6el, inhom.
   if ( false )
   {
      //MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid(Point3D({0, 0, 0}), Point3D({1, 1, 1}), 4, 4, 4);
      auto                         meshInfo = hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_6el.msh" );
      hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                                 walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
      auto storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, 1 );

      // cube_6el, inhom.
      WALBERLA_LOG_INFO_ON_ROOT( "### cube_6el, inhom. ###" );
      StokesConvergenceOrderTest< EGP0StokesOperator >(
          "EGP0StokesOp3D_cube_6el_inhom",
          std::make_tuple(
              []( const hyteg::Point3D& xx ) { return -real_c( 4 ) * std::cos( real_c( 4 ) * xx[2] ); },
              []( const hyteg::Point3D& xx ) { return real_c( 8 ) * std::cos( real_c( 8 ) * xx[0] ); },
              []( const hyteg::Point3D& xx ) { return -real_c( 2 ) * std::cos( real_c( 2 ) * xx[1] ); },
              []( const hyteg::Point3D& xx ) { return std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] ); } ),
          std::make_tuple(
              []( const hyteg::Point3D& xx ) {
                 return 4 * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] ) * std::cos( 4 * xx[0] ) - 64 * std::cos( 4 * xx[2] );
              },
              []( const hyteg::Point3D& xx ) {
                 return 8 * std::sin( 4 * xx[0] ) * std::sin( 2 * xx[2] ) * std::cos( 8 * xx[1] ) + 512 * std::cos( 8 * xx[0] );
              },
              []( const hyteg::Point3D& xx ) {
                 return 2 * std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::cos( 2 * xx[2] ) - 8 * std::cos( 2 * xx[1] );
              },
              []( const hyteg::Point3D& ) { return 0; } ),
          std::make_shared< EGP0StokesOperator >( storage, minLevel, maxLevel ),
          storage,
          minLevel,
          maxLevel,
          5,
          false,
          false );
   }

   // cube_6el, inhom., Nitsche BCs and rhs-integration
   if ( true )
   {
      //MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid(Point3D({0, 0, 0}), Point3D({1, 1, 1}), 4, 4, 4);
      auto                         meshInfo = hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_6el.msh" );
      hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                                 walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
      auto storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, 1 );
       auto resNormsP2P1 = {1e-5, 1e-6, 1e-6, 1e-5, 1e-6};

      // EG
      if ( true )
      {
         EGP0StokesOperatorNitscheBC EGP0StokesOp( storage, minLevel, maxLevel );
         WALBERLA_LOG_INFO_ON_ROOT( "### EG cube_6el, inhom.,  Nitsche BCs and rhs-integration ###" );
         StokesConvergenceOrderTest< EGP0StokesOperatorNitscheBC >(
             "EGP0StokesOperatorNitscheBC3D_cube_6el_inhom",
             std::make_tuple( []( const hyteg::Point3D& xx ) { return -real_c( 4 ) * std::cos( real_c( 4 ) * xx[2] ); },
                              []( const hyteg::Point3D& xx ) { return real_c( 8 ) * std::cos( real_c( 8 ) * xx[0] ); },
                              []( const hyteg::Point3D& xx ) { return -real_c( 2 ) * std::cos( real_c( 2 ) * xx[1] ); },
                              []( const hyteg::Point3D& xx ) {
                                 return std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] );
                              } ),
             std::make_tuple(
                 []( const hyteg::Point3D& xx ) {
                    return 4 * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] ) * std::cos( 4 * xx[0] ) - 64 * std::cos( 4 * xx[2] );
                 },
                 []( const hyteg::Point3D& xx ) {
                    return 8 * std::sin( 4 * xx[0] ) * std::sin( 2 * xx[2] ) * std::cos( 8 * xx[1] ) +
                           512 * std::cos( 8 * xx[0] );
                 },
                 []( const hyteg::Point3D& xx ) {
                    return 2 * std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::cos( 2 * xx[2] ) - 8 * std::cos( 2 * xx[1] );
                 },
                 []( const hyteg::Point3D& ) { return 0; } ),
             std::make_shared< EGP0StokesOperatorNitscheBC >( storage, minLevel, maxLevel ),
             storage,
             minLevel,
             maxLevel,
             2,
             false,
             false, NULL, std::make_shared<std::vector<real_t>>(
                         resNormsP2P1) );
      }

      // P1P0
      if constexpr ( false )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "### P1P0 cube_6el, inhom. ###" );
         StokesConvergenceOrderTest< hyteg::P1P0StokesOperator >(
             "P1P0StokesOperator3D_cube_6el_inhom",
             std::make_tuple( []( const hyteg::Point3D& xx ) { return -real_c( 4 ) * std::cos( real_c( 4 ) * xx[2] ); },
                              []( const hyteg::Point3D& xx ) { return real_c( 8 ) * std::cos( real_c( 8 ) * xx[0] ); },
                              []( const hyteg::Point3D& xx ) { return -real_c( 2 ) * std::cos( real_c( 2 ) * xx[1] ); },
                              []( const hyteg::Point3D& xx ) {
                                 return std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] );
                              } ),
             std::make_tuple(
                 []( const hyteg::Point3D& xx ) {
                    return 4 * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] ) * std::cos( 4 * xx[0] ) - 64 * std::cos( 4 * xx[2] );
                 },
                 []( const hyteg::Point3D& xx ) {
                    return 8 * std::sin( 4 * xx[0] ) * std::sin( 2 * xx[2] ) * std::cos( 8 * xx[1] ) +
                           512 * std::cos( 8 * xx[0] );
                 },
                 []( const hyteg::Point3D& xx ) {
                    return 2 * std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::cos( 2 * xx[2] ) - 8 * std::cos( 2 * xx[1] );
                 },
                 []( const hyteg::Point3D& ) { return 0; } ),
             std::make_shared< hyteg::P1P0StokesOperator >( storage, minLevel, maxLevel, 0.1 ),
             storage,
             minLevel,
             maxLevel,
             5,
             false,
             true );
      }

      // P2P1
      if ( true )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "### P2P1 cube_6el, inhom. ###" );
         StokesConvergenceOrderTest< hyteg::P2P1TaylorHoodStokesOperator >(
             "P2P1StokesOperator3D_cube_6el_inhom",
             std::make_tuple( []( const hyteg::Point3D& xx ) { return -real_c( 4 ) * std::cos( real_c( 4 ) * xx[2] ); },
                              []( const hyteg::Point3D& xx ) { return real_c( 8 ) * std::cos( real_c( 8 ) * xx[0] ); },
                              []( const hyteg::Point3D& xx ) { return -real_c( 2 ) * std::cos( real_c( 2 ) * xx[1] ); },
                              []( const hyteg::Point3D& xx ) {
                                 return std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] );
                              } ),
             std::make_tuple(
                 []( const hyteg::Point3D& xx ) {
                    return 4 * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] ) * std::cos( 4 * xx[0] ) - 64 * std::cos( 4 * xx[2] );
                 },
                 []( const hyteg::Point3D& xx ) {
                    return 8 * std::sin( 4 * xx[0] ) * std::sin( 2 * xx[2] ) * std::cos( 8 * xx[1] ) +
                           512 * std::cos( 8 * xx[0] );
                 },
                 []( const hyteg::Point3D& xx ) {
                    return 2 * std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::cos( 2 * xx[2] ) - 8 * std::cos( 2 * xx[1] );
                 },
                 []( const hyteg::Point3D& ) { return 0; } ),
             std::make_shared< hyteg::P2P1TaylorHoodStokesOperator >( storage, minLevel, maxLevel ),
             storage,
             minLevel,
             maxLevel,
             2,
             false,
             false, NULL, std::make_shared<std::vector<real_t>>(
                         resNormsP2P1) );
      }
   }

   // cube_24el, hom. solution
   if ( false )
   {
      // MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid( Point3D( { 0, 0, 0 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );
      auto meshInfo = hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_24el.msh" );

      hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                                 walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
      auto storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, 1 );

      // cube, hom. solution
      WALBERLA_LOG_INFO_ON_ROOT( "### cube_24el, hom. solution ###" );
      StokesConvergenceOrderTest< EGP0StokesOperator >(
          "EGP0StokesOp3D_cube_24el_hom",
          std::make_tuple(
              []( const Point3D& xx ) -> real_t {
                 const real_t x  = xx[0];
                 const real_t y  = xx[1];
                 const real_t z  = xx[2];
                 const real_t x0 = 2 * M_PI;
                 return std::sin( x * x0 ) * std::sin( x0 * y ) * std::sin( x0 * z );
              },
              []( const Point3D& xx ) -> real_t {
                 const real_t x  = xx[0];
                 const real_t y  = xx[1];
                 const real_t z  = xx[2];
                 const real_t x0 = 2 * M_PI;
                 return std::sin( x * x0 ) * std::sin( x0 * y ) * std::sin( x0 * z );
              },
              []( const Point3D& xx ) -> real_t {
                 const real_t x  = xx[0];
                 const real_t y  = xx[1];
                 const real_t z  = xx[2];
                 const real_t x0 = 2 * M_PI;
                 return std::sin( x * x0 ) * std::sin( x0 * y ) * std::sin( x0 * z );
              },
              []( const Point3D& xx ) -> real_t {
                 return std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] );
              } ),
          std::make_tuple(
              []( const Point3D& xx ) -> real_t {
                 const real_t x  = xx[0];
                 const real_t y  = xx[1];
                 const real_t z  = xx[2];
                 const real_t x0 = 2 * z;
                 const real_t x1 = 2 * M_PI;
                 return 4 * std::sin( x0 ) * std::sin( 8 * y ) * std::cos( 4 * x ) +
                        12 * std::pow( M_PI, 2 ) * std::sin( M_PI * x0 ) * std::sin( x * x1 ) * std::sin( x1 * y );
              },
              []( const Point3D& xx ) -> real_t {
                 const real_t x  = xx[0];
                 const real_t y  = xx[1];
                 const real_t z  = xx[2];
                 const real_t x0 = 2 * z;
                 const real_t x1 = 2 * M_PI;
                 return 8 * std::sin( 4 * x ) * std::sin( x0 ) * std::cos( 8 * y ) +
                        12 * std::pow( M_PI, 2 ) * std::sin( M_PI * x0 ) * std::sin( x * x1 ) * std::sin( x1 * y );
              },
              []( const Point3D& xx ) -> real_t {
                 const real_t x  = xx[0];
                 const real_t y  = xx[1];
                 const real_t z  = xx[2];
                 const real_t x0 = 2 * z;
                 const real_t x1 = 2 * M_PI;
                 return 2 * std::sin( 4 * x ) * std::sin( 8 * y ) * std::cos( x0 ) +
                        12 * std::pow( M_PI, 2 ) * std::sin( M_PI * x0 ) * std::sin( x * x1 ) * std::sin( x1 * y );
              },
              []( const Point3D& xx ) -> real_t {
                 const real_t x  = xx[0];
                 const real_t y  = xx[1];
                 const real_t z  = xx[2];
                 const real_t x0 = 2 * M_PI;
                 const real_t x1 = x0 * y;
                 const real_t x2 = std::sin( x1 );
                 const real_t x3 = x * x0;
                 const real_t x4 = std::sin( x3 );
                 const real_t x5 = x0 * z;
                 const real_t x6 = x0 * std::sin( x5 );
                 return -x0 * x2 * x4 * std::cos( x5 ) - x2 * x6 * std::cos( x3 ) - x4 * x6 * std::cos( x1 );
              } ),
          std::make_shared< EGP0StokesOperator >( storage, minLevel, maxLevel ),
          storage,
          minLevel,
          maxLevel - 1,
          2,
          true );
   }

   // IIPG: cube_24el, hom. solution
   if ( false )
   {
      // MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid( Point3D( { 0, 0, 0 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );
      auto meshInfo = hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_24el.msh" );

      hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                                 walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
      auto storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, 1 );

      EGP0IIPGStokesOperator EGP0StokesOp( storage, minLevel - 1, maxLevel );

      // cube, hom. solution
      WALBERLA_LOG_INFO_ON_ROOT( "### cube_24el, hom. solution ###" );
      StokesConvergenceOrderTest< EGP0IIPGStokesOperator >(
          "EGP0IIPGStokesOp3D_cube_24el_hom",
          std::make_tuple(
              []( const Point3D& xx ) -> real_t {
                 const real_t x  = xx[0];
                 const real_t y  = xx[1];
                 const real_t z  = xx[2];
                 const real_t x0 = 2 * M_PI;
                 return std::sin( x * x0 ) * std::sin( x0 * y ) * std::sin( x0 * z );
              },
              []( const Point3D& xx ) -> real_t {
                 const real_t x  = xx[0];
                 const real_t y  = xx[1];
                 const real_t z  = xx[2];
                 const real_t x0 = 2 * M_PI;
                 return std::sin( x * x0 ) * std::sin( x0 * y ) * std::sin( x0 * z );
              },
              []( const Point3D& xx ) -> real_t {
                 const real_t x  = xx[0];
                 const real_t y  = xx[1];
                 const real_t z  = xx[2];
                 const real_t x0 = 2 * M_PI;
                 return std::sin( x * x0 ) * std::sin( x0 * y ) * std::sin( x0 * z );
              },
              []( const Point3D& xx ) -> real_t {
                 return std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] );
              } ),
          std::make_tuple(
              []( const Point3D& xx ) -> real_t {
                 const real_t x  = xx[0];
                 const real_t y  = xx[1];
                 const real_t z  = xx[2];
                 const real_t x0 = 2 * z;
                 const real_t x1 = 2 * M_PI;
                 return 4 * std::sin( x0 ) * std::sin( 8 * y ) * std::cos( 4 * x ) +
                        12 * std::pow( M_PI, 2 ) * std::sin( M_PI * x0 ) * std::sin( x * x1 ) * std::sin( x1 * y );
              },
              []( const Point3D& xx ) -> real_t {
                 const real_t x  = xx[0];
                 const real_t y  = xx[1];
                 const real_t z  = xx[2];
                 const real_t x0 = 2 * z;
                 const real_t x1 = 2 * M_PI;
                 return 8 * std::sin( 4 * x ) * std::sin( x0 ) * std::cos( 8 * y ) +
                        12 * std::pow( M_PI, 2 ) * std::sin( M_PI * x0 ) * std::sin( x * x1 ) * std::sin( x1 * y );
              },
              []( const Point3D& xx ) -> real_t {
                 const real_t x  = xx[0];
                 const real_t y  = xx[1];
                 const real_t z  = xx[2];
                 const real_t x0 = 2 * z;
                 const real_t x1 = 2 * M_PI;
                 return 2 * std::sin( 4 * x ) * std::sin( 8 * y ) * std::cos( x0 ) +
                        12 * std::pow( M_PI, 2 ) * std::sin( M_PI * x0 ) * std::sin( x * x1 ) * std::sin( x1 * y );
              },
              []( const Point3D& xx ) -> real_t {
                 const real_t x  = xx[0];
                 const real_t y  = xx[1];
                 const real_t z  = xx[2];
                 const real_t x0 = 2 * M_PI;
                 const real_t x1 = x0 * y;
                 const real_t x2 = std::sin( x1 );
                 const real_t x3 = x * x0;
                 const real_t x4 = std::sin( x3 );
                 const real_t x5 = x0 * z;
                 const real_t x6 = x0 * std::sin( x5 );
                 return -x0 * x2 * x4 * std::cos( x5 ) - x2 * x6 * std::cos( x3 ) - x4 * x6 * std::cos( x1 );
              } ),
          std::make_shared< EGP0IIPGStokesOperator >( storage, minLevel, maxLevel ),

          storage,
          minLevel,
          maxLevel - 1,
          2,
          true );
   }

   // cube_24el, inhom. solution
   if ( false )
   {
      // MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid( Point3D( { 0, 0, 0 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );
      auto meshInfo = hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_24el.msh" );

      hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                                 walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
      auto storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, 1 );

      // cube, hom. solution
      WALBERLA_LOG_INFO_ON_ROOT( "### cube_24el, inhom. solution ###" );
      StokesConvergenceOrderTest< EGP0StokesOperator >(
          "EGP0StokesOp3D_cube_24el_inhom",
          std::make_tuple(
              []( const hyteg::Point3D& xx ) { return -real_c( 4 ) * std::cos( real_c( 4 ) * xx[2] ); },
              []( const hyteg::Point3D& xx ) { return real_c( 8 ) * std::cos( real_c( 8 ) * xx[0] ); },
              []( const hyteg::Point3D& xx ) { return -real_c( 2 ) * std::cos( real_c( 2 ) * xx[1] ); },
              []( const hyteg::Point3D& xx ) { return std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] ); } ),
          std::make_tuple(
              []( const hyteg::Point3D& xx ) {
                 return 4 * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] ) * std::cos( 4 * xx[0] ) - 64 * std::cos( 4 * xx[2] );
              },
              []( const hyteg::Point3D& xx ) {
                 return 8 * std::sin( 4 * xx[0] ) * std::sin( 2 * xx[2] ) * std::cos( 8 * xx[1] ) + 512 * std::cos( 8 * xx[0] );
              },
              []( const hyteg::Point3D& xx ) {
                 return 2 * std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::cos( 2 * xx[2] ) - 8 * std::cos( 2 * xx[1] );
              },
              []( const hyteg::Point3D& ) { return 0; } ),
          std::make_shared< EGP0StokesOperator >( storage, minLevel, maxLevel ),
          storage,
          minLevel,
          maxLevel - 1,
          2,
          true );
   }

   // IIPG: cube_24el, inhom. solution
   if ( false )
   {
      // MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid( Point3D( { 0, 0, 0 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );
      auto meshInfo = hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_24el.msh" );

      hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                                 walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
      auto storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, 1 );

      EGP0IIPGStokesOperator EGP0StokesOp( storage, minLevel - 1, maxLevel );

      // cube, hom. solution
      WALBERLA_LOG_INFO_ON_ROOT( "### cube_24el, inhom. solution ###" );
      StokesConvergenceOrderTest< EGP0IIPGStokesOperator >(
          "EGP0IIPGStokesOp3D_cube_24el_inhom",
          std::make_tuple(
              []( const hyteg::Point3D& xx ) { return -real_c( 4 ) * std::cos( real_c( 4 ) * xx[2] ); },
              []( const hyteg::Point3D& xx ) { return real_c( 8 ) * std::cos( real_c( 8 ) * xx[0] ); },
              []( const hyteg::Point3D& xx ) { return -real_c( 2 ) * std::cos( real_c( 2 ) * xx[1] ); },
              []( const hyteg::Point3D& xx ) { return std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] ); } ),
          std::make_tuple(
              []( const hyteg::Point3D& xx ) {
                 return 4 * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] ) * std::cos( 4 * xx[0] ) - 64 * std::cos( 4 * xx[2] );
              },
              []( const hyteg::Point3D& xx ) {
                 return 8 * std::sin( 4 * xx[0] ) * std::sin( 2 * xx[2] ) * std::cos( 8 * xx[1] ) + 512 * std::cos( 8 * xx[0] );
              },
              []( const hyteg::Point3D& xx ) {
                 return 2 * std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::cos( 2 * xx[2] ) - 8 * std::cos( 2 * xx[1] );
              },
              []( const hyteg::Point3D& ) { return 0; } ),
          std::make_shared< EGP0IIPGStokesOperator >( storage, minLevel, maxLevel ),
          storage,
          minLevel,
          maxLevel,
          2,
          true );
   }

   /*
     StokesConvergenceOrderTest< EGP0StokesOperator >(
         "EGP0StokesOp3D",
         std::make_tuple( []( const Point3D& xx ) -> real_t { return -real_c( 4 ) * std::cos( real_c( 4 ) * xx[2] ); },
                          []( const Point3D& xx ) -> real_t { return real_c( 8 ) * std::cos( real_c( 8 ) * xx[0] ); },
                          []( const Point3D& xx ) -> real_t { return -real_c( 2 ) * std::cos( real_c( 2 ) * xx[1] ); },
                          []( const Point3D& xx ) -> real_t {
                             return std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] );
                          } ),
         std::make_tuple(
             []( const Point3D& xx ) -> real_t {
                return 4 * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] ) * std::cos( 4 * xx[0] ) - 64 * std::cos( 4 * xx[2] );
             },
             []( const Point3D& xx ) -> real_t {
                return 8 * std::sin( 4 * xx[0] ) * std::sin( 2 * xx[2] ) * std::cos( 8 * xx[1] ) + 512 * std::cos( 8 * xx[0] );
             },
             []( const Point3D& xx ) -> real_t {
                return 2 * std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::cos( 2 * xx[2] ) - 8 * std::cos( 2 * xx[1] );
             },
             []( const Point3D& ) -> real_t { return 0; } ),
         EGP0StokesOp,
         storage,
         minLevel,
         maxLevel,
         2,
         true );

     StokesConvergenceOrderTest< hyteg::P2P1TaylorHoodStokesOperator >(
         "P2P1StokesOp3D",
         std::make_tuple( []( const Point3D& xx ) -> real_t { return -real_c( 4 ) * std::cos( real_c( 4 ) * xx[2] ); },
                          []( const Point3D& xx ) -> real_t { return real_c( 8 ) * std::cos( real_c( 8 ) * xx[0] ); },
                          []( const Point3D& xx ) -> real_t { return -real_c( 2 ) * std::cos( real_c( 2 ) * xx[1] ); },
                          []( const Point3D& xx ) -> real_t {
                             return std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] );
                          } ),
         std::make_tuple(
             []( const Point3D& xx ) -> real_t {
                return 4 * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] ) * std::cos( 4 * xx[0] ) - 64 * std::cos( 4 * xx[2] );
             },
             []( const Point3D& xx ) -> real_t {
                return 8 * std::sin( 4 * xx[0] ) * std::sin( 2 * xx[2] ) * std::cos( 8 * xx[1] ) + 512 * std::cos( 8 * xx[0] );
             },
             []( const Point3D& xx ) -> real_t {
                return 2 * std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::cos( 2 * xx[2] ) - 8 * std::cos( 2 * xx[1] );
             },
             []( const Point3D& ) -> real_t { return 0; } ),
         P2P1StokesOp,
         storage,
         minLevel,
         maxLevel,
         2 );*/
}

} // namespace eg
} // namespace dg
} // namespace hyteg
