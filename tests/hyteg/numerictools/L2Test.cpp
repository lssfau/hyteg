/*
 * Copyright (c) 2017-2023 Benjamin Mann.
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
#include "core/math/all.h"

#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/numerictools/L2Space.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VariableOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

using walberla::real_t;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

template < class MassOperator, class FE >
void testL2Dot( const std::string&                               testName,
                const SetupPrimitiveStorage&                     setupStorage,
                const std::function< real_t( const Point3D& ) >& u,
                const std::function< real_t( const Point3D& ) >& v,
                const std::function< real_t( const Point3D& ) >& f,
                const real_t&                                    uv_expected,
                const uint_t&                                    lvl,
                const uint_t&                                    quad,
                const real_t&                                    accuracy_uv,
                const real_t&                                    accuracy_rhs )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Test: " << testName );

   // storage
   auto              storage = std::make_shared< PrimitiveStorage >( setupStorage );
   L2Space< real_t > L2( storage, lvl, quad );

   // operators and functions
   MassOperator M( storage, lvl, lvl );
   FE           f_proj( "f", storage, lvl, lvl );
   FE           b_M( "b_M", storage, lvl, lvl );
   FE           b_q( "b_q", storage, lvl, lvl );
   FE           e( "b_q - b_M", storage, lvl, lvl );

   // compute (u,v)_0
   auto uv = L2.dot( u, v );
   WALBERLA_CHECK_LESS( std::abs( uv - uv_expected ), accuracy_uv );

   // compare (φ_i, f)_i with [Mf]_i
   f_proj.interpolate( f, lvl, All );
   M.apply( f_proj, b_M, lvl, All );
   L2.dot( f, b_q );
   e.assign( { 1.0, -1.0 }, { b_q, b_M }, lvl, All );
   auto err = std::sqrt( e.dotGlobal( e, lvl ) );
   WALBERLA_CHECK_LESS( err, accuracy_rhs );
}

int main( int argc, char* argv[] )
{
   // mpi
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   // mesh
   auto                  meshInfo = MeshInfo::emptyMeshInfo();
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   // functions
   std::function< real_t( const Point3D& ) > u, v, f;

   // -------
   //  TESTS
   // -------

   const real_t eps = ( std::is_same_v< real_t, double > ) ? 1e-16 : 1e-7;

   // todo implement tests

   // test annulus
   meshInfo     = MeshInfo::meshAnnulus( 2, 4, 6, 2 );
   setupStorage = SetupPrimitiveStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   AnnulusMap::setMap( setupStorage );
   u = [&]( const Point3D& x ) { return x.norm(); };                             // u(r,φ) = r
   v = [&]( const Point3D& x ) { return std::sin( std::atan2( x[1], x[0] ) ); }; // u(r,φ) = sin(φ)
   f = [&]( const Point3D& x ) { return u( x ) * v( x ); };                      // f(x) = u(x)v(x)
   testL2Dot( "Annulus1", setupStorage, u, v, f, 0.0, 4, 5, eps, 1e2 * eps );
   testL2Dot( "Annulus2", setupStorage, f, f, f, 60.0 * pi, 4, 5, 2e2 * eps, 1e2 * eps );

   // meshInfo = MeshInfo::meshRectangle( Point2D( 0.0, 0.0 ), Point2D( 1.0, 1.0 ), MeshInfo::CROSS, 1, 1 );

   // Point3D                                   c{ 1.0, 2.0, 3.0 };
   // std::function< real_t( const Point3D& ) > linear = [&]( const Point3D& x ) { return c.dot( x ); };
   // Matrix3r                                  C;

   return EXIT_SUCCESS;
}
