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

template < typename ValueType >
void testL2Dot( const std::shared_ptr< PrimitiveStorage >&       storage,
                const std::function< real_t( const Point3D& ) >& u,
                const std::function< real_t( const Point3D& ) >& v,
                const real_t&                                    uv_expected,
                const uint_t&                                    lvl,
                const uint_t&                                    quad,
                const real_t&                                    accuracy )
{
   L2Space< Undefined, ValueType > L2( storage, lvl, quad );

   Undefined b;
   L2.dot(u,b);

   // compute (u,v)_0
   auto uv  = L2.dot( u, v );
   auto err = std::abs( uv - uv_expected );
   WALBERLA_LOG_INFO_ON_ROOT( "|(u,v)_q - (u,v)_L2| = " << err );
   WALBERLA_CHECK_LESS( err, accuracy );
}

template < class MassOperator, class FE, typename ValueType >
void testRHS( const std::shared_ptr< PrimitiveStorage >&       storage,
              const std::function< real_t( const Point3D& ) >& f,
              const uint_t&                                    lvl,
              const uint_t&                                    quad,
              const real_t&                                    accuracy )
{
   L2Space< FE, ValueType > L2( storage, lvl, quad );

   // operators and functions
   MassOperator M( storage, lvl, lvl );
   FE           f_proj( "f", storage, lvl, lvl );
   FE           b_M( "b_M", storage, lvl, lvl );
   FE           b_q( "b_q", storage, lvl, lvl );
   FE           e( "b_q - b_M", storage, lvl, lvl );

   // number of global DoFs
   auto n_dof = numberOfGlobalDoFs< typename FE::Tag >( *storage, lvl );

   // compare (φ_i, f)_i with [Mf]_i
   f_proj.interpolate( f, lvl, All );
   M.apply( f_proj, b_M, lvl, All ); // compute rhs using mass matrix
   L2.dot( f, b_q );                 // compute rhs using quadrature
   e.assign( { 1.0, -1.0 }, { b_q, b_M }, lvl, All );
   auto err = std::sqrt( e.dotGlobal( e, lvl ) / real_c( n_dof ) ); // weighted error
   WALBERLA_LOG_INFO_ON_ROOT( "||(v_i, f)_q - [Mf]_i|| = " << err );
   WALBERLA_CHECK_LESS( err, real_c( n_dof ) * accuracy );
}

int main( int argc, char* argv[] )
{
   // mpi
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   // mesh
   auto                                meshInfo = MeshInfo::emptyMeshInfo();
   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage;

   // functions
   std::function< real_t( const Point3D& ) > u, v, f;

   // -------
   //  TESTS
   // -------

   const real_t eps = ( std::is_same_v< real_t, double > ) ? 5e-12 : 5e-5;

   using P1     = P1Function< real_t >;
   using P1Mass = P1BlendingMassOperator;

   uint_t quad = 5;
   uint_t lvl  = 5;

   // === test annulus ===

   meshInfo     = MeshInfo::meshAnnulus( 2, 4, MeshInfo::CRISS, 6, 2 );
   setupStorage = SetupPrimitiveStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   AnnulusMap::setMap( setupStorage );
   storage = std::make_shared< PrimitiveStorage >( setupStorage );

   u = [&]( const Point3D& x ) { return x.dot( x ); }; // u(r,φ) = r^2
   v = [&]( const Point3D& x ) {
      auto s = std::sin( std::atan2( x[1], x[0] ) );
      return s * s;
   };                                                       // v(r,φ) = sin^2(φ)
   f = [&]( const Point3D& x ) { return u( x ) * v( x ); }; // f(x) = u(x)v(x)

   WALBERLA_LOG_INFO_ON_ROOT( "Test: Annulus" );
   testL2Dot< real_t >( storage, u, v, 60.0 * pi, lvl, quad, eps );
   testRHS< P1Mass, P1, real_t >( storage, f, lvl, quad, 1e-5 );

   // === test polynomials 3D ===

   meshInfo     = MeshInfo::meshCuboid( { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 }, 1, 1, 1 );
   setupStorage = SetupPrimitiveStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   storage      = std::make_shared< PrimitiveStorage >( setupStorage );

   Point3D one{ 1.0, 1.0, 1.0 };
   Point3D c{ 1.0, 2.0, 3.0 };
   u = [&]( const Point3D& x ) { return one.dot( x ); }; // u(x) = x+y+z
   v = [&]( const Point3D& x ) { return c.dot( x ); };   // v(x) = c⋅x

   WALBERLA_LOG_INFO_ON_ROOT( "Test: Polynomial 3D" );
   testL2Dot< real_t >( storage, u, v, 5.0, lvl, quad, eps );
   testRHS< P1Mass, P1, real_t >( storage, v, lvl, quad, eps );

   // todo implement tests for other discretizations

   return EXIT_SUCCESS;
}
