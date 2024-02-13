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

#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
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

template < typename ValueType, uint_t Quad >
void testL2Dot( const std::shared_ptr< PrimitiveStorage >&          storage,
                const std::function< ValueType( const Point3D& ) >& u,
                const std::function< ValueType( const Point3D& ) >& v,
                const real_t&                                       uv_expected,
                const uint_t&                                       lvl,
                const real_t&                                       accuracy )
{
   L2Space< Quad, Undefined, ValueType > L2( storage, lvl );

   // compute (u,v)_0
   auto uv  = L2.dot( u, v );
   auto err = std::abs( uv - uv_expected );
   WALBERLA_LOG_INFO_ON_ROOT( "|(u,v)_q - (u,v)_L2| = " << err );
   WALBERLA_CHECK_LESS( err, accuracy );
}

template < class MassOperator, class FE, typename ValueType, uint_t Quad >
void testRHS( const std::shared_ptr< PrimitiveStorage >&          storage,
              const std::function< ValueType( const Point3D& ) >& f,
              const uint_t&                                       lvl,
              const real_t&                                       min_error,
              const real_t&                                       max_error )
{
   L2Space< Quad, FE, ValueType > L2( storage, lvl );

   // operators and functions
   MassOperator M( storage, lvl, lvl );
   FE           f_proj( "f", storage, lvl, lvl );
   FE           b_M( "b_M", storage, lvl, lvl );
   FE           b_q( "b_q", storage, lvl, lvl );
   FE           e( "b_q - b_M", storage, lvl, lvl );

   // compare (φ_i, f)_i with [Mf]_i
   f_proj.interpolate( f, lvl, All );
   M.apply( f_proj, b_M, lvl, All );                  // compute rhs using mass matrix
   L2.dot( f, b_q );                                  // compute rhs using quadrature
   e.assign( { 1.0, -1.0 }, { b_q, b_M }, lvl, All ); // compute error
   auto err = e.getMaxMagnitude( lvl, All, true );    // max norm of error
   WALBERLA_LOG_INFO_ON_ROOT( "max_i |(v_i, f)_q - [Mf]_i| = " << err );
   WALBERLA_CHECK_LESS_EQUAL( min_error, err );
   WALBERLA_CHECK_LESS_EQUAL( err, max_error );
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
   std::function< real_t( const Point3D& ) >  u, v, f;
   std::function< Point3D( const Point3D& ) > u3, v3, f3;

   // -------
   //  TESTS
   // -------

   const real_t eps = real_c( std::is_same_v< real_t, double > ? 5e-12 : 5e-5 );

   using P1     = P1Function< real_t >;
   using P2     = P2Function< real_t >;
   using P1Mass = P1BlendingMassOperator;
   using P2Mass = P2ElementwiseBlendingMassOperator;

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

   u3 = [&]( const Point3D& x ) { return std::sin( std::atan2( x[1], x[0] ) ) * x; }; //u3(x) = sin(φ) x

   WALBERLA_LOG_INFO_ON_ROOT( "Test: Annulus" );
   testL2Dot< real_t, 5 >( storage, u, v, 60.0 * pi, 5, eps );
   testL2Dot< Point3D, 5 >( storage, u3, u3, 60.0 * pi, 5, eps );
   testRHS< P1Mass, P1, real_t, 5 >( storage, f, 5, 5e-6, 2e-5 );
   testRHS< P2Mass, P2, real_t, 7 >( storage, f, 5, 2e-9, real_c( std::is_same_v< real_t, double > ? 2e-8 : 6e-6 ) );

   // === test polynomials 3D ===

   meshInfo     = MeshInfo::meshCuboid( { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 }, 1, 1, 1 );
   setupStorage = SetupPrimitiveStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   storage      = std::make_shared< PrimitiveStorage >( setupStorage );

   Matrix3r C;
   C << 2.0, 3.0, 4.0, //
       3.0, 4.0, 5.0,  //
       4.0, 5.0, 6.0;
   Point3D c = C.diagonal();
   Point3D one{ 1.0, 1.0, 1.0 };
   u  = [&]( const Point3D& x ) { return one.dot( x ); };     // u(x) = x+y+z
   v  = [&]( const Point3D& x ) { return c.dot( x ); };       // v(x) = c⋅x
   u3 = [&]( const Point3D& x ) { return x; };                //u3(x) = x
   v3 = [&]( const Point3D& x ) -> Point3D { return C * x; }; //v3(x) = Cx

   WALBERLA_LOG_INFO_ON_ROOT( "Test: Polynomial 3D" );
   testL2Dot< real_t, 5 >( storage, u, v, 10.0, 4, eps );
   testL2Dot< Point3D, 5 >( storage, u3, v3, 10.0, 4, eps );
   testRHS< P1Mass, P1, real_t, 5 >( storage, v, 4, 0, eps );
   testRHS< P2Mass, P2, real_t, 7 >( storage, v, 4, 0, eps );

   return EXIT_SUCCESS;
}
