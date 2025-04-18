/*
 * Copyright (c) 2024 Andreas Burkhart, Marcus Mohr, Ponsuganth Ilangovan P.
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

#include <random>

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/geometry/AffineMap2D.hpp"
#include "hyteg/geometry/AffineMap3D.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/geometry/IdentityMap.hpp"
#include "hyteg/geometry/PolarCoordsMap.hpp"
#include "hyteg/geometry/SphericalCoordsMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"

using walberla::real_t;
using walberla::math::pi;

using namespace hyteg;

auto createAnnulusStorage( uint_t nTan, uint_t nRad )
{
   const real_t rMin = real_c( 0.5 );
   const real_t rMax = real_c( 1.5 );

   MeshInfo meshInfo = MeshInfo::meshAnnulus( rMin, rMax, MeshInfo::CRISS, nTan, nRad );

   auto setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   AnnulusMap::setMap( *setupStorage );

   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage, 3 );

   return storage;
}

auto createPolarCoordsStorage( uint_t nSubCells )
{
   MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( real_c( 0.5 ), real_c( 0 ) ),
                                                Point2D( real_c( 1 ), real_c( 7.0 / 4.0 * pi ) ),
                                                MeshInfo::CRISS,
                                                nSubCells,
                                                nSubCells );

   auto setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   PolarCoordsMap::setMap( *setupStorage );

   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage, 3 );

   return storage;
}

auto createThickSphericalShellStorage( uint_t nTan, uint_t nRad )
{
   const real_t rMin = real_c( 0.5 ), rMax = real_c( 1.5 );
   MeshInfo     meshInfo = MeshInfo::meshSphericalShell( nTan, nRad, rMin, rMax );

   auto setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   IcosahedralShellMap::setMap( *setupStorage );
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage, 3 );

   return storage;
}

auto createSphericalCoordsStorage( uint_t nSubCells )
{
   MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid( Point3D( real_c( 1 ), real_c( 0 ), real_c( 0 ) ),
                                                      Point3D( real_c( 2 ), real_c( pi ), real_c( 2 * pi ) ),
                                                      nSubCells,
                                                      nSubCells,
                                                      nSubCells );

   auto setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   SphericalCoordsMap::setMap( *setupStorage );
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage, 3 );

   return storage;
}

real_t firstDerivativeTestTaylor2D( const std::shared_ptr< PrimitiveStorage >& storage )
{
   std::shared_ptr< Face > face;
   real_t                  maximum = std::numeric_limits< real_t >::min();
   real_t                  h       = real_c( 1e-2 );

   for ( auto& it : storage->getFaces() )
   {
      face = it.second;

      std::array< hyteg::Point3D, 3UL > faceCoordinates = face->getCoordinates();

      Point3D a = ( faceCoordinates[0] + faceCoordinates[1] + faceCoordinates[2] ) / real_c( 3.0 );

      Point3D dx( walberla::math::realRandom(), walberla::math::realRandom(), real_c( 0.0 ) );
      dx.normalize();
      dx = h * dx;

      Point3D x = a + dx;

      auto geometryMap = face->getGeometryMap();

      Point3D  Fa;
      Matrix2r Dfa;

      Point3D Tx, TxTaylor;

      geometryMap->evalF( a, Fa );
      geometryMap->evalF( x, Tx );

      geometryMap->evalDF( a, Dfa );

      TxTaylor( 0 ) = Fa( 0 ) + ( dx( 0 ) * Dfa( 0, 0 ) + dx( 1 ) * Dfa( 1, 0 ) );
      TxTaylor( 1 ) = Fa( 1 ) + ( dx( 0 ) * Dfa( 0, 1 ) + dx( 1 ) * Dfa( 1, 1 ) );

      const real_t error = ( Tx - TxTaylor ).norm();

      if ( error > maximum )
      {
         maximum = error;
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT( "firstDerivativeTestTaylor2D: " << maximum );

   return maximum;
}

real_t firstDerivativeTestTaylor3D( const std::shared_ptr< PrimitiveStorage >& storage )
{
   std::shared_ptr< Cell > cell;
   real_t                  maximum = std::numeric_limits< real_t >::min();
   real_t                  h       = real_c( 1e-2 );

   for ( auto& it : storage->getCells() )
   {
      cell = it.second;

      std::array< hyteg::Point3D, 4UL > cellCoordinates = cell->getCoordinates();

      Point3D a = ( cellCoordinates[0] + cellCoordinates[1] + cellCoordinates[2] + cellCoordinates[3] ) / real_c( 4.0 );

      Point3D dx( walberla::math::realRandom(), walberla::math::realRandom(), walberla::math::realRandom() );
      dx.normalize();
      dx = h * dx;

      Point3D x = a + dx;

      auto geometryMap = cell->getGeometryMap();

      Point3D  Fa;
      Matrix3r Dfa;

      Point3D Tx, TxTaylor;

      geometryMap->evalF( a, Fa );
      geometryMap->evalF( x, Tx );

      geometryMap->evalDF( a, Dfa );

      TxTaylor( 0 ) = Fa( 0 ) + ( dx( 0 ) * Dfa( 0, 0 ) + dx( 1 ) * Dfa( 1, 0 ) + dx( 2 ) * Dfa( 2, 0 ) );
      TxTaylor( 1 ) = Fa( 1 ) + ( dx( 0 ) * Dfa( 0, 1 ) + dx( 1 ) * Dfa( 1, 1 ) + dx( 2 ) * Dfa( 2, 1 ) );
      TxTaylor( 2 ) = Fa( 2 ) + ( dx( 0 ) * Dfa( 0, 2 ) + dx( 1 ) * Dfa( 1, 2 ) + dx( 2 ) * Dfa( 2, 2 ) );

      const real_t error = ( Tx - TxTaylor ).norm();

      if ( error > maximum )
      {
         maximum = error;
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT( "firstDerivativeTestTaylor3D: " << maximum );

   return maximum;
}

real_t secondDerivativesFiniteDifferenceTest2D( const std::shared_ptr< PrimitiveStorage >& storage )
{
   std::shared_ptr< Face > face;
   real_t                  maximum = std::numeric_limits< real_t >::min();

   for ( auto& it : storage->getFaces() )
   {
      face = it.second;

      std::array< hyteg::Point3D, 3UL > faceCoordinates = face->getCoordinates();

      Point3D x = ( faceCoordinates[0] + faceCoordinates[1] + faceCoordinates[2] ) / real_c( 3.0 );
      Point3D evalX;

      Matrix2r DF, DFinv;

      std::shared_ptr< GeometryMap > geometryMap = face->getGeometryMap();

      Matrixr< 2, 4 > DFInvDFxNumerical, DFInvDFx;

      //WALBERLA_LOG_INFO_ON_ROOT( "EvalDFinvDf:" );
      geometryMap->evalDFinvDF( x, DFInvDFx );

      Point3D dx( real_c( 0.0 ), real_c( 0.0 ), real_c( 0.0 ) );
      real_t  h = real_c( 1e-2 );

      Matrix2r DInvFx_plus_dx, DInvFx_minus_dx;

      dx( 0 ) = h;
      geometryMap->evalDFinv( x + dx, DInvFx_plus_dx );
      geometryMap->evalDFinv( x - dx, DInvFx_minus_dx );

      DFInvDFxNumerical( 0, 0 ) = ( DInvFx_plus_dx( 0, 0 ) - DInvFx_minus_dx( 0, 0 ) ) / ( real_c( 2.0 ) * dx( 0 ) );
      DFInvDFxNumerical( 1, 0 ) = ( DInvFx_plus_dx( 1, 0 ) - DInvFx_minus_dx( 1, 0 ) ) / ( real_c( 2.0 ) * dx( 0 ) );
      DFInvDFxNumerical( 0, 1 ) = ( DInvFx_plus_dx( 0, 1 ) - DInvFx_minus_dx( 0, 1 ) ) / ( real_c( 2.0 ) * dx( 0 ) );
      DFInvDFxNumerical( 1, 1 ) = ( DInvFx_plus_dx( 1, 1 ) - DInvFx_minus_dx( 1, 1 ) ) / ( real_c( 2.0 ) * dx( 0 ) );

      dx( 0 ) = real_c( 0.0 );
      dx( 1 ) = h;

      geometryMap->evalDFinv( x + dx, DInvFx_plus_dx );
      geometryMap->evalDFinv( x - dx, DInvFx_minus_dx );

      DFInvDFxNumerical( 0, 2 ) = ( DInvFx_plus_dx( 0, 0 ) - DInvFx_minus_dx( 0, 0 ) ) / ( real_c( 2.0 ) * dx( 1 ) );
      DFInvDFxNumerical( 1, 2 ) = ( DInvFx_plus_dx( 1, 0 ) - DInvFx_minus_dx( 1, 0 ) ) / ( real_c( 2.0 ) * dx( 1 ) );
      DFInvDFxNumerical( 0, 3 ) = ( DInvFx_plus_dx( 0, 1 ) - DInvFx_minus_dx( 0, 1 ) ) / ( real_c( 2.0 ) * dx( 1 ) );
      DFInvDFxNumerical( 1, 3 ) = ( DInvFx_plus_dx( 1, 1 ) - DInvFx_minus_dx( 1, 1 ) ) / ( real_c( 2.0 ) * dx( 1 ) );

      Matrixr< 2, 4 > Difference = DFInvDFx - DFInvDFxNumerical;

      real_t frobenius = std::sqrt( Difference.squaredNorm() );

      if ( frobenius > maximum )
      {
         maximum = frobenius;
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT( "secondDerivativesFiniteDifferenceTest2D: " << maximum );

   return maximum;
}

real_t secondDerivativesFiniteDifferenceTest3D( const std::shared_ptr< PrimitiveStorage >& storage )
{
   std::shared_ptr< Cell > cell;
   real_t                  maximum = std::numeric_limits< real_t >::min();
   for ( auto& it : storage->getCells() )
   {
      cell = it.second;

      std::array< hyteg::Point3D, 4UL > cellCoordinates = cell->getCoordinates();

      Point3D x = ( cellCoordinates[0] + cellCoordinates[1] + cellCoordinates[2] + cellCoordinates[3] ) / real_c( 4.0 );

      Point3D evalX;

      Matrix3r DF, DFinv;

      std::shared_ptr< GeometryMap > geometryMap = cell->getGeometryMap();

      Matrixr< 3, 9 > DFInvDFxNumerical, DFInvDFx;
      Matrix3r        DInvFx_plus_dx, DInvFx_minus_dx;

      Point3D dx( real_c( 0.0 ), real_c( 0.0 ), real_c( 0.0 ) );
      real_t  h = real_c( 1e-2 );

      dx( 0 ) = h;

      geometryMap->evalDF( x + dx, DInvFx_plus_dx );
      DInvFx_plus_dx = DInvFx_plus_dx.inverse().eval();

      geometryMap->evalDF( x - dx, DInvFx_minus_dx );
      DInvFx_minus_dx = DInvFx_minus_dx.inverse().eval();

      DFInvDFxNumerical( 0, 0 ) = ( DInvFx_plus_dx( 0, 0 ) - DInvFx_minus_dx( 0, 0 ) ) / ( real_c( 2.0 ) * dx( 0 ) );
      DFInvDFxNumerical( 1, 0 ) = ( DInvFx_plus_dx( 1, 0 ) - DInvFx_minus_dx( 1, 0 ) ) / ( real_c( 2.0 ) * dx( 0 ) );
      DFInvDFxNumerical( 2, 0 ) = ( DInvFx_plus_dx( 2, 0 ) - DInvFx_minus_dx( 2, 0 ) ) / ( real_c( 2.0 ) * dx( 0 ) );
      DFInvDFxNumerical( 0, 1 ) = ( DInvFx_plus_dx( 0, 1 ) - DInvFx_minus_dx( 0, 1 ) ) / ( real_c( 2.0 ) * dx( 0 ) );
      DFInvDFxNumerical( 1, 1 ) = ( DInvFx_plus_dx( 1, 1 ) - DInvFx_minus_dx( 1, 1 ) ) / ( real_c( 2.0 ) * dx( 0 ) );
      DFInvDFxNumerical( 2, 1 ) = ( DInvFx_plus_dx( 2, 1 ) - DInvFx_minus_dx( 2, 1 ) ) / ( real_c( 2.0 ) * dx( 0 ) );
      DFInvDFxNumerical( 0, 2 ) = ( DInvFx_plus_dx( 0, 2 ) - DInvFx_minus_dx( 0, 2 ) ) / ( real_c( 2.0 ) * dx( 0 ) );
      DFInvDFxNumerical( 1, 2 ) = ( DInvFx_plus_dx( 1, 2 ) - DInvFx_minus_dx( 1, 2 ) ) / ( real_c( 2.0 ) * dx( 0 ) );
      DFInvDFxNumerical( 2, 2 ) = ( DInvFx_plus_dx( 2, 2 ) - DInvFx_minus_dx( 2, 2 ) ) / ( real_c( 2.0 ) * dx( 0 ) );

      dx( 0 ) = real_c( 0.0 );
      dx( 1 ) = h;

      geometryMap->evalDF( x + dx, DInvFx_plus_dx );
      DInvFx_plus_dx = DInvFx_plus_dx.inverse().eval();

      geometryMap->evalDF( x - dx, DInvFx_minus_dx );
      DInvFx_minus_dx = DInvFx_minus_dx.inverse().eval();

      DFInvDFxNumerical( 0, 3 ) = ( DInvFx_plus_dx( 0, 0 ) - DInvFx_minus_dx( 0, 0 ) ) / ( real_c( 2.0 ) * dx( 1 ) );
      DFInvDFxNumerical( 1, 3 ) = ( DInvFx_plus_dx( 1, 0 ) - DInvFx_minus_dx( 1, 0 ) ) / ( real_c( 2.0 ) * dx( 1 ) );
      DFInvDFxNumerical( 2, 3 ) = ( DInvFx_plus_dx( 2, 0 ) - DInvFx_minus_dx( 2, 0 ) ) / ( real_c( 2.0 ) * dx( 1 ) );
      DFInvDFxNumerical( 0, 4 ) = ( DInvFx_plus_dx( 0, 1 ) - DInvFx_minus_dx( 0, 1 ) ) / ( real_c( 2.0 ) * dx( 1 ) );
      DFInvDFxNumerical( 1, 4 ) = ( DInvFx_plus_dx( 1, 1 ) - DInvFx_minus_dx( 1, 1 ) ) / ( real_c( 2.0 ) * dx( 1 ) );
      DFInvDFxNumerical( 2, 4 ) = ( DInvFx_plus_dx( 2, 1 ) - DInvFx_minus_dx( 2, 1 ) ) / ( real_c( 2.0 ) * dx( 1 ) );
      DFInvDFxNumerical( 0, 5 ) = ( DInvFx_plus_dx( 0, 2 ) - DInvFx_minus_dx( 0, 2 ) ) / ( real_c( 2.0 ) * dx( 1 ) );
      DFInvDFxNumerical( 1, 5 ) = ( DInvFx_plus_dx( 1, 2 ) - DInvFx_minus_dx( 1, 2 ) ) / ( real_c( 2.0 ) * dx( 1 ) );
      DFInvDFxNumerical( 2, 5 ) = ( DInvFx_plus_dx( 2, 2 ) - DInvFx_minus_dx( 2, 2 ) ) / ( real_c( 2.0 ) * dx( 1 ) );

      dx( 1 ) = real_c( 0.0 );
      dx( 2 ) = h;

      geometryMap->evalDF( x + dx, DInvFx_plus_dx );
      DInvFx_plus_dx = DInvFx_plus_dx.inverse().eval();

      geometryMap->evalDF( x - dx, DInvFx_minus_dx );
      DInvFx_minus_dx = DInvFx_minus_dx.inverse().eval();

      DFInvDFxNumerical( 0, 6 ) = ( DInvFx_plus_dx( 0, 0 ) - DInvFx_minus_dx( 0, 0 ) ) / ( real_c( 2.0 ) * dx( 2 ) );
      DFInvDFxNumerical( 1, 6 ) = ( DInvFx_plus_dx( 1, 0 ) - DInvFx_minus_dx( 1, 0 ) ) / ( real_c( 2.0 ) * dx( 2 ) );
      DFInvDFxNumerical( 2, 6 ) = ( DInvFx_plus_dx( 2, 0 ) - DInvFx_minus_dx( 2, 0 ) ) / ( real_c( 2.0 ) * dx( 2 ) );
      DFInvDFxNumerical( 0, 7 ) = ( DInvFx_plus_dx( 0, 1 ) - DInvFx_minus_dx( 0, 1 ) ) / ( real_c( 2.0 ) * dx( 2 ) );
      DFInvDFxNumerical( 1, 7 ) = ( DInvFx_plus_dx( 1, 1 ) - DInvFx_minus_dx( 1, 1 ) ) / ( real_c( 2.0 ) * dx( 2 ) );
      DFInvDFxNumerical( 2, 7 ) = ( DInvFx_plus_dx( 2, 1 ) - DInvFx_minus_dx( 2, 1 ) ) / ( real_c( 2.0 ) * dx( 2 ) );
      DFInvDFxNumerical( 0, 8 ) = ( DInvFx_plus_dx( 0, 2 ) - DInvFx_minus_dx( 0, 2 ) ) / ( real_c( 2.0 ) * dx( 2 ) );
      DFInvDFxNumerical( 1, 8 ) = ( DInvFx_plus_dx( 1, 2 ) - DInvFx_minus_dx( 1, 2 ) ) / ( real_c( 2.0 ) * dx( 2 ) );
      DFInvDFxNumerical( 2, 8 ) = ( DInvFx_plus_dx( 2, 2 ) - DInvFx_minus_dx( 2, 2 ) ) / ( real_c( 2.0 ) * dx( 2 ) );

      geometryMap->evalDFinvDF( x, DFInvDFx );

      Matrixr< 3, 9 > Difference = DFInvDFx - DFInvDFxNumerical;

      real_t frobenius = std::sqrt( Difference.squaredNorm() );

      if ( frobenius > maximum )
      {
         maximum = frobenius;
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT( "secondDerivativesFiniteDifferenceTest3D: " << maximum );

   return maximum;
}

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   walberla::math::seedRandomGenerator( 42 );

   real_t maximum;

   WALBERLA_LOG_INFO_ON_ROOT( "----------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "Running tests for AnnulusMap" );
   WALBERLA_LOG_INFO_ON_ROOT( "----------------------------" );
   auto annulusStorage = createAnnulusStorage( 12u, 4u );

   maximum = firstDerivativeTestTaylor2D( annulusStorage );
   WALBERLA_CHECK_LESS( maximum, real_c( 1e-2 ) );

   maximum = secondDerivativesFiniteDifferenceTest2D( annulusStorage );
   WALBERLA_CHECK_LESS( maximum, real_c( 1e-3 ) );

   WALBERLA_LOG_INFO_ON_ROOT( "-------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "Running tests for IcosahedralShellMap" );
   WALBERLA_LOG_INFO_ON_ROOT( "-------------------------------------" );
   auto thickShellStorage = createThickSphericalShellStorage( 3u, 2u );

   maximum = firstDerivativeTestTaylor3D( thickShellStorage );
   WALBERLA_CHECK_LESS( maximum, real_c( 1e-2 ) );

   maximum = secondDerivativesFiniteDifferenceTest3D( thickShellStorage );
   WALBERLA_CHECK_LESS( maximum, real_c( 1e-3 ) );

   WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "Running tests for PolarCoordsMap" );
   WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------" );
   auto polarCoordsStorage = createPolarCoordsStorage( 3u );

   maximum = firstDerivativeTestTaylor2D( polarCoordsStorage );
   WALBERLA_CHECK_LESS( maximum, real_c( 2e-2 ) );

   maximum = secondDerivativesFiniteDifferenceTest2D( polarCoordsStorage );
   WALBERLA_CHECK_LESS( maximum, real_c( 1.2e-3 ) );

   WALBERLA_LOG_INFO_ON_ROOT( "------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "Running tests for SphericalCoordsMap" );
   WALBERLA_LOG_INFO_ON_ROOT( "------------------------------------" );
   auto sphericalCoordsStorage = createSphericalCoordsStorage( 2u );

   maximum = firstDerivativeTestTaylor3D( sphericalCoordsStorage );
   WALBERLA_CHECK_LESS( maximum, real_c( 5e-2 ) );

   maximum = secondDerivativesFiniteDifferenceTest3D( sphericalCoordsStorage );
   WALBERLA_CHECK_LESS( maximum, real_c( 7e-2 ) );

   WALBERLA_LOG_INFO_ON_ROOT( "----------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "Running test for AffineMap2D" );
   WALBERLA_LOG_INFO_ON_ROOT( "----------------------------" );
   {
      Matrix2r mat;
      real_t   phi = 2.0 / 9.0 * pi;
      mat( 0, 0 )  = +std::cos( phi );
      mat( 0, 1 )  = -std::sin( phi );
      mat( 1, 0 )  = +std::sin( phi ) * 2.25;
      mat( 1, 1 )  = +std::cos( phi ) * 2.25;
      Point2D         vec( -7.0, 3.0 );
      AffineMap2D     map( mat, vec );
      Point3D         point( { real_c( 0 ), real_c( 0 ), real_c( 0 ) } );
      Matrixr< 2, 4 > DFinvDFx = Matrixr< 2, 4 >::Random();
      ;
      map.evalDFinvDF( point, DFinvDFx );
      WALBERLA_CHECK_FLOAT_EQUAL( DFinvDFx.norm(), real_c( 0 ) );
   }
   WALBERLA_LOG_INFO_ON_ROOT( "passed" );

   WALBERLA_LOG_INFO_ON_ROOT( "----------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "Running test for AffineMap3D" );
   WALBERLA_LOG_INFO_ON_ROOT( "----------------------------" );
   {
      Matrix3r mat;
      mat( 0, 0 ) = +8.660254037844387e-01;
      mat( 0, 1 ) = -1.545084971874737e-01;
      mat( 0, 2 ) = -9.510565162951534e-01;
      mat( 1, 0 ) = +0.000000000000000e+00;
      mat( 1, 1 ) = +9.510565162951535e-01;
      mat( 1, 2 ) = -6.180339887498948e-01;
      mat( 2, 0 ) = +4.999999999999999e-01;
      mat( 2, 1 ) = +2.676165673298174e-01;
      mat( 2, 2 ) = +1.647278207092664e+00;
      Point3D         vec( -7.0, 3.0, 2.0 );
      AffineMap3D     map( mat, vec );
      Point3D         point( { real_c( 0 ), real_c( 0 ), real_c( 0 ) } );
      Matrixr< 3, 9 > DFinvDFx = Matrixr< 3, 9 >::Random();
      map.evalDFinvDF( point, DFinvDFx );
      WALBERLA_CHECK_FLOAT_EQUAL( DFinvDFx.norm(), real_c( 0 ) );
   }
   WALBERLA_LOG_INFO_ON_ROOT( "passed" );

   WALBERLA_LOG_INFO_ON_ROOT( "----------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "Running test for IdentityMap" );
   WALBERLA_LOG_INFO_ON_ROOT( "----------------------------" );
   {
      IdentityMap     map;
      Point3D         point( { real_c( 0 ), real_c( 0 ), real_c( 0 ) } );
      Matrixr< 2, 4 > derivs2D = Matrixr< 2, 4 >::Random();
      Matrixr< 3, 9 > derivs3D = Matrixr< 3, 9 >::Random();

      map.evalDFinvDF( point, derivs2D );
      WALBERLA_CHECK_FLOAT_EQUAL( derivs2D.norm(), real_c( 0 ) );
      WALBERLA_LOG_INFO_ON_ROOT( "2D passed" );

      map.evalDFinvDF( point, derivs3D );
      WALBERLA_CHECK_FLOAT_EQUAL( derivs3D.norm(), real_c( 0 ) );
      WALBERLA_LOG_INFO_ON_ROOT( "3D passed" );
   }
   return 0;
}
