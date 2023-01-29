/*
 * Copyright (c) 2023 Daniel Bauer.
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

#include <cmath>
#include <sstream>

#include "core/DataTypes.h"
#include "core/math/Constants.h"
#include "core/math/Random.h"
#include "core/mpi/Environment.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/edgedofspace/EdgeDoFOrientation.hpp"
#include "hyteg/eigen/typeAliases.hpp"
#include "hyteg/n1e1functionspace/N1E1VectorFunction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

#include "Table.hpp"

using namespace hyteg;
using walberla::real_t;
using walberla::uint_t;

class Quadrature
{
   std::vector< std::tuple< real_t, real_t > > pointWeights; // wrt [-1, 1]
   Quadrature( std::vector< std::tuple< real_t, real_t > > weights )
   : pointWeights( weights ){};

 public:
   static Quadrature midPointRule() { return Quadrature( { { 0.0, 2.0 } } ); }
   static Quadrature gaussLegendre2()
   {
      return Quadrature( { { -std::sqrt( 1.0 / 3.0 ), 1.0 }, { std::sqrt( 1.0 / 3.0 ), 1.0 } } );
   }
   static Quadrature gaussLegendre3()
   {
      return Quadrature( { { 0.0, 8.0 / 9.0 }, { -std::sqrt( 3.0 / 5.0 ), 5.0 / 9.0 }, { std::sqrt( 3.0 / 5.0 ), 5.0 / 9.0 } } );
   }

   real_t integrate( const Eigen::Vector3r                                  start,
                     const Eigen::Vector3r                                  end,
                     const std::function< real_t( const Eigen::Vector3r ) > expr ) const
   {
      const real_t len = ( end - start ).norm();
      real_t       res = 0.0;
      for ( const auto& [point, weight] : pointWeights )
      {
         res += 0.5 * weight * len * expr( 0.5 * ( point + 1.0 ) * end - 0.5 * ( point - 1.0 ) * start );
      }
      return res;
   }
};

inline Eigen::Vector3r edgeTangent( const Cell& cell, const edgedof::EdgeDoFOrientation& orientation )
{
   const Eigen::Vector3r xDir = ( cell.getCoordinates()[1].vector_ - cell.getCoordinates()[0].vector_ );
   const Eigen::Vector3r yDir = ( cell.getCoordinates()[2].vector_ - cell.getCoordinates()[0].vector_ );
   const Eigen::Vector3r zDir = ( cell.getCoordinates()[3].vector_ - cell.getCoordinates()[0].vector_ );

   switch ( orientation )
   {
   case edgedof::EdgeDoFOrientation::X:
      return xDir.normalized();
   case edgedof::EdgeDoFOrientation::Y:
      return yDir.normalized();
   case edgedof::EdgeDoFOrientation::Z:
      return zDir.normalized();
   case edgedof::EdgeDoFOrientation::XY:
      return ( yDir - xDir ).normalized();
   case edgedof::EdgeDoFOrientation::XZ:
      return ( zDir - xDir ).normalized();
   case edgedof::EdgeDoFOrientation::YZ:
      return ( zDir - yDir ).normalized();
   case edgedof::EdgeDoFOrientation::XYZ:
      return ( zDir - yDir + xDir ).normalized();
   default:
      WALBERLA_ABORT( "wrong orienation" )
   }
}

void interpolate( const Quadrature&                                               quadrature,
                  const Cell&                                                     cell,
                  const PrimitiveDataID< FunctionMemory< real_t >, Cell >&        cellMemoryId,
                  const std::function< Eigen::Vector3r( const Eigen::Vector3r ) > expr,
                  const uint_t                                                    level )
{
   auto cellData = cell.getData( cellMemoryId )->getPointer( level );

   for ( const auto& it : edgedof::macrocell::Iterator( level ) )
   {
      const Point3D v0 = vertexdof::macrocell::coordinateFromIndex( level, cell, it );
      const Point3D v1 = vertexdof::macrocell::coordinateFromIndex( level, cell, it + indexing::IndexIncrement{ 1, 0, 0 } );
      const Point3D v2 = vertexdof::macrocell::coordinateFromIndex( level, cell, it + indexing::IndexIncrement{ 0, 1, 0 } );
      const Point3D v3 = vertexdof::macrocell::coordinateFromIndex( level, cell, it + indexing::IndexIncrement{ 0, 0, 1 } );

      // x ↦ ∫ₑ x·t dΓ
      const real_t dofScalarX  = quadrature.integrate( v0.vector_, v1.vector_, [&]( const Eigen::Vector3r x ) {
         return expr( x ).dot( edgeTangent( cell, edgedof::EdgeDoFOrientation::X ) );
      } );
      const real_t dofScalarY  = quadrature.integrate( v0.vector_, v2.vector_, [&]( const Eigen::Vector3r x ) {
         return expr( x ).dot( edgeTangent( cell, edgedof::EdgeDoFOrientation::Y ) );
      } );
      const real_t dofScalarZ  = quadrature.integrate( v0.vector_, v3.vector_, [&]( const Eigen::Vector3r x ) {
         return expr( x ).dot( edgeTangent( cell, edgedof::EdgeDoFOrientation::Z ) );
      } );
      const real_t dofScalarXY = quadrature.integrate( v1.vector_, v2.vector_, [&]( const Eigen::Vector3r x ) {
         return expr( x ).dot( edgeTangent( cell, edgedof::EdgeDoFOrientation::XY ) );
      } );
      const real_t dofScalarXZ = quadrature.integrate( v1.vector_, v3.vector_, [&]( const Eigen::Vector3r x ) {
         return expr( x ).dot( edgeTangent( cell, edgedof::EdgeDoFOrientation::XZ ) );
      } );
      const real_t dofScalarYZ = quadrature.integrate( v2.vector_, v3.vector_, [&]( const Eigen::Vector3r x ) {
         return expr( x ).dot( edgeTangent( cell, edgedof::EdgeDoFOrientation::YZ ) );
      } );

      cellData[edgedof::macrocell::xIndex( level, it.x(), it.y(), it.z() )]  = dofScalarX;
      cellData[edgedof::macrocell::yIndex( level, it.x(), it.y(), it.z() )]  = dofScalarY;
      cellData[edgedof::macrocell::zIndex( level, it.x(), it.y(), it.z() )]  = dofScalarZ;
      cellData[edgedof::macrocell::xyIndex( level, it.x(), it.y(), it.z() )] = dofScalarXY;
      cellData[edgedof::macrocell::xzIndex( level, it.x(), it.y(), it.z() )] = dofScalarXZ;
      cellData[edgedof::macrocell::yzIndex( level, it.x(), it.y(), it.z() )] = dofScalarYZ;
   }

   for ( const auto& it : edgedof::macrocell::IteratorXYZ( level ) )
   {
      const Point3D v0 = vertexdof::macrocell::coordinateFromIndex( level, cell, it + indexing::IndexIncrement{ 0, 1, 0 } );
      const Point3D v1 = vertexdof::macrocell::coordinateFromIndex( level, cell, it + indexing::IndexIncrement{ 1, 0, 1 } );

      const real_t dofScalarXYZ = quadrature.integrate( v0.vector_, v1.vector_, [&]( const Eigen::Vector3r x ) {
         return expr( x ).dot( edgeTangent( cell, edgedof::EdgeDoFOrientation::XYZ ) );
      } );

      cellData[edgedof::macrocell::xyzIndex( level, it.x(), it.y(), it.z() )] = dofScalarXYZ;
   }
}

void interpolate( const Quadrature&                                               quadrature,
                  const std::shared_ptr< PrimitiveStorage >                       storage,
                  n1e1::N1E1VectorFunction< real_t >&                             fun,
                  const std::function< Eigen::Vector3r( const Eigen::Vector3r ) > expr,
                  const uint_t                                                    level )
{
   for ( const auto& cellID : storage->getCellIDs() )
   {
      Cell& cell = *storage->getCell( cellID );
      interpolate( quadrature, cell, fun.getDoFs()->getCellDataID(), expr, level );
   }
}

real_t test( const Quadrature&                                               quadrature,
             const std::function< Eigen::Vector3r( const Eigen::Vector3r ) > expr,
             const uint_t                                                    level,
             const bool                                                      writeVTK = false )
{
   using namespace n1e1;

   MeshInfo              meshInfo = MeshInfo::meshSymmetricCuboid( Point3D( { 0, 0, 0 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   N1E1VectorFunction< real_t > u( "u", storage, level, level );
   interpolate( quadrature, storage, u, expr, level );

   if ( writeVTK )
   {
      VTKOutput vtk( "./output", "projectionQuadrature", storage );
      vtk.add( u );
      vtk.write( level );
   }

   walberla::math::seedRandomGenerator( 42 );
   const uint_t numRandomEvaluations = 1000;
   real_t       error                = 0.0;

   for ( uint_t i = 0; i < numRandomEvaluations; ++i )
   {
      Point3D coordinates;
      coordinates[0] = walberla::math::realRandom( 0.0, 1.0 );
      coordinates[1] = walberla::math::realRandom( 0.0, 1.0 );
      coordinates[2] = walberla::math::realRandom( 0.0, 1.0 );

      Eigen::Vector3r eval;
      auto            success = u.evaluate( coordinates, level, eval );
      WALBERLA_CHECK( success );

      error += ( eval - expr( coordinates.vector_ ) ).norm();
   }

   return error / numRandomEvaluations;
}

void projectionQuadrature()
{
   const uint_t minLevel = 3;
   const uint_t maxLevel = 7;

   auto expr = []( const Eigen::Vector3r p ) {
      using std::sin;
      using std::cos;
      using walberla::math::pi;

      const real_t kp = 4.0 * pi;

      const real_t x = p[0];
      const real_t y = p[1];
      const real_t z = p[2];

      return Eigen::Vector3r{ sin( x * kp ) * cos( y * kp ) * sin( z * kp ) - sin( x * kp ) * sin( y * kp ) * cos( z * kp ),
                              sin( x * kp ) * sin( y * kp ) * cos( z * kp ) - cos( x * kp ) * sin( y * kp ) * sin( z * kp ),
                              cos( x * kp ) * sin( y * kp ) * sin( z * kp ) - sin( x * kp ) * cos( y * kp ) * sin( z * kp ) };
   };

   Table< 4 > table( { "level", "gl1", "gl2", "gl3" } );

   for ( uint_t level = minLevel; level <= maxLevel; ++level )
   {
      const real_t gl1 = test( Quadrature::midPointRule(), expr, level );
      const real_t gl2 = test( Quadrature::gaussLegendre2(), expr, level );
      const real_t gl3 = test( Quadrature::gaussLegendre3(), expr, level );

      table.addElement( level - minLevel, 0, level );
      table.addElement( level - minLevel, 1, gl1 );
      table.addElement( level - minLevel, 2, gl2 );
      table.addElement( level - minLevel, 3, gl3 );
   }

   WALBERLA_LOG_INFO_ON_ROOT( table )
   table.write( "output", "projectionQuadrature" );
}

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::mpi::MPIManager::instance()->useWorldComm();

   projectionQuadrature();
}
