/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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
#include "VertexDoFFunction.hpp"

#include <typeinfo>
#include <utility>

#include "core/OpenMP.h"

#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/facedofspace/FaceDoFFunction.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/geometry/Intersection.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/p1functionspace/VertexDoFAdditivePackInfo.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/p1functionspace/VertexDoFPackInfo.hpp"
#include "hyteg/p1functionspace/generatedKernels/add_2D_macroface_vertexdof_1_rhsfunction.hpp"
#include "hyteg/p1functionspace/generatedKernels/add_2D_macroface_vertexdof_2_rhsfunctions.hpp"
#include "hyteg/p1functionspace/generatedKernels/add_2D_macroface_vertexdof_3_rhsfunctions.hpp"
#include "hyteg/p1functionspace/generatedKernels/add_3D_macrocell_vertexdof_1_rhsfunction.hpp"
#include "hyteg/p1functionspace/generatedKernels/assign_2D_macroface_vertexdof_1_rhsfunction.hpp"
#include "hyteg/p1functionspace/generatedKernels/assign_2D_macroface_vertexdof_2_rhsfunctions.hpp"
#include "hyteg/p1functionspace/generatedKernels/assign_2D_macroface_vertexdof_3_rhsfunctions.hpp"
#include "hyteg/p1functionspace/generatedKernels/assign_3D_macrocell_vertexdof_1_rhsfunction.hpp"
#include "hyteg/p1functionspace/generatedKernels/assign_3D_macrocell_vertexdof_2_rhsfunctions.hpp"
#include "hyteg/p1functionspace/generatedKernels/assign_3D_macrocell_vertexdof_3_rhsfunctions.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"

namespace hyteg {
namespace vertexdof {

using walberla::int_c;

template < typename ValueType >
VertexDoFFunction< ValueType >::VertexDoFFunction( const std::string& name, const std::shared_ptr< PrimitiveStorage >& storage )
: Function< VertexDoFFunction< ValueType > >( name, storage )
, vertexDataID_( storage->generateInvalidPrimitiveDataID< MemoryDataHandling< FunctionMemory< ValueType >, Vertex >, Vertex >() )
, edgeDataID_( storage->generateInvalidPrimitiveDataID< MemoryDataHandling< FunctionMemory< ValueType >, Edge >, Edge >() )
, faceDataID_( storage->generateInvalidPrimitiveDataID< MemoryDataHandling< FunctionMemory< ValueType >, Face >, Face >() )
, cellDataID_( storage->generateInvalidPrimitiveDataID< MemoryDataHandling< FunctionMemory< ValueType >, Cell >, Cell >() )
{}

template < typename ValueType >
VertexDoFFunction< ValueType >::VertexDoFFunction( const std::string&                         name,
                                                   const std::shared_ptr< PrimitiveStorage >& storage,
                                                   uint_t                                     minLevel,
                                                   uint_t                                     maxLevel )
: VertexDoFFunction( name, storage, minLevel, maxLevel, BoundaryCondition::create0123BC() )
{}

template < typename ValueType >
VertexDoFFunction< ValueType >::VertexDoFFunction( const std::string&                         name,
                                                   const std::shared_ptr< PrimitiveStorage >& storage,
                                                   uint_t                                     minLevel,
                                                   uint_t                                     maxLevel,
                                                   BoundaryCondition                          boundaryCondition )
: Function< VertexDoFFunction< ValueType > >( name, storage, minLevel, maxLevel )
, boundaryCondition_( std::move( boundaryCondition ) )
{
   auto cellVertexDoFFunctionMemoryDataHandling = std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Cell > >();
   auto faceVertexDoFFunctionMemoryDataHandling = std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Face > >();
   auto edgeVertexDoFFunctionMemoryDataHandling = std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Edge > >();
   auto vertexVertexDoFFunctionMemoryDataHandling = std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Vertex > >();

   storage->addCellData( cellDataID_, cellVertexDoFFunctionMemoryDataHandling, name );
   storage->addFaceData( faceDataID_, faceVertexDoFFunctionMemoryDataHandling, name );
   storage->addEdgeData( edgeDataID_, edgeVertexDoFFunctionMemoryDataHandling, name );
   storage->addVertexData( vertexDataID_, vertexVertexDoFFunctionMemoryDataHandling, name );

   for ( uint_t level = minLevel; level <= maxLevel; ++level )
   {
      for ( const auto & it : storage->getVertices() )
      {
         allocateMemory( level, *it.second );
      }
      for ( const auto & it : storage->getEdges() )
      {
         allocateMemory( level, *it.second );
      }
      for ( const auto & it : storage->getFaces() )
      {
         allocateMemory( level, *it.second );
      }
      for ( const auto & it : storage->getCells() )
      {
         allocateMemory( level, *it.second );
      }

      communicators_[level]->addPackInfo( std::make_shared< VertexDoFPackInfo< ValueType > >(
          level, vertexDataID_, edgeDataID_, faceDataID_, cellDataID_, this->getStorage() ) );
      additiveCommunicators_[level]->addPackInfo( std::make_shared< VertexDoFAdditivePackInfo< ValueType > >(
          level, vertexDataID_, edgeDataID_, faceDataID_, cellDataID_, this->getStorage() ) );
   }
}

template < typename ValueType >
bool VertexDoFFunction< ValueType >::hasMemoryAllocated( const uint_t & level, const Vertex & vertex ) const
{
   WALBERLA_CHECK( this->getStorage()->vertexExistsLocally( vertex.getID() ) );
   return vertex.hasData( getVertexDataID() ) && vertex.getData( getVertexDataID() )->hasLevel( level );
}

template < typename ValueType >
bool VertexDoFFunction< ValueType >::hasMemoryAllocated( const uint_t & level, const Edge & edge ) const
{
   WALBERLA_CHECK( this->getStorage()->edgeExistsLocally( edge.getID() ) );
   return edge.hasData( getEdgeDataID() ) && edge.getData( getEdgeDataID() )->hasLevel( level );
}

template < typename ValueType >
bool VertexDoFFunction< ValueType >::hasMemoryAllocated( const uint_t & level, const Face & face ) const
{
   WALBERLA_CHECK( this->getStorage()->faceExistsLocally( face.getID() ) );
   return face.hasData( getFaceDataID() ) && face.getData( getFaceDataID() )->hasLevel( level );
}

template < typename ValueType >
bool VertexDoFFunction< ValueType >::hasMemoryAllocated( const uint_t & level, const Cell & cell ) const
{
   WALBERLA_CHECK( this->getStorage()->cellExistsLocally( cell.getID() ) );
   return cell.hasData( getCellDataID() ) && cell.getData( getCellDataID() )->hasLevel( level );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::allocateMemory( const uint_t & level, const Vertex & vertex )
{
   WALBERLA_CHECK( this->getStorage()->vertexExistsLocally( vertex.getID() ) );
   WALBERLA_CHECK( vertex.hasData( getVertexDataID() ) )
   if ( hasMemoryAllocated( level, vertex ) )
      return;
   vertex.getData( getVertexDataID() )->addData( level, vertexDoFMacroVertexFunctionMemorySize( level, vertex ), 0 );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::allocateMemory( const uint_t & level, const Edge & edge )
{
   WALBERLA_CHECK( this->getStorage()->edgeExistsLocally( edge.getID() ) );
   WALBERLA_CHECK( edge.hasData( getEdgeDataID() ) )
   if ( hasMemoryAllocated( level, edge ) )
      return;
   edge.getData( getEdgeDataID() )->addData( level, vertexDoFMacroEdgeFunctionMemorySize( level, edge ), 0 );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::allocateMemory( const uint_t & level, const Face & face )
{
   WALBERLA_CHECK( this->getStorage()->faceExistsLocally( face.getID() ) );
   WALBERLA_CHECK( face.hasData( getFaceDataID() ) )
   if ( hasMemoryAllocated( level, face ) )
      return;
   face.getData( getFaceDataID() )->addData( level, vertexDoFMacroFaceFunctionMemorySize( level, face ), 0 );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::allocateMemory( const uint_t & level, const Cell & cell )
{
   WALBERLA_CHECK( this->getStorage()->cellExistsLocally( cell.getID() ) );
   WALBERLA_CHECK( cell.hasData( getCellDataID() ) )
   if ( hasMemoryAllocated( level, cell ) )
      return;
   cell.getData( getCellDataID() )->addData( level, vertexDoFMacroCellFunctionMemorySize( level, cell ), 0 );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::deleteMemory( const uint_t & level, const Vertex & vertex )
{
   WALBERLA_CHECK( this->getStorage()->vertexExistsLocally( vertex.getID() ) );
   if ( !hasMemoryAllocated( level, vertex ) )
      return;
   vertex.getData( getVertexDataID() )->deleteData( level );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::deleteMemory( const uint_t & level, const Edge & edge )
{
   WALBERLA_CHECK( this->getStorage()->edgeExistsLocally( edge.getID() ) );
   if ( !hasMemoryAllocated( level, edge ) )
      return;
   edge.getData( getEdgeDataID() )->deleteData( level );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::deleteMemory( const uint_t & level, const Face & face )
{
   WALBERLA_CHECK( this->getStorage()->faceExistsLocally( face.getID() ) );
   if ( !hasMemoryAllocated( level, face ) )
      return;
   face.getData( getFaceDataID() )->deleteData( level );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::deleteMemory( const uint_t & level, const Cell & cell )
{
   WALBERLA_CHECK( this->getStorage()->cellExistsLocally( cell.getID() ) );
   if ( !hasMemoryAllocated( level, cell ) )
      return;
   cell.getData( getCellDataID() )->deleteData( level );
}


template < typename ValueType >
BoundaryCondition VertexDoFFunction< ValueType >::getBoundaryCondition() const
{
   return boundaryCondition_;
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::setBoundaryCondition( BoundaryCondition bc )
{
   boundaryCondition_ = bc;
}


template < typename ValueType >
void VertexDoFFunction< ValueType >::interpolate( ValueType constant, uint_t level, DoFType flag ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "Interpolate" );

   interpolateByPrimitiveType< Vertex >( constant, level, flag );
   interpolateByPrimitiveType< Edge >( constant, level, flag );
   interpolateByPrimitiveType< Face >( constant, level, flag );
   interpolateByPrimitiveType< Cell >( constant, level, flag );

   this->stopTiming( "Interpolate" );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::interpolate( ValueType constant, uint_t level, BoundaryUID boundaryUID ) const
{
  interpolate( [constant]( const Point3D& ){ return constant; }, level, boundaryUID );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::interpolate( const std::function< ValueType( const Point3D& ) >& expr,
                                                  uint_t                                              level,
                                                  DoFType                                             flag ) const
{
   if ( isDummy() )
   {
      return;
   }
   std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) > exprExtended =
       [&expr]( const hyteg::Point3D& x, const std::vector< ValueType >& ) { return expr( x ); };
   interpolate( exprExtended, {}, level, flag );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::interpolate( const std::function< ValueType( const Point3D& ) >& expr,
                                                  uint_t                                              level,
                                                  BoundaryUID                                         boundaryUID ) const
{
   if ( isDummy() )
   {
      return;
   }
   std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) > exprExtended =
       [&expr]( const hyteg::Point3D& x, const std::vector< ValueType >& ) { return expr( x ); };
   interpolate( exprExtended, {}, level, boundaryUID );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::interpolate(
    const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
    const std::vector< std::reference_wrapper< const VertexDoFFunction< ValueType > > >& srcFunctions,
    uint_t                                                                               level,
    DoFType                                                                              flag ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "Interpolate" );
   // Collect all source IDs in a vector
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > > srcVertexIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >   srcEdgeIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >   srcFaceIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > >   srcCellIDs;

   for ( const VertexDoFFunction& function : srcFunctions )
   {
      srcVertexIDs.push_back( function.vertexDataID_ );
      srcEdgeIDs.push_back( function.edgeDataID_ );
      srcFaceIDs.push_back( function.faceDataID_ );
      srcCellIDs.push_back( function.cellDataID_ );
   }

   std::vector< PrimitiveID > vertexIDs = this->getStorage()->getVertexIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( vertexIDs.size() ); i++ )
   {
      Vertex& vertex = *this->getStorage()->getVertex( vertexIDs[uint_c(i)] );

      if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrovertex::interpolate( vertex, vertexDataID_, srcVertexIDs, expr, level );
      }
   }

   std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
   {
      Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c(i)] );

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroedge::interpolate< ValueType >( level, edge, edgeDataID_, srcEdgeIDs, expr );
      }
   }

   std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
   {
      Face& face = *this->getStorage()->getFace( faceIDs[uint_c(i)] );

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroface::interpolate< ValueType >( level, face, faceDataID_, srcFaceIDs, expr );
      }
   }

   std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
   {
      Cell& cell = *this->getStorage()->getCell( cellIDs[ uint_c( i ) ] );

      if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrocell::interpolate< ValueType >( level, cell, cellDataID_, srcCellIDs, expr );
      }
   }
   this->stopTiming( "Interpolate" );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::interpolate(
    const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
    const std::vector< std::reference_wrapper< const VertexDoFFunction< ValueType > > >& srcFunctions,
    uint_t                                                                               level,
    BoundaryUID                                                                          boundaryUID ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "Interpolate" );
   // Collect all source IDs in a vector
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > > srcVertexIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >   srcEdgeIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >   srcFaceIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > >   srcCellIDs;

   for ( const VertexDoFFunction& function : srcFunctions )
   {
      srcVertexIDs.push_back( function.vertexDataID_ );
      srcEdgeIDs.push_back( function.edgeDataID_ );
      srcFaceIDs.push_back( function.faceDataID_ );
      srcCellIDs.push_back( function.cellDataID_ );
   }

   std::vector< PrimitiveID > vertexIDs = this->getStorage()->getVertexIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( vertexIDs.size() ); i++ )
   {
      Vertex& vertex = *this->getStorage()->getVertex( vertexIDs[ uint_c( i ) ] );

      if ( boundaryCondition_.getBoundaryUIDFromMeshFlag( vertex.getMeshBoundaryFlag() ) == boundaryUID )
      {
         vertexdof::macrovertex::interpolate( vertex, vertexDataID_, srcVertexIDs, expr, level );
      }
   }

   std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
   {
      Edge& edge = *this->getStorage()->getEdge( edgeIDs[ uint_c( i ) ] );

      if ( boundaryCondition_.getBoundaryUIDFromMeshFlag( edge.getMeshBoundaryFlag() ) == boundaryUID )
      {
         vertexdof::macroedge::interpolate< ValueType >( level, edge, edgeDataID_, srcEdgeIDs, expr );
      }
   }

   std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
   {
      Face& face = *this->getStorage()->getFace( faceIDs[ uint_c( i ) ] );

      if ( boundaryCondition_.getBoundaryUIDFromMeshFlag( face.getMeshBoundaryFlag() ) == boundaryUID )
      {
         vertexdof::macroface::interpolate< ValueType >( level, face, faceDataID_, srcFaceIDs, expr );
      }
   }

   std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
   {
      Cell& cell = *this->getStorage()->getCell( cellIDs[ uint_c( i ) ] );

      if ( boundaryCondition_.getBoundaryUIDFromMeshFlag( cell.getMeshBoundaryFlag() ) == boundaryUID )
      {
         vertexdof::macrocell::interpolate< ValueType >( level, cell, cellDataID_, srcCellIDs, expr );
      }
   }
   this->stopTiming( "Interpolate" );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::setToZero( uint_t level ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "setToZero" );

   for ( const auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;
      vertex.getData( vertexDataID_ )->setToZero( level );
   }

   for ( const auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;
      edge.getData( edgeDataID_ )->setToZero( level );
   }

   for ( const auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;
      face.getData( faceDataID_ )->setToZero( level );
   }

   for ( const auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      cell.getData( cellDataID_ )->setToZero( level );
   }

   this->stopTiming( "setToZero" );
}

template < typename ValueType >
bool VertexDoFFunction< ValueType >::evaluate( const Point3D& coordinates,
                                               uint_t         level,
                                               ValueType&     value,
                                               real_t         searchToleranceRadius ) const
{
   WALBERLA_UNUSED( coordinates );
   WALBERLA_UNUSED( level );
   WALBERLA_UNUSED( value );
   WALBERLA_UNUSED( searchToleranceRadius );
   WALBERLA_ABORT( "VertexDoFFunction< ValueType >::evaluate not implemented for requested template parameter" );
}

template <>
bool VertexDoFFunction< real_t >::evaluate( const Point3D& coordinates,
                                            uint_t         level,
                                            real_t&        value,
                                            real_t         searchToleranceRadius ) const
{
   if ( !this->getStorage()->hasGlobalCells() )
   {
      Point2D coordinates2D( {coordinates[0], coordinates[1]} );

      for ( auto& it : this->getStorage()->getFaces() )
      {
         Face& face = *it.second;

         Point2D faceCoodinates0( {face.getCoordinates()[0][0], face.getCoordinates()[0][1]} );
         Point2D faceCoodinates1( {face.getCoordinates()[1][0], face.getCoordinates()[1][1]} );
         Point2D faceCoodinates2( {face.getCoordinates()[2][0], face.getCoordinates()[2][1]} );

         if ( isPointInTriangle( coordinates2D, faceCoodinates0, faceCoodinates1, faceCoodinates2 ) )
         {
            value = vertexdof::macroface::evaluate< real_t >( level, face, coordinates, faceDataID_ );
            return true;
         }
      }

      if ( searchToleranceRadius > 0 )
      {
         for ( auto& it : this->getStorage()->getFaces() )
         {
            Face& face = *it.second;

            Point2D faceCoodinates0( {face.getCoordinates()[0][0], face.getCoordinates()[0][1]} );
            Point2D faceCoodinates1( {face.getCoordinates()[1][0], face.getCoordinates()[1][1]} );
            Point2D faceCoodinates2( {face.getCoordinates()[2][0], face.getCoordinates()[2][1]} );

            if ( circleTriangleIntersection( coordinates2D,
                                             searchToleranceRadius,
                                             faceCoodinates0,
                                             faceCoodinates1,
                                             faceCoodinates2 ) )
            {
               value = vertexdof::macroface::evaluate< real_t >( level, face, coordinates, faceDataID_ );
               return true;
            }
         }
      }
   }
   else
   {
      for ( auto& it : this->getStorage()->getCells() )
      {
         Cell& cell = *it.second;

         if ( isPointInTetrahedron( coordinates,

                                    cell.getCoordinates()[0],
                                    cell.getCoordinates()[1],
                                    cell.getCoordinates()[2],
                                    cell.getCoordinates()[3],
                                    cell.getFaceInwardNormal( 0 ),
                                    cell.getFaceInwardNormal( 1 ),
                                    cell.getFaceInwardNormal( 2 ),
                                    cell.getFaceInwardNormal( 3 ) ) )
         {
            value = vertexdof::macrocell::evaluate< real_t >( level, cell, coordinates, cellDataID_ );
            return true;
         }
      }

      if ( searchToleranceRadius > 0 )
      {
         for ( auto& it : this->getStorage()->getCells() )
         {
            Cell& cell = *it.second;

            if ( sphereTetrahedronIntersection( coordinates,
                                                searchToleranceRadius,
                                                cell.getCoordinates()[0],
                                                cell.getCoordinates()[1],
                                                cell.getCoordinates()[2],
                                                cell.getCoordinates()[3] ) )
            {
               value = vertexdof::macrocell::evaluate< real_t >( level, cell, coordinates, cellDataID_ );
               return true;
            }
         }
      }
   }

   return false;
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::evaluateGradient( const Point3D& coordinates, uint_t level, Point3D& gradient ) const
{
   // Check if 2D or 3D function
   if ( !this->getStorage()->hasGlobalCells() )
   {
      for ( auto& it : this->getStorage()->getFaces() )
      {
         Face& face = *it.second;

         if ( sphereTriangleIntersection(
                  coordinates, 0.0, face.getCoordinates()[0], face.getCoordinates()[1], face.getCoordinates()[2] ) )
         {
            vertexdof::macroface::evaluateGradient< ValueType >( level, face, coordinates, faceDataID_, gradient );
            return;
         }
      }
   }
   else
   {
      for ( auto& it : this->getStorage()->getCells() )
      {
         Cell& cell = *it.second;

         if ( isPointInTetrahedron( coordinates,
                                    cell.getCoordinates()[0],
                                    cell.getCoordinates()[1],
                                    cell.getCoordinates()[2],
                                    cell.getCoordinates()[3] ) )
         {
            WALBERLA_ABORT( "Not implemented." )
         }
      }
   }

   WALBERLA_ABORT( "There is no local macro element including a point at the given coordinates " << coordinates )
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::swap( const VertexDoFFunction< ValueType >& other,
                                           const uint_t&                         level,
                                           const DoFType&                        flag ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "Swap" );

   for ( auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrovertex::swap< ValueType >( level, vertex, other.getVertexDataID(), vertexDataID_ );
      }
   }

   for ( auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroedge::swap< ValueType >( level, edge, other.getEdgeDataID(), edgeDataID_ );
      }
   }

   for ( auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroface::swap< ValueType >( level, face, other.getFaceDataID(), faceDataID_ );
      }
   }

   for ( auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrocell::swap< ValueType >( level, cell, other.getCellDataID(), cellDataID_ );
      }
   }

   this->stopTiming( "Swap" );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::copyFrom( const VertexDoFFunction< ValueType >& other, const uint_t& level ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "Copy" );

   for ( auto& it : this->getStorage()->getVertices() )
   {
      auto primitiveID = it.first;
      WALBERLA_ASSERT( other.getStorage()->vertexExistsLocally( primitiveID ) )
      this->getStorage()
          ->getVertex( primitiveID )
          ->getData( vertexDataID_ )
          ->copyFrom( *other.getStorage()->getVertex( primitiveID )->getData( other.getVertexDataID() ), level );
   }

   for ( auto& it : this->getStorage()->getEdges() )
   {
      auto primitiveID = it.first;
      WALBERLA_ASSERT( other.getStorage()->edgeExistsLocally( primitiveID ) )
      this->getStorage()
          ->getEdge( primitiveID )
          ->getData( edgeDataID_ )
          ->copyFrom( *other.getStorage()->getEdge( primitiveID )->getData( other.getEdgeDataID() ), level );
   }

   for ( auto& it : this->getStorage()->getFaces() )
   {
      auto primitiveID = it.first;
      WALBERLA_ASSERT( other.getStorage()->faceExistsLocally( primitiveID ) )
      this->getStorage()
          ->getFace( primitiveID )
          ->getData( faceDataID_ )
          ->copyFrom( *other.getStorage()->getFace( primitiveID )->getData( other.getFaceDataID() ), level );
   }

   for ( auto& it : this->getStorage()->getCells() )
   {
      auto primitiveID = it.first;
      WALBERLA_ASSERT( other.getStorage()->cellExistsLocally( primitiveID ) )
      this->getStorage()
          ->getCell( primitiveID )
          ->getData( cellDataID_ )
          ->copyFrom( *other.getStorage()->getCell( primitiveID )->getData( other.getCellDataID() ), level );
   }

   this->stopTiming( "Copy" );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::copyFrom( const VertexDoFFunction< ValueType >&          other,
                                               const uint_t&                                  level,
                                               const std::map< PrimitiveID::IDType, uint_t >& localPrimitiveIDsToRank,
                                               const std::map< PrimitiveID::IDType, uint_t >& otherPrimitiveIDsToRank ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "Copy" );

   walberla::mpi::BufferSystem bufferSystem( walberla::mpi::MPIManager::instance()->comm(), 9563 );
   std::set< walberla::mpi::MPIRank > receiverRanks;
   for ( auto it : localPrimitiveIDsToRank )
   {
      receiverRanks.insert( walberla::mpi::MPIRank( it.second ) );
   }
   bufferSystem.setReceiverInfo( receiverRanks, true );

   for ( auto& it : other.getStorage()->getVertices() )
   {
      PrimitiveID::IDType otherPrimitiveID = it.first;
      WALBERLA_CHECK_GREATER( otherPrimitiveIDsToRank.count( otherPrimitiveID ), 0 );
      auto otherData     = it.second->getData( other.getVertexDataID() )->getPointer( level );
      auto otherDataSize = it.second->getData( other.getVertexDataID() )->getSize( level );
      auto targetRank    = otherPrimitiveIDsToRank.at( otherPrimitiveID );
      bufferSystem.sendBuffer( targetRank ) << otherPrimitiveID;
      bufferSystem.sendBuffer( targetRank ) << uint_c( 0 );
      bufferSystem.sendBuffer( targetRank ) << otherDataSize;
      for ( uint_t i = 0; i < otherDataSize; i++ )
         bufferSystem.sendBuffer( targetRank ) << otherData[i];
   }

   for ( auto& it : other.getStorage()->getEdges() )
   {
      PrimitiveID::IDType otherPrimitiveID = it.first;
      WALBERLA_CHECK_GREATER( otherPrimitiveIDsToRank.count( otherPrimitiveID ), 0 );
      auto otherData     = it.second->getData( other.getEdgeDataID() )->getPointer( level );
      auto otherDataSize = it.second->getData( other.getEdgeDataID() )->getSize( level );
      auto targetRank    = otherPrimitiveIDsToRank.at( otherPrimitiveID );
      bufferSystem.sendBuffer( targetRank ) << otherPrimitiveID;
      bufferSystem.sendBuffer( targetRank ) << uint_c( 1 );
      bufferSystem.sendBuffer( targetRank ) << otherDataSize;
      for ( uint_t i = 0; i < otherDataSize; i++ )
         bufferSystem.sendBuffer( targetRank ) << otherData[i];
   }

   for ( auto& it : other.getStorage()->getFaces() )
   {
      PrimitiveID::IDType otherPrimitiveID = it.first;
      WALBERLA_CHECK_GREATER( otherPrimitiveIDsToRank.count( otherPrimitiveID ), 0 );
      auto otherData     = it.second->getData( other.getFaceDataID() )->getPointer( level );
      auto otherDataSize = it.second->getData( other.getFaceDataID() )->getSize( level );
      auto targetRank    = otherPrimitiveIDsToRank.at( otherPrimitiveID );
      bufferSystem.sendBuffer( targetRank ) << otherPrimitiveID;
      bufferSystem.sendBuffer( targetRank ) << uint_c( 2 );
      bufferSystem.sendBuffer( targetRank ) << otherDataSize;
      for ( uint_t i = 0; i < otherDataSize; i++ )
         bufferSystem.sendBuffer( targetRank ) << otherData[i];
   }

   for ( auto& it : other.getStorage()->getCells() )
   {
      PrimitiveID::IDType otherPrimitiveID = it.first;
      WALBERLA_CHECK_GREATER( otherPrimitiveIDsToRank.count( otherPrimitiveID ), 0 );
      auto otherData     = it.second->getData( other.getCellDataID() )->getPointer( level );
      auto otherDataSize = it.second->getData( other.getCellDataID() )->getSize( level );
      auto targetRank    = otherPrimitiveIDsToRank.at( otherPrimitiveID );
      bufferSystem.sendBuffer( targetRank ) << otherPrimitiveID;
      bufferSystem.sendBuffer( targetRank ) << uint_c( 3 );
      bufferSystem.sendBuffer( targetRank ) << otherDataSize;
      for ( uint_t i = 0; i < otherDataSize; i++ )
         bufferSystem.sendBuffer( targetRank ) << otherData[i];
   }

   bufferSystem.sendAll();

   for ( auto pkg = bufferSystem.begin(); pkg != bufferSystem.end(); ++pkg )
   {
      while ( !pkg.buffer().isEmpty() )
      {
         PrimitiveID::IDType otherID;
         uint_t              primitiveType = 4;
         uint_t              dataSize      = 0;
         ValueType           value;
         ValueType*          dstPointer;

         pkg.buffer() >> otherID;
         pkg.buffer() >> primitiveType;
         pkg.buffer() >> dataSize;

         WALBERLA_CHECK( this->getStorage()->primitiveExistsLocally( PrimitiveID( otherID ) ) );

         switch ( primitiveType )
         {
         case 0:
            dstPointer = this->getStorage()->getVertex( PrimitiveID( otherID ) )->getData( vertexDataID_ )->getPointer( level );
            break;
         case 1:
            dstPointer = this->getStorage()->getEdge( PrimitiveID( otherID ) )->getData( edgeDataID_ )->getPointer( level );
            break;
         case 2:
            dstPointer = this->getStorage()->getFace( PrimitiveID( otherID ) )->getData( faceDataID_ )->getPointer( level );
            break;
         case 3:
            dstPointer = this->getStorage()->getCell( PrimitiveID( otherID ) )->getData( cellDataID_ )->getPointer( level );
            break;
         default:
            WALBERLA_ABORT( "Invalid primitive type" )
         }

         for ( uint_t i = 0; i < dataSize; i++ )
         {
            pkg.buffer() >> value;
            dstPointer[i] = value;
         }
      }
   }

   this->stopTiming( "Copy" );
}

template < typename ValueType >
void macroFaceAssign( const uint_t&                                                              level,
                      Face&                                                                      face,
                      const std::vector< ValueType >&                                            scalars,
                      const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >& srcFaceIDs,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                dstFaceID,
                      const PrimitiveStorage& )
{
   vertexdof::macroface::assign< ValueType >( level, face, scalars, srcFaceIDs, dstFaceID );
}

template <>
void macroFaceAssign< double >( const uint_t&                                                           level,
                                Face&                                                                   face,
                                const std::vector< double >&                                            scalars,
                                const std::vector< PrimitiveDataID< FunctionMemory< double >, Face > >& srcFaceIDs,
                                const PrimitiveDataID< FunctionMemory< double >, Face >&                dstFaceID,
                                const PrimitiveStorage&                                                 storage )
{
   if ( hyteg::globalDefines::useGeneratedKernels && scalars.size() == 1 )
   {
      WALBERLA_NON_OPENMP_SECTION() { storage.getTimingTree()->start( "1 RHS function" ); }
      auto dstData = face.getData( dstFaceID )->getPointer( level );
      auto srcData = face.getData( srcFaceIDs.at( 0 ) )->getPointer( level );
      auto scalar  = scalars.at( 0 );
      vertexdof::macroface::generated::assign_2D_macroface_vertexdof_1_rhs_function(
          dstData, srcData, scalar, static_cast< int32_t >( level ) );
      WALBERLA_NON_OPENMP_SECTION() { storage.getTimingTree()->stop( "1 RHS function" ); }
   }
   else if ( hyteg::globalDefines::useGeneratedKernels && scalars.size() == 2 )
   {
      WALBERLA_NON_OPENMP_SECTION() { storage.getTimingTree()->start( "2 RHS functions" ); }
      auto dstData  = face.getData( dstFaceID )->getPointer( level );
      auto srcData0 = face.getData( srcFaceIDs.at( 0 ) )->getPointer( level );
      auto srcData1 = face.getData( srcFaceIDs.at( 1 ) )->getPointer( level );
      auto scalar0  = scalars.at( 0 );
      auto scalar1  = scalars.at( 1 );
      vertexdof::macroface::generated::assign_2D_macroface_vertexdof_2_rhs_functions(
          dstData, srcData0, srcData1, scalar0, scalar1, static_cast< int32_t >( level ) );
         WALBERLA_NON_OPENMP_SECTION() { storage.getTimingTree()->stop( "2 RHS functions" ); }
   }
   else if ( hyteg::globalDefines::useGeneratedKernels && scalars.size() == 3 )
   {
      WALBERLA_NON_OPENMP_SECTION() { storage.getTimingTree()->start( "3 RHS functions" ); }
      auto dstData  = face.getData( dstFaceID )->getPointer( level );
      auto srcData0 = face.getData( srcFaceIDs.at( 0 ) )->getPointer( level );
      auto srcData1 = face.getData( srcFaceIDs.at( 1 ) )->getPointer( level );
      auto srcData2 = face.getData( srcFaceIDs.at( 2 ) )->getPointer( level );
      auto scalar0  = scalars.at( 0 );
      auto scalar1  = scalars.at( 1 );
      auto scalar2  = scalars.at( 2 );
      vertexdof::macroface::generated::assign_2D_macroface_vertexdof_3_rhs_functions(
          dstData, srcData0, srcData1, srcData2, scalar0, scalar1, scalar2, static_cast< int32_t >( level ) );
         WALBERLA_NON_OPENMP_SECTION() { storage.getTimingTree()->stop( "3 RHS functions" ); }
   }
   else
   {
      vertexdof::macroface::assign< double >( level, face, scalars, srcFaceIDs, dstFaceID );
   }
}

template < typename ValueType >
void macroCellAssign( const uint_t&                                                              level,
                      Cell&                                                                      cell,
                      const std::vector< ValueType >&                                            scalars,
                      const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > >& srcCellIDs,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Cell >&                dstCellID )
{
   vertexdof::macrocell::assign< ValueType >( level, cell, scalars, srcCellIDs, dstCellID );
}

template <>
void macroCellAssign< double >( const uint_t&                                                           level,
                                Cell&                                                                   cell,
                                const std::vector< double >&                                            scalars,
                                const std::vector< PrimitiveDataID< FunctionMemory< double >, Cell > >& srcCellIDs,
                                const PrimitiveDataID< FunctionMemory< double >, Cell >&                dstCellID )
{
   if ( hyteg::globalDefines::useGeneratedKernels && scalars.size() == 1 )
   {
      auto dstData = cell.getData( dstCellID )->getPointer( level );
      auto srcData = cell.getData( srcCellIDs.at( 0 ) )->getPointer( level );
      auto scalar  = scalars.at( 0 );
      vertexdof::macrocell::generated::assign_3D_macrocell_vertexdof_1_rhs_function(
          dstData, srcData, scalar, static_cast< int32_t >( level ) );
   }
   else if ( hyteg::globalDefines::useGeneratedKernels && scalars.size() == 2 )
   {
      auto dstData  = cell.getData( dstCellID )->getPointer( level );
      auto srcData0 = cell.getData( srcCellIDs.at( 0 ) )->getPointer( level );
      auto srcData1 = cell.getData( srcCellIDs.at( 1 ) )->getPointer( level );
      auto scalar0  = scalars.at( 0 );
      auto scalar1  = scalars.at( 1 );
      vertexdof::macrocell::generated::assign_3D_macrocell_vertexdof_2_rhs_functions(
          dstData, srcData0, srcData1, scalar0, scalar1, static_cast< int32_t >( level ) );
   }
   else if ( hyteg::globalDefines::useGeneratedKernels && scalars.size() == 3 )
   {
      auto dstData  = cell.getData( dstCellID )->getPointer( level );
      auto srcData0 = cell.getData( srcCellIDs.at( 0 ) )->getPointer( level );
      auto srcData1 = cell.getData( srcCellIDs.at( 1 ) )->getPointer( level );
      auto srcData2 = cell.getData( srcCellIDs.at( 2 ) )->getPointer( level );
      auto scalar0  = scalars.at( 0 );
      auto scalar1  = scalars.at( 1 );
      auto scalar2  = scalars.at( 2 );
      vertexdof::macrocell::generated::assign_3D_macrocell_vertexdof_3_rhs_functions(
          dstData, srcData0, srcData1, srcData2, scalar0, scalar1, scalar2, static_cast< int32_t >( level ) );
   }
   else
   {
      vertexdof::macrocell::assign< double >( level, cell, scalars, srcCellIDs, dstCellID );
   }
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::assign(
    const std::vector< ValueType >&                                                      scalars,
    const std::vector< std::reference_wrapper< const VertexDoFFunction< ValueType > > >& functions,
    size_t                                                                               level,
    DoFType                                                                              flag ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "Assign" );

   WALBERLA_ASSERT_EQUAL( scalars.size(), functions.size() )

   // Collect all source IDs in a vector
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > > srcVertexIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >   srcEdgeIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >   srcFaceIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > >   srcCellIDs;

   for ( const VertexDoFFunction< ValueType >& function : functions )
   {
      srcVertexIDs.push_back( function.vertexDataID_ );
      srcEdgeIDs.push_back( function.edgeDataID_ );
      srcFaceIDs.push_back( function.faceDataID_ );
      srcCellIDs.push_back( function.cellDataID_ );
   }
   this->getStorage()->getTimingTree()->start( "Vertex" );
   
   std::vector< PrimitiveID > vertexIDs = this->getStorage()->getVertexIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( vertexIDs.size() ); i++ )
   {
      Vertex& vertex = *this->getStorage()->getVertex( vertexIDs[uint_c(i)] );

      if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrovertex::assign< ValueType >( vertex, scalars, srcVertexIDs, vertexDataID_, level );
      }
   }
   
   this->getStorage()->getTimingTree()->stop( "Vertex" );
   this->getStorage()->getTimingTree()->start( "Edge" );
   
   std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
   {
      Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c(i)] );

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroedge::assign< ValueType >( level, edge, scalars, srcEdgeIDs, edgeDataID_ );
      }
   }
   
   this->getStorage()->getTimingTree()->stop( "Edge" );
   this->getStorage()->getTimingTree()->start( "Face" );

   std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
   {
      Face& face = *this->getStorage()->getFace( faceIDs[uint_c(i)] );

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         macroFaceAssign< ValueType >( level, face, scalars, srcFaceIDs, faceDataID_, *this->getStorage() );
      }
   }

   this->getStorage()->getTimingTree()->stop( "Face" );
   this->getStorage()->getTimingTree()->start( "Cell" );

   std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
   {
      Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c(i)] );
      if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         macroCellAssign< ValueType >( level, cell, scalars, srcCellIDs, cellDataID_ );
      }
   }
   this->getStorage()->getTimingTree()->stop( "Cell" );
   this->stopTiming( "Assign" );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::add( ValueType scalar, uint_t level, DoFType flag ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "Add" );

   std::vector< PrimitiveID > vertexIDs = this->getStorage()->getVertexIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( vertexIDs.size() ); i++ )
   {
      Vertex& vertex = *this->getStorage()->getVertex( vertexIDs[uint_c(i)] );

      if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrovertex::add< ValueType >( vertex, scalar, vertexDataID_, level );
      }
   }

   std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
   {
      Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c(i)] );

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroedge::add< ValueType >( level, edge, scalar, edgeDataID_ );
      }
   }

   std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
   {
      Face& face = *this->getStorage()->getFace( faceIDs[uint_c(i)] );

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroface::add< ValueType >( level, face, scalar, faceDataID_ );
      }
   }

   std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
   {
      Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c(i)] );
      if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrocell::add< ValueType >( level, cell, scalar, cellDataID_ );
      }
   }

   this->stopTiming( "Add" );
}

template < typename ValueType >
void macroFaceAdd( const uint_t&                                                              level,
                   Face&                                                                      face,
                   const std::vector< ValueType >&                                            scalars,
                   const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >& srcFaceIDs,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                dstFaceID,
                   const PrimitiveStorage& )
{
   vertexdof::macroface::add< ValueType >( level, face, scalars, srcFaceIDs, dstFaceID );
}

template <>
void macroFaceAdd< double >( const uint_t&                                                           level,
                             Face&                                                                   face,
                             const std::vector< double >&                                            scalars,
                             const std::vector< PrimitiveDataID< FunctionMemory< double >, Face > >& srcFaceIDs,
                             const PrimitiveDataID< FunctionMemory< double >, Face >&                dstFaceID,
                             const PrimitiveStorage&                                                 storage )
{
   if ( hyteg::globalDefines::useGeneratedKernels && scalars.size() == 1 )
   {
      WALBERLA_NON_OPENMP_SECTION() { storage.getTimingTree()->start( "1 RHS function" ); }
      auto dstData = face.getData( dstFaceID )->getPointer( level );
      auto srcData = face.getData( srcFaceIDs.at( 0 ) )->getPointer( level );
      auto scalar  = scalars.at( 0 );
      vertexdof::macroface::generated::add_2D_macroface_vertexdof_1_rhs_function(
          dstData, srcData, scalar, static_cast< int32_t >( level ) );
      WALBERLA_NON_OPENMP_SECTION() { storage.getTimingTree()->stop( "1 RHS function" ); }
   }
   else if ( hyteg::globalDefines::useGeneratedKernels && scalars.size() == 2 )
   {
      WALBERLA_NON_OPENMP_SECTION() { storage.getTimingTree()->start( "2 RHS functions" ); }
      auto dstData  = face.getData( dstFaceID )->getPointer( level );
      auto srcData0 = face.getData( srcFaceIDs.at( 0 ) )->getPointer( level );
      auto srcData1 = face.getData( srcFaceIDs.at( 1 ) )->getPointer( level );
      auto scalar0  = scalars.at( 0 );
      auto scalar1  = scalars.at( 1 );
      vertexdof::macroface::generated::add_2D_macroface_vertexdof_2_rhs_functions(
          dstData, srcData0, srcData1, scalar0, scalar1, static_cast< int32_t >( level ) );
      WALBERLA_NON_OPENMP_SECTION() { storage.getTimingTree()->stop( "2 RHS functions" ); }
   }
   else if ( hyteg::globalDefines::useGeneratedKernels && scalars.size() == 3 )
   {
      WALBERLA_NON_OPENMP_SECTION() { storage.getTimingTree()->start( "3 RHS functions" ); }
      auto dstData  = face.getData( dstFaceID )->getPointer( level );
      auto srcData0 = face.getData( srcFaceIDs.at( 0 ) )->getPointer( level );
      auto srcData1 = face.getData( srcFaceIDs.at( 1 ) )->getPointer( level );
      auto srcData2 = face.getData( srcFaceIDs.at( 2 ) )->getPointer( level );
      auto scalar0  = scalars.at( 0 );
      auto scalar1  = scalars.at( 1 );
      auto scalar2  = scalars.at( 2 );
      vertexdof::macroface::generated::add_2D_macroface_vertexdof_3_rhs_functions(
          dstData, srcData0, srcData1, srcData2, scalar0, scalar1, scalar2, static_cast< int32_t >( level ) );
      WALBERLA_NON_OPENMP_SECTION() { storage.getTimingTree()->stop( "3 RHS functions" ); }
   }
   else
   {
      vertexdof::macroface::add< double >( level, face, scalars, srcFaceIDs, dstFaceID );
   }
}

template < typename ValueType >
void macroCellAdd( const uint_t&                                                              level,
                   Cell&                                                                      cell,
                   const std::vector< ValueType >&                                            scalars,
                   const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > >& srcCellIDs,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Cell >&                dstCellID )
{
   vertexdof::macrocell::add< ValueType >( level, cell, scalars, srcCellIDs, dstCellID );
}

template <>
void macroCellAdd< double >( const uint_t&                                                           level,
                             Cell&                                                                   cell,
                             const std::vector< double >&                                            scalars,
                             const std::vector< PrimitiveDataID< FunctionMemory< double >, Cell > >& srcCellIDs,
                             const PrimitiveDataID< FunctionMemory< double >, Cell >&                dstCellID )
{
   if ( hyteg::globalDefines::useGeneratedKernels && scalars.size() == 1 )
   {
      auto dstData = cell.getData( dstCellID )->getPointer( level );
      auto srcData = cell.getData( srcCellIDs.at( 0 ) )->getPointer( level );
      auto scalar  = scalars.at( 0 );
      vertexdof::macrocell::generated::add_3D_macrocell_vertexdof_1_rhs_function(
          dstData, srcData, scalar, static_cast< int32_t >( level ) );
   }
   else
   {
      vertexdof::macrocell::add< double >( level, cell, scalars, srcCellIDs, dstCellID );
   }
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::add(
    const std::vector< ValueType >&                                                      scalars,
    const std::vector< std::reference_wrapper< const VertexDoFFunction< ValueType > > >& functions,
    size_t                                                                               level,
    DoFType                                                                              flag ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "Add" );
   // Collect all source IDs in a vector
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > > srcVertexIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >   srcEdgeIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >   srcFaceIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > >   srcCellIDs;

   for ( const VertexDoFFunction< ValueType >& function : functions )
   {
      srcVertexIDs.push_back( function.vertexDataID_ );
      srcEdgeIDs.push_back( function.edgeDataID_ );
      srcFaceIDs.push_back( function.faceDataID_ );
      srcCellIDs.push_back( function.cellDataID_ );
   }

   std::vector< PrimitiveID > vertexIDs = this->getStorage()->getVertexIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( vertexIDs.size() ); i++ )
   {
      Vertex& vertex = *this->getStorage()->getVertex( vertexIDs[uint_c(i)] );

      if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrovertex::add( vertex, scalars, srcVertexIDs, vertexDataID_, level );
      }
   }

   std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
   {
      Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c(i)] );

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroedge::add< ValueType >( level, edge, scalars, srcEdgeIDs, edgeDataID_ );
      }
   }

   std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
   {
      Face& face = *this->getStorage()->getFace( faceIDs[uint_c(i)] );

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         macroFaceAdd< ValueType >( level, face, scalars, srcFaceIDs, faceDataID_, *this->getStorage() );
      }
   }

   std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
   {
      Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c(i)] );
      if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         macroCellAdd< ValueType >( level, cell, scalars, srcCellIDs, cellDataID_ );
      }
   }
   this->stopTiming( "Add" );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::multElementwise(
    const std::vector< std::reference_wrapper< const VertexDoFFunction< ValueType > > >& functions,
    const uint_t                                                                         level,
    const DoFType                                                                        flag ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "Multiply elementwise" );
   // Collect all source IDs in a vector
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > > srcVertexIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >   srcEdgeIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >   srcFaceIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > >   srcCellIDs;

   for ( const VertexDoFFunction& function : functions )
   {
      srcVertexIDs.push_back( function.vertexDataID_ );
      srcEdgeIDs.push_back( function.edgeDataID_ );
      srcFaceIDs.push_back( function.faceDataID_ );
      srcCellIDs.push_back( function.cellDataID_ );
   }

   std::vector< PrimitiveID > vertexIDs = this->getStorage()->getVertexIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( vertexIDs.size() ); i++ )
   {
      Vertex& vertex = *this->getStorage()->getVertex( vertexIDs[uint_c(i)] );

      if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrovertex::multElementwise< ValueType >( level, vertex, srcVertexIDs, vertexDataID_ );
      }
   }

   std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
   {
      Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c(i)] );

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroedge::multElementwise< ValueType >( level, edge, srcEdgeIDs, edgeDataID_ );
      }
   }

   std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
   {
      Face& face = *this->getStorage()->getFace( faceIDs[uint_c(i)] );

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroface::multElementwise< ValueType >( level, face, srcFaceIDs, faceDataID_ );
      }
   }

   std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
   {
      Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c(i)] );
      
      if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrocell::multElementwise< ValueType >( level, cell, srcCellIDs, cellDataID_ );
      }
   }
   this->stopTiming( "Multiply elementwise" );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::invertElementwise( uint_t level, DoFType flag, bool workOnHalos ) const
{
   WALBERLA_UNUSED( level );
   WALBERLA_UNUSED( flag );
   WALBERLA_UNUSED( workOnHalos );
   WALBERLA_ABORT( "VertexDoFFunction< ValueType >::invertElementwise not available for requested ValueType" );
}

template < typename ValueType >
ValueType VertexDoFFunction< ValueType >::dotGlobal( const VertexDoFFunction< ValueType >& rhs, size_t level, DoFType flag ) const
{
   ValueType scalarProduct = dotLocal( rhs, level, flag );
   this->startTiming( "Dot (reduce)" );
   walberla::mpi::allReduceInplace( scalarProduct, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
   this->stopTiming( "Dot (reduce)" );
   return scalarProduct;
}

template < typename ValueType >
ValueType VertexDoFFunction< ValueType >::dotLocal( const VertexDoFFunction< ValueType >& rhs, size_t level, DoFType flag ) const
{
   if ( isDummy() )
   {
      return ValueType( 0 );
   }
   this->startTiming( "Dot (local)" );
   auto scalarProduct = ValueType( 0 );

   ValueType scalarProductVertices = 0;
   std::vector< PrimitiveID > vertexIDs = this->getStorage()->getVertexIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for reduction(+: scalarProductVertices)
   #endif
   for ( int i = 0; i < int_c( vertexIDs.size() ); i++ )
   {
      Vertex& vertex = *this->getStorage()->getVertex( vertexIDs[ uint_c( i ) ] );

      if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         scalarProductVertices += vertexdof::macrovertex::dot( vertex, vertexDataID_, rhs.vertexDataID_, level );
      }
   }
   scalarProduct += scalarProductVertices;

   if ( level >= 1 )
   {
      ValueType scalarProductEdges = 0;
      std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
      #ifdef WALBERLA_BUILD_WITH_OPENMP
      #pragma omp parallel for reduction(+: scalarProductEdges)
      #endif
      for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
      {
         Edge& edge = *this->getStorage()->getEdge( edgeIDs[ uint_c( i ) ] );

         if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
         {
            scalarProductEdges += vertexdof::macroedge::dot< ValueType >( level, edge, edgeDataID_, rhs.edgeDataID_ );
         }
      }
      scalarProduct += scalarProductEdges;

      ValueType scalarProductFaces = 0;
      std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
      #ifdef WALBERLA_BUILD_WITH_OPENMP
      #pragma omp parallel for reduction(+: scalarProductFaces)
      #endif
      for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
      {
         Face& face = *this->getStorage()->getFace( faceIDs[ uint_c( i ) ] );

         if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
         {
            scalarProductFaces += vertexdof::macroface::dot< ValueType >( level, face, faceDataID_, rhs.faceDataID_ );
         }
      }
      scalarProduct += scalarProductFaces;

      ValueType scalarProductCells = 0;
      std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
      #ifdef WALBERLA_BUILD_WITH_OPENMP
      #pragma omp parallel for reduction(+: scalarProductCells)
      #endif
      for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
      {
         Cell& cell = *this->getStorage()->getCell( cellIDs[ uint_c( i ) ] );

         if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
         {
            scalarProductCells += vertexdof::macrocell::dot< ValueType >( level, cell, cellDataID_, rhs.cellDataID_ );
         }
      }
      scalarProduct += scalarProductCells;
   }
   this->stopTiming( "Dot (local)" );
   return scalarProduct;
}

template < typename ValueType >
ValueType VertexDoFFunction< ValueType >::sumGlobal( const uint_t& level, const DoFType& flag, const bool& absolute ) const
{
   ValueType sum = sumLocal( level, flag, absolute );
   this->startTiming( "Sum (reduce)" );
   walberla::mpi::allReduceInplace( sum, walberla::mpi::SUM, this->getStorage()->getSplitCommunicatorByPrimitiveDistribution() );
   this->stopTiming( "Sum (reduce)" );
   return sum;
}

template < typename ValueType >
ValueType VertexDoFFunction< ValueType >::sumLocal( const uint_t& level, const DoFType& flag, const bool& absolute ) const
{
   if ( isDummy() )
   {
      return ValueType( 0 );
   }
   this->startTiming( "Sum (local)" );
   auto sum = ValueType( 0 );

   for ( const auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         sum += vertexdof::macrovertex::sum( level, vertex, vertexDataID_, absolute );
      }
   }

   for ( const auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         sum += vertexdof::macroedge::sum< ValueType >( level, edge, edgeDataID_, absolute );
      }
   }

   for ( const auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         sum += vertexdof::macroface::sum< ValueType >( level, face, faceDataID_, absolute );
      }
   }

   for ( const auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         sum += vertexdof::macrocell::sum< ValueType >( level, cell, cellDataID_, absolute );
      }
   }
   this->stopTiming( "Sum (local)" );
   return sum;
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::enumerate( uint_t level ) const
{
   if ( isDummy() )
   {
      return;
   }

   this->startTiming( "Enumerate" );

   uint_t counter = hyteg::numberOfLocalDoFs< VertexDoFFunctionTag >( *( this->getStorage() ), level );

   std::vector< uint_t > dofs_per_rank = walberla::mpi::allGather( counter );

   auto startOnRank = ValueType( 0 );

   for ( uint_t i = 0; i < uint_c( walberla::MPIManager::instance()->rank() ); ++i )
   {
      startOnRank += static_cast< ValueType >( dofs_per_rank[i] );
   }
   enumerate( level, startOnRank );
   this->stopTiming( "Enumerate" );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::enumerate( uint_t level, ValueType& offset ) const
{
   if ( isDummy() )
   {
      return;
   }

   for ( auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;
      vertexdof::macrovertex::enumerate( level, vertex, vertexDataID_, offset );
   }

   for ( auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;
      vertexdof::macroedge::enumerate< ValueType >( level, edge, edgeDataID_, offset );
   }

   if ( level >= 2 )
   {
      for ( auto& it : this->getStorage()->getFaces() )
      {
         Face& face = *it.second;
         vertexdof::macroface::enumerate< ValueType >( level, face, faceDataID_, offset );
      }
   }

   if ( level >= 2 )
   {
      for ( auto& it : this->getStorage()->getCells() )
      {
         Cell& cell = *it.second;
         vertexdof::macrocell::enumerate< ValueType >( level, cell, cellDataID_, offset );
      }
   }

   /// in contrast to other methods in the function class enumerate needs to communicate due to its usage in the PETSc solvers
   communication::syncFunctionBetweenPrimitives( *this, level );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::integrateDG( FaceDoFFunction< ValueType >&        rhs,
                                                  VertexDoFFunction< ValueType >& rhsP1,
                                                  uint_t                          level,
                                                  DoFType                         flag )
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "integrateDG" );

   rhsP1.startCommunication< Edge, Vertex >( level );
   rhsP1.startCommunication< Face, Edge >( level );

   rhs.template startCommunication< Face, Edge >( level );
   rhs.template endCommunication< Face, Edge >( level );

   rhs.template startCommunication< Edge, Vertex >( level );
   rhs.template endCommunication< Edge, Vertex >( level );

   rhsP1.endCommunication< Edge, Vertex >( level );

   for ( auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrovertex::integrateDG< ValueType >(
             vertex, this->getStorage(), rhs.getVertexDataID(), rhsP1.getVertexDataID(), vertexDataID_, level );
      }
   }

   communicators_[level]->template startCommunication< Vertex, Edge >();
   rhsP1.endCommunication< Face, Edge >( level );

   for ( auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroedge::integrateDG< ValueType >(
             level, edge, this->getStorage(), rhs.getEdgeDataID(), rhsP1.getEdgeDataID(), edgeDataID_ );
      }
   }

   communicators_[level]->template endCommunication< Vertex, Edge >();
   communicators_[level]->template startCommunication< Edge, Face >();

   for ( auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroface::integrateDG< ValueType >( level, face, rhs.getFaceDataID(), rhsP1.getFaceDataID(), faceDataID_ );
      }
   }

   communicators_[level]->template endCommunication< Edge, Face >();

   this->stopTiming( "integrateDG" );
}

template < typename ValueType >
ValueType VertexDoFFunction< ValueType >::getMaxValue( uint_t level, DoFType flag, bool mpiReduce ) const
{
   if ( isDummy() )
   {
      return ValueType( 0 );
   }
   auto localMax = -std::numeric_limits< ValueType >::max();

   for ( auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      localMax   = std::max( localMax, vertexdof::macrocell::getMaxValue< ValueType >( level, cell, cellDataID_ ) );
   }

   for ( auto& it : this->getStorage()->getFaces() )
   {
      Face&         face   = *it.second;
      const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         localMax = std::max( localMax, vertexdof::macroface::getMaxValue< ValueType >( level, face, faceDataID_ ) );
      }
   }

   for ( auto& it : this->getStorage()->getEdges() )
   {
      Edge&         edge   = *it.second;
      const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         localMax = std::max( localMax, vertexdof::macroedge::getMaxValue< ValueType >( level, edge, edgeDataID_ ) );
      }
   }

   for ( auto& it : this->getStorage()->getVertices() )
   {
      Vertex&       vertex   = *it.second;
      const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         localMax = std::max( localMax, vertexdof::macrovertex::getMaxValue< ValueType >( level, vertex, vertexDataID_ ) );
      }
   }

   ValueType globalMax = localMax;
   if ( mpiReduce )
   {
      globalMax = walberla::mpi::allReduce( localMax, walberla::mpi::MAX );
   }

   return globalMax;
}

template < typename ValueType >
ValueType VertexDoFFunction< ValueType >::getMaxMagnitude( uint_t level, DoFType flag, bool mpiReduce ) const
{
   if ( isDummy() )
   {
      return ValueType( 0 );
   }
   auto localMax = ValueType( 0.0 );

   for ( auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      localMax   = std::max( localMax, vertexdof::macrocell::getMaxMagnitude< ValueType >( level, cell, cellDataID_ ) );
   }

   for ( auto& it : this->getStorage()->getFaces() )
   {
      Face&         face   = *it.second;
      const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         localMax = std::max( localMax, vertexdof::macroface::getMaxMagnitude< ValueType >( level, face, faceDataID_ ) );
      }
   }

   for ( auto& it : this->getStorage()->getEdges() )
   {
      Edge&         edge   = *it.second;
      const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         localMax = std::max( localMax, vertexdof::macroedge::getMaxMagnitude< ValueType >( level, edge, edgeDataID_ ) );
      }
   }

   for ( auto& it : this->getStorage()->getVertices() )
   {
      Vertex&       vertex   = *it.second;
      const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         localMax = std::max( localMax, vertexdof::macrovertex::getMaxMagnitude< ValueType >( level, vertex, vertexDataID_ ) );
      }
   }

   ValueType globalMax = localMax;
   if ( mpiReduce )
   {
      globalMax = walberla::mpi::allReduce( localMax, walberla::mpi::MAX );
   }

   return globalMax;
}

template < typename ValueType >
ValueType VertexDoFFunction< ValueType >::getMinValue( uint_t level, DoFType flag, bool mpiReduce ) const
{
   if ( isDummy() )
   {
      return ValueType( 0 );
   }
   auto localMin = std::numeric_limits< ValueType >::max();

   for ( auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      localMin   = std::min( localMin, vertexdof::macrocell::getMinValue< ValueType >( level, cell, cellDataID_ ) );
   }

   for ( auto& it : this->getStorage()->getFaces() )
   {
      Face&         face   = *it.second;
      const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         localMin = std::min( localMin, vertexdof::macroface::getMinValue< ValueType >( level, face, faceDataID_ ) );
      }
   }

   for ( auto& it : this->getStorage()->getEdges() )
   {
      Edge&         edge   = *it.second;
      const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         localMin = std::min( localMin, vertexdof::macroedge::getMinValue< ValueType >( level, edge, edgeDataID_ ) );
      }
   }

   for ( auto& it : this->getStorage()->getVertices() )
   {
      Vertex&       vertex   = *it.second;
      const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         localMin = std::min( localMin, vertexdof::macrovertex::getMinValue< ValueType >( level, vertex, vertexDataID_ ) );
      }
   }

   ValueType globalMin = localMin;
   if ( mpiReduce )
   {
      globalMin = -walberla::mpi::allReduce( -localMin, walberla::mpi::MAX );
   }

   return globalMin;
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::setLocalCommunicationMode(
    const communication::BufferedCommunicator::LocalCommunicationMode& localCommunicationMode )
{
   if ( isDummy() )
   {
      return;
   }
   for ( auto& communicator : communicators_ )
   {
      communicator.second->setLocalCommunicationMode( localCommunicationMode );
   }
   for ( auto& communicator : additiveCommunicators_ )
   {
      communicator.second->setLocalCommunicationMode( localCommunicationMode );
   }
}

template < typename ValueType >
template < typename PrimitiveType >
void VertexDoFFunction< ValueType >::interpolateByPrimitiveType( const ValueType& constant, uint_t level, DoFType flag ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "Interpolate" );

   if ( std::is_same< PrimitiveType, Vertex >::value )
   {
      std::vector< PrimitiveID > vertexIDs = this->getStorage()->getVertexIDs();
      #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
      for ( int i = 0; i < int_c( vertexIDs.size() ); i++ )
      {
         Vertex& vertex = *this->getStorage()->getVertex( vertexIDs[uint_c(i)] );

         if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
         {
            vertexdof::macrovertex::interpolate( level, vertex, vertexDataID_, constant );
         }
      }
   }
   else if ( std::is_same< PrimitiveType, Edge >::value )
   {
      std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
      #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
      for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
      {
         Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c(i)] );

         if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
         {
            vertexdof::macroedge::interpolate( level, edge, edgeDataID_, constant );
         }
      }
   }
   else if ( std::is_same< PrimitiveType, Face >::value )
   {
      std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
      #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
      for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
      {
         Face& face = *this->getStorage()->getFace( faceIDs[uint_c(i)] );

         if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
         {
            vertexdof::macroface::interpolate( level, face, faceDataID_, constant );
         }
      }
   }
   else if ( std::is_same< PrimitiveType, Cell >::value )
   {
      std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
      #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
      for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
      {
         Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c(i)] );

         if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
         {
            vertexdof::macrocell::interpolate( level, cell, cellDataID_, constant );
         }
      }
   }

   this->stopTiming( "Interpolate" );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::toVector( const VertexDoFFunction< idx_t >&     numerator,
                                               const std::shared_ptr< VectorProxy >& vec,
                                               uint_t                                level,
                                               DoFType                               flag ) const
{
   if constexpr ( !std::is_same< ValueType, real_t >::value )
   {
      WALBERLA_UNUSED( numerator );
      WALBERLA_UNUSED( vec );
      WALBERLA_UNUSED( level );
      WALBERLA_UNUSED( flag );
      WALBERLA_ABORT( "VertexDoFFunction< T >::toVector() not implemented for T = " << typeid( ValueType ).name() );
   }
   else
   {
      for ( auto& it : this->getStorage()->getVertices() )
      {
         Vertex& vertex = *it.second;

         const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
         if ( testFlag( vertexBC, flag ) )
         {
            vertexdof::macrovertex::createVectorFromFunction< real_t >(
                vertex, this->getVertexDataID(), numerator.getVertexDataID(), vec, level );
         }
      }

      for ( auto& it : this->getStorage()->getEdges() )
      {
         Edge& edge = *it.second;

         const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if ( testFlag( edgeBC, flag ) )
         {
            vertexdof::macroedge::createVectorFromFunction< real_t >(
                level, edge, this->getEdgeDataID(), numerator.getEdgeDataID(), vec );
         }
      }

      for ( auto& it : this->getStorage()->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if ( testFlag( faceBC, flag ) )
         {
            vertexdof::macroface::createVectorFromFunction< real_t >(
                level, face, this->getFaceDataID(), numerator.getFaceDataID(), vec );
         }
      }

      for ( auto& it : this->getStorage()->getCells() )
      {
         Cell& cell = *it.second;

         const DoFType cellBC = this->getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
         if ( testFlag( cellBC, flag ) )
         {
            vertexdof::macrocell::createVectorFromFunction< real_t >(
                level, cell, this->getCellDataID(), numerator.getCellDataID(), vec );
         }
      }
   }
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::fromVector( const VertexDoFFunction< idx_t >&     numerator,
                                                 const std::shared_ptr< VectorProxy >& vec,
                                                 uint_t                                level,
                                                 DoFType                               flag ) const
{
   if constexpr ( !std::is_same< ValueType, real_t >::value )
   {
      WALBERLA_UNUSED( numerator );
      WALBERLA_UNUSED( vec );
      WALBERLA_UNUSED( level );
      WALBERLA_UNUSED( flag );
      WALBERLA_ABORT( "VertexDoFFunction< T >::fromVector() not implemented for T = " << typeid( ValueType ).name() );
   }
   else
   {
      for ( auto& it : this->getStorage()->getVertices() )
      {
         Vertex& vertex = *it.second;

         const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
         if ( testFlag( vertexBC, flag ) )
         {
            vertexdof::macrovertex::createFunctionFromVector< real_t >(
                vertex, this->getVertexDataID(), numerator.getVertexDataID(), vec, level );
         }
      }

      this->startCommunication< Vertex, Edge >( level );
      this->endCommunication< Vertex, Edge >( level );

      for ( auto& it : this->getStorage()->getEdges() )
      {
         Edge& edge = *it.second;

         const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if ( testFlag( edgeBC, flag ) )
         {
            vertexdof::macroedge::createFunctionFromVector< real_t >(
                level, edge, this->getEdgeDataID(), numerator.getEdgeDataID(), vec );
         }
      }

      this->startCommunication< Edge, Face >( level );
      this->endCommunication< Edge, Face >( level );

      for ( auto& it : this->getStorage()->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if ( testFlag( faceBC, flag ) )
         {
            vertexdof::macroface::createFunctionFromVector< real_t >(
                level, face, this->getFaceDataID(), numerator.getFaceDataID(), vec );
         }
      }

      this->startCommunication< Face, Cell >( level );
      this->endCommunication< Face, Cell >( level );

      for ( auto& it : this->getStorage()->getCells() )
      {
         Cell& cell = *it.second;

         const DoFType cellBC = this->getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
         if ( testFlag( cellBC, flag ) )
         {
            vertexdof::macrocell::createFunctionFromVector< real_t >(
                level, cell, this->getCellDataID(), numerator.getCellDataID(), vec );
         }
      }
   }
}

// =================
//  specialisations
// =================
template <>
void VertexDoFFunction< real_t >::invertElementwise( uint_t level, DoFType flag, bool workOnHalos ) const
{
   if ( isDummy() )
   {
      return;
   }

   this->startTiming( "Invert elementwise" );

   if ( workOnHalos )
   {
      for ( const auto& it : this->getStorage()->getVertices() )
      {
         Vertex& vertex = *it.second;

         if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
         {
            real_t* data = vertex.getData( vertexDataID_ )->getPointer( level );
            uint_t  size = vertex.getData( vertexDataID_ )->getSize( level );
            for ( uint_t k = 0; k < size; ++k )
            {
               data[k] = real_c( 1.0 ) / data[k];
               // data[0]      = real_c( 1.0 ) / data[0];
            }
         }
      }

      for ( const auto& it : this->getStorage()->getEdges() )
      {
         Edge& edge = *it.second;

         if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
         {
            real_t* data = edge.getData( edgeDataID_ )->getPointer( level );
            uint_t  size = edge.getData( edgeDataID_ )->getSize( level );
            for ( uint_t k = 0; k < size; ++k )
            {
               data[k] = real_c( 1.0 ) / data[k];
            }
         }
      }

      for ( const auto& it : this->getStorage()->getFaces() )
      {
         Face& face = *it.second;

         if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
         {
            real_t* data = face.getData( faceDataID_ )->getPointer( level );
            uint_t  size = face.getData( faceDataID_ )->getSize( level );
            for ( uint_t k = 0; k < size; ++k )
            {
               data[k] = real_c( 1.0 ) / data[k];
            }
         }
      }

      for ( const auto& it : this->getStorage()->getCells() )
      {
         Cell& cell = *it.second;

         if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
         {
            real_t* data = cell.getData( cellDataID_ )->getPointer( level );
            uint_t  size = cell.getData( cellDataID_ )->getSize( level );
            for ( uint_t k = 0; k < size; ++k )
            {
               data[k] = real_c( 1.0 ) / data[k];
            }
         }
      }
   }

   // do not work on halos
   else
   {
      for ( const auto& it : this->getStorage()->getVertices() )
      {
         Vertex& vertex = *it.second;

         if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
         {
            real_t* data = vertex.getData( vertexDataID_ )->getPointer( level );
            data[0]      = real_c( 1.0 ) / data[0];
         }
      }

      for ( const auto& it : this->getStorage()->getEdges() )
      {
         Edge& edge = *it.second;

         if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
         {
            real_t* data = edge.getData( edgeDataID_ )->getPointer( level );
            for ( const auto& iter : vertexdof::macroedge::Iterator( level, 1 ) )
            {
               const uint_t idx = vertexdof::macroedge::indexFromVertex( level, iter.x(), stencilDirection::VERTEX_C );
               data[idx]        = real_c( 1.0 ) / data[idx];
            }
         }
      }

      for ( const auto& it : this->getStorage()->getFaces() )
      {
         Face& face = *it.second;

         if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
         {
            real_t* data = face.getData( faceDataID_ )->getPointer( level );
            for ( const auto& iter : vertexdof::macroface::Iterator( level, 1 ) )
            {
               const uint_t idx =
                   vertexdof::macroface::indexFromVertex( level, iter.col(), iter.row(), stencilDirection::VERTEX_C );
               data[idx] = real_c( 1.0 ) / data[idx];
            }
         }
      }

      for ( const auto& it : this->getStorage()->getCells() )
      {
         Cell& cell = *it.second;

         if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
         {
            real_t* data = cell.getData( cellDataID_ )->getPointer( level );
            for ( const auto& iter : vertexdof::macrocell::Iterator( level, 1 ) )
            {
               const uint_t idx =
                   vertexdof::macrocell::indexFromVertex( level, iter.x(), iter.y(), iter.z(), stencilDirection::VERTEX_C );
               data[idx] = real_c( 1.0 ) / data[idx];
            }
         }
      }
   }

   this->stopTiming( "Invert elementwise" );
}

// ========================
//  explicit instantiation
// ========================
template class VertexDoFFunction< double >;
// template class VertexDoFFunction< float >;
template class VertexDoFFunction< int32_t >;
template class VertexDoFFunction< int64_t >;

template void VertexDoFFunction< double >::interpolateByPrimitiveType< hyteg::Vertex >( const double& constant,
                                                                                        uint_t        level,
                                                                                        DoFType       flag ) const;

template void VertexDoFFunction< double >::interpolateByPrimitiveType< hyteg::Edge >( const double& constant,
                                                                                      uint_t        level,
                                                                                      DoFType       flag ) const;

template void VertexDoFFunction< double >::interpolateByPrimitiveType< hyteg::Face >( const double& constant,
                                                                                      uint_t        level,
                                                                                      DoFType       flag ) const;

template void VertexDoFFunction< double >::interpolateByPrimitiveType< hyteg::Cell >( const double& constant,
                                                                                      uint_t        level,
                                                                                      DoFType       flag ) const;

} // namespace vertexdof
} // namespace hyteg
