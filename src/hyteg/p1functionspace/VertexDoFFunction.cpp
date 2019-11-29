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

#include <utility>

#include "hyteg/Function.hpp"
#include "hyteg/FunctionMemory.hpp"
#include "hyteg/FunctionProperties.hpp"
#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/geometry/Intersection.hpp"
#include "hyteg/p1functionspace/VertexDoFAdditivePackInfo.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/p1functionspace/VertexDoFPackInfo.hpp"
#include "hyteg/p1functionspace/generatedKernels/all.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"

namespace hyteg {
namespace vertexdof {

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
: VertexDoFFunction( name, storage, minLevel, maxLevel, BoundaryCondition::create012BC() )
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
   auto cellVertexDoFFunctionMemoryDataHandling = std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Cell > >(
       minLevel, maxLevel, vertexDoFMacroCellFunctionMemorySize );
   auto faceVertexDoFFunctionMemoryDataHandling = std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Face > >(
       minLevel, maxLevel, vertexDoFMacroFaceFunctionMemorySize );
   auto edgeVertexDoFFunctionMemoryDataHandling = std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Edge > >(
       minLevel, maxLevel, vertexDoFMacroEdgeFunctionMemorySize );
   auto vertexVertexDoFFunctionMemoryDataHandling = std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Vertex > >(
       minLevel, maxLevel, vertexDoFMacroVertexFunctionMemorySize );

   storage->addCellData( cellDataID_, cellVertexDoFFunctionMemoryDataHandling, name );
   storage->addFaceData( faceDataID_, faceVertexDoFFunctionMemoryDataHandling, name );
   storage->addEdgeData( edgeDataID_, edgeVertexDoFFunctionMemoryDataHandling, name );
   storage->addVertexData( vertexDataID_, vertexVertexDoFFunctionMemoryDataHandling, name );

   for( uint_t level = minLevel; level <= maxLevel; ++level )
   {
      communicators_[level]->addPackInfo( std::make_shared< VertexDoFPackInfo< ValueType > >(
          level, vertexDataID_, edgeDataID_, faceDataID_, cellDataID_, this->getStorage() ) );
      additiveCommunicators_[level]->addPackInfo(
          std::make_shared< VertexDoFAdditivePackInfo< ValueType > >( level,
                                                                      vertexDataID_,
                                                                      edgeDataID_,
                                                                      faceDataID_,
                                                                      cellDataID_,
                                                                      this->getStorage()) );
   }
}

template < typename ValueType >
BoundaryCondition VertexDoFFunction< ValueType >::getBoundaryCondition() const
{
   return boundaryCondition_;
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::interpolate( const ValueType& constant, uint_t level, DoFType flag ) const
{
   if( isDummy() )
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
void VertexDoFFunction< ValueType >::interpolate( const std::function< ValueType( const Point3D& ) >& expr,
                                                  uint_t                                              level,
                                                  DoFType                                             flag ) const
{
   if( isDummy() )
   {
      return;
   }
   std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) > exprExtended =
       [&expr]( const hyteg::Point3D& x, const std::vector< ValueType >& ) { return expr( x ); };
   interpolateExtended( exprExtended, {}, level, flag );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::interpolate( const std::function< ValueType( const Point3D& ) >& expr,
                                                  uint_t                                              level,
                                                  BoundaryUID                                         boundaryUID ) const
{
   if( isDummy() )
   {
      return;
   }
   std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) > exprExtended =
       [&expr]( const hyteg::Point3D& x, const std::vector< ValueType >& ) { return expr( x ); };
   interpolateExtended( exprExtended, {}, level, boundaryUID );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::interpolateExtended(
    const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
    const std::vector< std::reference_wrapper< const VertexDoFFunction< ValueType > > >& srcFunctions,
    uint_t                                                                               level,
    DoFType                                                                              flag ) const
{
   if( isDummy() )
   {
      return;
   }
   this->startTiming( "Interpolate" );
   // Collect all source IDs in a vector
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > > srcVertexIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >   srcEdgeIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >   srcFaceIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > >   srcCellIDs;

   for( const VertexDoFFunction& function : srcFunctions )
   {
      srcVertexIDs.push_back( function.vertexDataID_ );
      srcEdgeIDs.push_back( function.edgeDataID_ );
      srcFaceIDs.push_back( function.faceDataID_ );
      srcCellIDs.push_back( function.cellDataID_ );
   }

   for( const auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrovertex::interpolate( vertex, vertexDataID_, srcVertexIDs, expr, level );
      }
   }

   for( const auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroedge::interpolate< ValueType >( level, edge, edgeDataID_, srcEdgeIDs, expr );
      }
   }

   for( auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroface::interpolate< ValueType >( level, face, faceDataID_, srcFaceIDs, expr );
      }
   }

   for( const auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrocell::interpolate< ValueType >( level, cell, cellDataID_, srcCellIDs, expr );
      }
   }
   this->stopTiming( "Interpolate" );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::interpolateExtended(
    const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
    const std::vector< std::reference_wrapper< const VertexDoFFunction< ValueType > > >& srcFunctions,
    uint_t                                                                               level,
    BoundaryUID                                                                          boundaryUID ) const
{
   if( isDummy() )
   {
      return;
   }
   this->startTiming( "Interpolate" );
   // Collect all source IDs in a vector
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > > srcVertexIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >   srcEdgeIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >   srcFaceIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > >   srcCellIDs;

   for( const VertexDoFFunction& function : srcFunctions )
   {
      srcVertexIDs.push_back( function.vertexDataID_ );
      srcEdgeIDs.push_back( function.edgeDataID_ );
      srcFaceIDs.push_back( function.faceDataID_ );
      srcCellIDs.push_back( function.cellDataID_ );
   }

   for( const auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      if( boundaryCondition_.getBoundaryUIDFromMeshFlag( vertex.getMeshBoundaryFlag() ) == boundaryUID )
      {
         vertexdof::macrovertex::interpolate( vertex, vertexDataID_, srcVertexIDs, expr, level );
      }
   }

   for( const auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      if( boundaryCondition_.getBoundaryUIDFromMeshFlag( edge.getMeshBoundaryFlag() ) == boundaryUID )
      {
         vertexdof::macroedge::interpolate< ValueType >( level, edge, edgeDataID_, srcEdgeIDs, expr );
      }
   }

   for( auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      if( boundaryCondition_.getBoundaryUIDFromMeshFlag( face.getMeshBoundaryFlag() ) == boundaryUID )
      {
         vertexdof::macroface::interpolate< ValueType >( level, face, faceDataID_, srcFaceIDs, expr );
      }
   }

   for( const auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;

      if( boundaryCondition_.getBoundaryUIDFromMeshFlag( cell.getMeshBoundaryFlag() ) == boundaryUID )
      {
         vertexdof::macrocell::interpolate< ValueType >( level, cell, cellDataID_, srcCellIDs, expr );
      }
   }
   this->stopTiming( "Interpolate" );
}

template < typename ValueType >
real_t VertexDoFFunction< ValueType >::evaluate( const Point3D& coordinates, uint_t level ) const
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
        return vertexdof::macroface::evaluate< ValueType >( level, face, coordinates, faceDataID_ );
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
        WALBERLA_ABORT("Not implemented.")
      }
    }
  }

  WALBERLA_ABORT( "There is no local macro element including a point at the given coordinates " << coordinates )
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::evaluateGradient( const Point3D& coordinates, uint_t level,
                                                       Point3D& gradient ) const
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
        WALBERLA_ABORT("Not implemented.")
      }
    }
  }

  WALBERLA_ABORT( "There is no local macro element including a point at the given coordinates " << coordinates )
}

template< typename ValueType >
void VertexDoFFunction< ValueType >::swap( const VertexDoFFunction< ValueType > & other,
                                           const uint_t & level,
                                           const DoFType & flag ) const
{
   if( isDummy() )
   {
      return;
   }
   this->startTiming( "Swap" );

   for( auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrovertex::swap< ValueType >( level, vertex, other.getVertexDataID(), vertexDataID_ );
      }
   }

   for( auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroedge::swap< ValueType >( level, edge, other.getEdgeDataID(), edgeDataID_ );
      }
   }

   for( auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroface::swap< ValueType >( level, face, other.getFaceDataID(), faceDataID_ );
      }
   }

   for( auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrocell::swap< ValueType >( level, cell, other.getCellDataID(), cellDataID_ );
      }
   }

   this->stopTiming( "Swap" );
}


template < typename ValueType >
void macroFaceAssign( const uint_t & level, Face & face, const std::vector< ValueType > & scalars,
                      const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > > & srcFaceIDs,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Face > & dstFaceID,
                      const PrimitiveStorage & )
{
  vertexdof::macroface::assign< ValueType >( level, face, scalars, srcFaceIDs, dstFaceID ); 
}

template<>
void macroFaceAssign< double >( const uint_t & level, Face & face, const std::vector< double > & scalars,
                                const std::vector< PrimitiveDataID< FunctionMemory< double >, Face > > & srcFaceIDs,
                                const PrimitiveDataID< FunctionMemory< double >, Face > & dstFaceID,
                                const PrimitiveStorage & storage )
{
  if ( hyteg::globalDefines::useGeneratedKernels && scalars.size() == 1 )
  {
     storage.getTimingTree()->start( "1 RHS function" );
     auto dstData = face.getData( dstFaceID )->getPointer( level );
     auto srcData = face.getData( srcFaceIDs.at( 0 ) )->getPointer( level );
     auto scalar  = scalars.at( 0 );
     vertexdof::macroface::generated::assign_2D_macroface_vertexdof_1_rhs_function( dstData, srcData, scalar, static_cast< int32_t >( level ) );
     storage.getTimingTree()->stop( "1 RHS function" );
  }
  else if ( hyteg::globalDefines::useGeneratedKernels && scalars.size() == 2 )
  {
     storage.getTimingTree()->start( "2 RHS functions" );
     auto dstData  = face.getData( dstFaceID )->getPointer( level );
     auto srcData0 = face.getData( srcFaceIDs.at( 0 ) )->getPointer( level );
     auto srcData1 = face.getData( srcFaceIDs.at( 1 ) )->getPointer( level );
     auto scalar0  = scalars.at( 0 );
     auto scalar1  = scalars.at( 1 );
     vertexdof::macroface::generated::assign_2D_macroface_vertexdof_2_rhs_functions( dstData, srcData0, srcData1, scalar0, scalar1, static_cast< int32_t >( level ) );
    storage.getTimingTree()->stop( "2 RHS functions" );
  }
  else if ( hyteg::globalDefines::useGeneratedKernels && scalars.size() == 3 )
  {
     storage.getTimingTree()->start( "3 RHS functions" );
     auto dstData  = face.getData( dstFaceID )->getPointer( level );
     auto srcData0 = face.getData( srcFaceIDs.at( 0 ) )->getPointer( level );
     auto srcData1 = face.getData( srcFaceIDs.at( 1 ) )->getPointer( level );
     auto srcData2 = face.getData( srcFaceIDs.at( 2 ) )->getPointer( level );
     auto scalar0  = scalars.at( 0 );
     auto scalar1  = scalars.at( 1 );
     auto scalar2  = scalars.at( 2 );
     vertexdof::macroface::generated::assign_2D_macroface_vertexdof_3_rhs_functions( dstData, srcData0, srcData1, srcData2, scalar0, scalar1, scalar2, static_cast< int32_t >( level ) );
     storage.getTimingTree()->stop( "3 RHS functions" );
  }
  else
  {
     vertexdof::macroface::assign< double >( level, face, scalars, srcFaceIDs, dstFaceID );
  }
}


template < typename ValueType >
void macroCellAssign( const uint_t & level, Cell & cell, const std::vector< ValueType > & scalars,
                      const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > > & srcCellIDs,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & dstCellID )
{
  vertexdof::macrocell::assign< ValueType >( level, cell, scalars, srcCellIDs, dstCellID );
}

template<>
void macroCellAssign< double >( const uint_t & level, Cell & cell, const std::vector< double > & scalars,
                                const std::vector< PrimitiveDataID< FunctionMemory< double >, Cell > > & srcCellIDs,
                                const PrimitiveDataID< FunctionMemory< double >, Cell > & dstCellID )
{
   if ( hyteg::globalDefines::useGeneratedKernels && scalars.size() == 1 )
   {
      auto dstData = cell.getData( dstCellID )->getPointer( level );
      auto srcData = cell.getData( srcCellIDs.at( 0 ) )->getPointer( level );
      auto scalar  = scalars.at( 0 );
      if ( hyteg::globalDefines::useP1Coloring )
      {
         std::map< uint_t, uint_t > groupFirstIdx;
         groupFirstIdx[0] = vertexdof::macrocell::index( level, 0, 0, 0 );
         groupFirstIdx[1] = vertexdof::macrocell::index( level, 1, 0, 0 );
         groupFirstIdx[2] = vertexdof::macrocell::index( level, 0, 1, 0 );
         groupFirstIdx[3] = vertexdof::macrocell::index( level, 1, 1, 0 );
         groupFirstIdx[4] = vertexdof::macrocell::index( level, 0, 0, 1 );
         groupFirstIdx[5] = vertexdof::macrocell::index( level, 1, 0, 1 );
         groupFirstIdx[6] = vertexdof::macrocell::index( level, 0, 1, 1 );
         groupFirstIdx[7] = vertexdof::macrocell::index( level, 1, 1, 1 );

         vertexdof::macrocell::generated::assign_3D_macrocell_vertexdof_1_rhs_function_colored( &dstData[groupFirstIdx[0]],
                                                                                                &dstData[groupFirstIdx[1]],
                                                                                                &dstData[groupFirstIdx[2]],
                                                                                                &dstData[groupFirstIdx[3]],
                                                                                                &dstData[groupFirstIdx[4]],
                                                                                                &dstData[groupFirstIdx[5]],
                                                                                                &dstData[groupFirstIdx[6]],
                                                                                                &dstData[groupFirstIdx[7]],
                                                                                                &srcData[groupFirstIdx[0]],
                                                                                                &srcData[groupFirstIdx[1]],
                                                                                                &srcData[groupFirstIdx[2]],
                                                                                                &srcData[groupFirstIdx[3]],
                                                                                                &srcData[groupFirstIdx[4]],
                                                                                                &srcData[groupFirstIdx[5]],
                                                                                                &srcData[groupFirstIdx[6]],
                                                                                                &srcData[groupFirstIdx[7]],
                                                                                                scalar,
                                                                                                static_cast< int32_t >( level ) );
      }
      else
      {
         vertexdof::macrocell::generated::assign_3D_macrocell_vertexdof_1_rhs_function(
             dstData, srcData, scalar, static_cast< int32_t >( level ) );
      }
   }
   else if ( hyteg::globalDefines::useGeneratedKernels && scalars.size() == 2 )
   {
      auto dstData  = cell.getData( dstCellID )->getPointer( level );
      auto srcData0 = cell.getData( srcCellIDs.at( 0 ) )->getPointer( level );
      auto srcData1 = cell.getData( srcCellIDs.at( 1 ) )->getPointer( level );
      auto scalar0  = scalars.at( 0 );
      auto scalar1  = scalars.at( 1 );
      if ( hyteg::globalDefines::useP1Coloring )
      {
         std::map< uint_t, uint_t > groupFirstIdx;
         groupFirstIdx[0] = vertexdof::macrocell::index( level, 0, 0, 0 );
         groupFirstIdx[1] = vertexdof::macrocell::index( level, 1, 0, 0 );
         groupFirstIdx[2] = vertexdof::macrocell::index( level, 0, 1, 0 );
         groupFirstIdx[3] = vertexdof::macrocell::index( level, 1, 1, 0 );
         groupFirstIdx[4] = vertexdof::macrocell::index( level, 0, 0, 1 );
         groupFirstIdx[5] = vertexdof::macrocell::index( level, 1, 0, 1 );
         groupFirstIdx[6] = vertexdof::macrocell::index( level, 0, 1, 1 );
         groupFirstIdx[7] = vertexdof::macrocell::index( level, 1, 1, 1 );

         vertexdof::macrocell::generated::assign_3D_macrocell_vertexdof_2_rhs_functions_colored(
             &dstData[groupFirstIdx[0]],
             &dstData[groupFirstIdx[1]],
             &dstData[groupFirstIdx[2]],
             &dstData[groupFirstIdx[3]],
             &dstData[groupFirstIdx[4]],
             &dstData[groupFirstIdx[5]],
             &dstData[groupFirstIdx[6]],
             &dstData[groupFirstIdx[7]],
             &srcData0[groupFirstIdx[0]],
             &srcData0[groupFirstIdx[1]],
             &srcData0[groupFirstIdx[2]],
             &srcData0[groupFirstIdx[3]],
             &srcData0[groupFirstIdx[4]],
             &srcData0[groupFirstIdx[5]],
             &srcData0[groupFirstIdx[6]],
             &srcData0[groupFirstIdx[7]],
             &srcData1[groupFirstIdx[0]],
             &srcData1[groupFirstIdx[1]],
             &srcData1[groupFirstIdx[2]],
             &srcData1[groupFirstIdx[3]],
             &srcData1[groupFirstIdx[4]],
             &srcData1[groupFirstIdx[5]],
             &srcData1[groupFirstIdx[6]],
             &srcData1[groupFirstIdx[7]],
             scalar0,
             scalar1,
             static_cast< int32_t >( level ) );
      }
      else
      {
         vertexdof::macrocell::generated::assign_3D_macrocell_vertexdof_2_rhs_functions(
             dstData, srcData0, srcData1, scalar0, scalar1, static_cast< int32_t >( level ) );
      }
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
      if ( hyteg::globalDefines::useP1Coloring )
      {
         std::map< uint_t, uint_t > groupFirstIdx;
         groupFirstIdx[0] = vertexdof::macrocell::index( level, 0, 0, 0 );
         groupFirstIdx[1] = vertexdof::macrocell::index( level, 1, 0, 0 );
         groupFirstIdx[2] = vertexdof::macrocell::index( level, 0, 1, 0 );
         groupFirstIdx[3] = vertexdof::macrocell::index( level, 1, 1, 0 );
         groupFirstIdx[4] = vertexdof::macrocell::index( level, 0, 0, 1 );
         groupFirstIdx[5] = vertexdof::macrocell::index( level, 1, 0, 1 );
         groupFirstIdx[6] = vertexdof::macrocell::index( level, 0, 1, 1 );
         groupFirstIdx[7] = vertexdof::macrocell::index( level, 1, 1, 1 );

         vertexdof::macrocell::generated::assign_3D_macrocell_vertexdof_3_rhs_functions_colored(
             &dstData[groupFirstIdx[0]],
             &dstData[groupFirstIdx[1]],
             &dstData[groupFirstIdx[2]],
             &dstData[groupFirstIdx[3]],
             &dstData[groupFirstIdx[4]],
             &dstData[groupFirstIdx[5]],
             &dstData[groupFirstIdx[6]],
             &dstData[groupFirstIdx[7]],
             &srcData0[groupFirstIdx[0]],
             &srcData0[groupFirstIdx[1]],
             &srcData0[groupFirstIdx[2]],
             &srcData0[groupFirstIdx[3]],
             &srcData0[groupFirstIdx[4]],
             &srcData0[groupFirstIdx[5]],
             &srcData0[groupFirstIdx[6]],
             &srcData0[groupFirstIdx[7]],
             &srcData1[groupFirstIdx[0]],
             &srcData1[groupFirstIdx[1]],
             &srcData1[groupFirstIdx[2]],
             &srcData1[groupFirstIdx[3]],
             &srcData1[groupFirstIdx[4]],
             &srcData1[groupFirstIdx[5]],
             &srcData1[groupFirstIdx[6]],
             &srcData1[groupFirstIdx[7]],
             &srcData2[groupFirstIdx[0]],
             &srcData2[groupFirstIdx[1]],
             &srcData2[groupFirstIdx[2]],
             &srcData2[groupFirstIdx[3]],
             &srcData2[groupFirstIdx[4]],
             &srcData2[groupFirstIdx[5]],
             &srcData2[groupFirstIdx[6]],
             &srcData2[groupFirstIdx[7]],
             scalar0,
             scalar1,
             scalar2,
             static_cast< int32_t >( level ) );
      }
      else
      {
         vertexdof::macrocell::generated::assign_3D_macrocell_vertexdof_3_rhs_functions(
             dstData, srcData0, srcData1, srcData2, scalar0, scalar1, scalar2, static_cast< int32_t >( level ) );
      }
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
   if( isDummy() )
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

   for( const VertexDoFFunction< ValueType >& function : functions )
   {
      srcVertexIDs.push_back( function.vertexDataID_ );
      srcEdgeIDs.push_back( function.edgeDataID_ );
      srcFaceIDs.push_back( function.faceDataID_ );
      srcCellIDs.push_back( function.cellDataID_ );
   }
   this->getStorage()->getTimingTree()->start( "Vertex" );
   for( const auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrovertex::assign< ValueType >( vertex, scalars, srcVertexIDs, vertexDataID_, level );
      }
   }
   this->getStorage()->getTimingTree()->stop( "Vertex" );
   this->getStorage()->getTimingTree()->start( "Edge" );
   for( const auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroedge::assign< ValueType >( level, edge, scalars, srcEdgeIDs, edgeDataID_ );
      }
   }
   this->getStorage()->getTimingTree()->stop( "Edge" );
   this->getStorage()->getTimingTree()->start( "Face" );
   for( const auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
        macroFaceAssign< ValueType >( level, face, scalars, srcFaceIDs, faceDataID_, *this->getStorage() );
      }
    }

   this->getStorage()->getTimingTree()->stop( "Face" );
   this->getStorage()->getTimingTree()->start( "Cell" );
   for( const auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      if( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
        macroCellAssign< ValueType >( level, cell, scalars, srcCellIDs, cellDataID_ );
      }
   }
   this->getStorage()->getTimingTree()->stop( "Cell" );
   this->stopTiming( "Assign" );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::assign( const P2Function< ValueType >& src,
                                             const uint_t&                  P1Level,
                                             const DoFType&                 flag ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "Assign (P2 -> P1)" );

   const auto P2Level = P1Level - 1;

   WALBERLA_CHECK_GREATER_EQUAL( P2Level, src.getMinLevel() )
   WALBERLA_CHECK_LESS_EQUAL( P2Level, src.getMaxLevel() )

   for ( const auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         auto P1Data = vertex.getData( vertexDataID_ )->getPointer( P1Level );
         auto P2Data = vertex.getData( src.getVertexDoFFunction().getVertexDataID() )->getPointer( P2Level );
         P1Data[0]   = P2Data[0];
      }
   }

   for ( const auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         auto P1Data   = edge.getData( edgeDataID_ )->getPointer( P1Level );
         auto P2Data_v = edge.getData( src.getVertexDoFFunction().getEdgeDataID() )->getPointer( P2Level );
         auto P2Data_e = edge.getData( src.getEdgeDoFFunction().getEdgeDataID() )->getPointer( P2Level );

         for ( auto itIdx : vertexdof::macroedge::Iterator( P2Level ) )
         {
            P1Data[vertexdof::macroedge::index( P1Level, itIdx.x() * 2 )] =
                P2Data_v[vertexdof::macroedge::index( P2Level, itIdx.x() )];
         }
         for ( auto itIdx : hyteg::edgedof::macroedge::Iterator( P2Level ) )
         {
            P1Data[vertexdof::macroedge::index( P1Level, itIdx.x() * 2 + 1 )] =
                P2Data_e[hyteg::edgedof::macroedge::index( P2Level, itIdx.x() )];
         }
      }
   }

   for ( const auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         auto P1Data   = face.getData( faceDataID_ )->getPointer( P1Level );
         auto P2Data_v = face.getData( src.getVertexDoFFunction().getFaceDataID() )->getPointer( P2Level );
         auto P2Data_e = face.getData( src.getEdgeDoFFunction().getFaceDataID() )->getPointer( P2Level );

         for ( auto itIdx : vertexdof::macroface::Iterator( P2Level ) )
         {
            P1Data[vertexdof::macroface::index( P1Level, itIdx.x() * 2, itIdx.y() * 2 )] =
                P2Data_v[vertexdof::macroface::index( P2Level, itIdx.x(), itIdx.y() )];
         }
         for ( auto itIdx : edgedof::macroface::Iterator( P2Level ) )
         {
            P1Data[vertexdof::macroface::index( P1Level, itIdx.x() * 2 + 1, itIdx.y() * 2 )] =
                P2Data_e[edgedof::macroface::index( P2Level, itIdx.x(), itIdx.y(), edgedof::EdgeDoFOrientation::X )];
            P1Data[vertexdof::macroface::index( P1Level, itIdx.x() * 2 + 1, itIdx.y() * 2 + 1 )] =
                P2Data_e[edgedof::macroface::index( P2Level, itIdx.x(), itIdx.y(), edgedof::EdgeDoFOrientation::XY )];
            P1Data[vertexdof::macroface::index( P1Level, itIdx.x() * 2, itIdx.y() * 2 + 1 )] =
                P2Data_e[edgedof::macroface::index( P2Level, itIdx.x(), itIdx.y(), edgedof::EdgeDoFOrientation::Y )];
         }
      }
   }

   for ( const auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         auto P1Data   = cell.getData( cellDataID_ )->getPointer( P1Level );
         auto P2Data_v = cell.getData( src.getVertexDoFFunction().getCellDataID() )->getPointer( P2Level );
         auto P2Data_e = cell.getData( src.getEdgeDoFFunction().getCellDataID() )->getPointer( P2Level );

         for ( auto itIdx : vertexdof::macrocell::Iterator( P2Level ) )
         {
            P1Data[vertexdof::macrocell::index( P1Level, itIdx.x() * 2, itIdx.y() * 2, itIdx.z() * 2 )] =
            P2Data_v[vertexdof::macrocell::index( P2Level, itIdx.x(), itIdx.y(), itIdx.z() )];
         }
         for ( auto itIdx : edgedof::macrocell::Iterator( P2Level ) )
         {
            P1Data[vertexdof::macrocell::index( P1Level, itIdx.x() * 2 + 1, itIdx.y() * 2, itIdx.z() * 2 )] =
            P2Data_e[edgedof::macrocell::index( P2Level, itIdx.x(), itIdx.y(), itIdx.z(), edgedof::EdgeDoFOrientation::X )];

            P1Data[vertexdof::macrocell::index( P1Level, itIdx.x() * 2, itIdx.y() * 2 + 1, itIdx.z() * 2 )] =
            P2Data_e[edgedof::macrocell::index( P2Level, itIdx.x(), itIdx.y(), itIdx.z(), edgedof::EdgeDoFOrientation::Y )];

            P1Data[vertexdof::macrocell::index( P1Level, itIdx.x() * 2, itIdx.y() * 2, itIdx.z() * 2 + 1 )] =
            P2Data_e[edgedof::macrocell::index( P2Level, itIdx.x(), itIdx.y(), itIdx.z(), edgedof::EdgeDoFOrientation::Z )];

            P1Data[vertexdof::macrocell::index( P1Level, itIdx.x() * 2 + 1, itIdx.y() * 2 + 1, itIdx.z() * 2 )] =
            P2Data_e[edgedof::macrocell::index( P2Level, itIdx.x(), itIdx.y(), itIdx.z(), edgedof::EdgeDoFOrientation::XY )];

            P1Data[vertexdof::macrocell::index( P1Level, itIdx.x() * 2 + 1, itIdx.y() * 2, itIdx.z() * 2 + 1)] =
            P2Data_e[edgedof::macrocell::index( P2Level, itIdx.x(), itIdx.y(), itIdx.z(), edgedof::EdgeDoFOrientation::XZ )];

            P1Data[vertexdof::macrocell::index( P1Level, itIdx.x() * 2, itIdx.y() * 2 + 1, itIdx.z() * 2 + 1)] =
            P2Data_e[edgedof::macrocell::index( P2Level, itIdx.x(), itIdx.y(), itIdx.z(), edgedof::EdgeDoFOrientation::YZ )];
         }

         for ( auto itIdx : edgedof::macrocell::IteratorXYZ( P2Level ) )
         {
            P1Data[vertexdof::macrocell::index( P1Level, itIdx.x() * 2 + 1, itIdx.y() * 2 + 1, itIdx.z() * 2 + 1 )] =
            P2Data_e[edgedof::macrocell::index( P2Level, itIdx.x(), itIdx.y(), itIdx.z(), edgedof::EdgeDoFOrientation::XYZ )];
         }
      }
   }

   this->stopTiming( "Assign (P2 -> P1)" );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::add( const ValueType& scalar, const uint_t& level, DoFType flag ) const
{
   if( isDummy() )
   {
      return;
   }
   this->startTiming( "Add" );

   for( const auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrovertex::add< ValueType >( vertex, scalar, vertexDataID_, level );
      }
   }

   for( const auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroedge::add< ValueType >( level, edge, scalar, edgeDataID_ );
      }
   }

   for( const auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroface::add< ValueType >( level, face, scalar, faceDataID_ );
      }
   }

   for( const auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      if( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrocell::add< ValueType >( level, cell, scalar, cellDataID_ );
      }
   }

   this->stopTiming( "Add" );
}


template < typename ValueType >
void macroFaceAdd( const uint_t & level, Face & face, const std::vector< ValueType > & scalars,
                   const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > > & srcFaceIDs,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Face > & dstFaceID,
                   const PrimitiveStorage & )
{
  vertexdof::macroface::add< ValueType >( level, face, scalars, srcFaceIDs, dstFaceID );
}

template<>
void macroFaceAdd< double >( const uint_t & level, Face & face, const std::vector< double > & scalars,
                                const std::vector< PrimitiveDataID< FunctionMemory< double >, Face > > & srcFaceIDs,
                                const PrimitiveDataID< FunctionMemory< double >, Face > & dstFaceID,
                                const PrimitiveStorage & storage )
{
  if ( hyteg::globalDefines::useGeneratedKernels && scalars.size() == 1 )
  {
    storage.getTimingTree()->start( "1 RHS function" );
    auto dstData = face.getData( dstFaceID )->getPointer( level );
    auto srcData = face.getData( srcFaceIDs.at( 0 ) )->getPointer( level );
    auto scalar  = scalars.at( 0 );
    vertexdof::macroface::generated::add_2D_macroface_vertexdof_1_rhs_function( dstData, srcData, scalar, static_cast< int32_t >( level ) );
    storage.getTimingTree()->stop( "1 RHS function" );
  }
  else if ( hyteg::globalDefines::useGeneratedKernels && scalars.size() == 2 )
  {
    storage.getTimingTree()->start( "2 RHS functions" );
    auto dstData  = face.getData( dstFaceID )->getPointer( level );
    auto srcData0 = face.getData( srcFaceIDs.at( 0 ) )->getPointer( level );
    auto srcData1 = face.getData( srcFaceIDs.at( 1 ) )->getPointer( level );
    auto scalar0  = scalars.at( 0 );
    auto scalar1  = scalars.at( 1 );
    vertexdof::macroface::generated::add_2D_macroface_vertexdof_2_rhs_functions( dstData, srcData0, srcData1, scalar0, scalar1, static_cast< int32_t >( level ) );
    storage.getTimingTree()->stop( "2 RHS functions" );
  }
  else if ( hyteg::globalDefines::useGeneratedKernels && scalars.size() == 3 )
  {
    storage.getTimingTree()->start( "3 RHS functions" );
    auto dstData  = face.getData( dstFaceID )->getPointer( level );
    auto srcData0 = face.getData( srcFaceIDs.at( 0 ) )->getPointer( level );
    auto srcData1 = face.getData( srcFaceIDs.at( 1 ) )->getPointer( level );
    auto srcData2 = face.getData( srcFaceIDs.at( 2 ) )->getPointer( level );
    auto scalar0  = scalars.at( 0 );
    auto scalar1  = scalars.at( 1 );
    auto scalar2  = scalars.at( 2 );
    vertexdof::macroface::generated::add_2D_macroface_vertexdof_3_rhs_functions( dstData, srcData0, srcData1, srcData2, scalar0, scalar1, scalar2, static_cast< int32_t >( level ) );
    storage.getTimingTree()->stop( "3 RHS functions" );
  }
  else
  {
    vertexdof::macroface::add< double >( level, face, scalars, srcFaceIDs, dstFaceID );
  }
}


template < typename ValueType >
void macroCellAdd( const uint_t & level, Cell & cell, const std::vector< ValueType > & scalars,
                      const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > > & srcCellIDs,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & dstCellID )
{
  vertexdof::macrocell::add< ValueType >( level, cell, scalars, srcCellIDs, dstCellID );
}

template<>
void macroCellAdd< double >( const uint_t & level, Cell & cell, const std::vector< double > & scalars,
                                const std::vector< PrimitiveDataID< FunctionMemory< double >, Cell > > & srcCellIDs,
                                const PrimitiveDataID< FunctionMemory< double >, Cell > & dstCellID )
{
  if ( hyteg::globalDefines::useGeneratedKernels && scalars.size() == 1 )
  {
    auto dstData = cell.getData( dstCellID )->getPointer( level );
    auto srcData = cell.getData( srcCellIDs.at( 0 ) )->getPointer( level );
    auto scalar  = scalars.at( 0 );
    if ( hyteg::globalDefines::useP1Coloring )
    {
      std::map< uint_t, uint_t > groupFirstIdx;
      groupFirstIdx[0] = vertexdof::macrocell::index( level, 0, 0, 0 );
      groupFirstIdx[1] = vertexdof::macrocell::index( level, 1, 0, 0 );
      groupFirstIdx[2] = vertexdof::macrocell::index( level, 0, 1, 0 );
      groupFirstIdx[3] = vertexdof::macrocell::index( level, 1, 1, 0 );
      groupFirstIdx[4] = vertexdof::macrocell::index( level, 0, 0, 1 );
      groupFirstIdx[5] = vertexdof::macrocell::index( level, 1, 0, 1 );
      groupFirstIdx[6] = vertexdof::macrocell::index( level, 0, 1, 1 );
      groupFirstIdx[7] = vertexdof::macrocell::index( level, 1, 1, 1 );

      vertexdof::macrocell::generated::add_3D_macrocell_vertexdof_1_rhs_function_colored( &dstData[groupFirstIdx[0]],
                                                                                             &dstData[groupFirstIdx[1]],
                                                                                             &dstData[groupFirstIdx[2]],
                                                                                             &dstData[groupFirstIdx[3]],
                                                                                             &dstData[groupFirstIdx[4]],
                                                                                             &dstData[groupFirstIdx[5]],
                                                                                             &dstData[groupFirstIdx[6]],
                                                                                             &dstData[groupFirstIdx[7]],
                                                                                             &srcData[groupFirstIdx[0]],
                                                                                             &srcData[groupFirstIdx[1]],
                                                                                             &srcData[groupFirstIdx[2]],
                                                                                             &srcData[groupFirstIdx[3]],
                                                                                             &srcData[groupFirstIdx[4]],
                                                                                             &srcData[groupFirstIdx[5]],
                                                                                             &srcData[groupFirstIdx[6]],
                                                                                             &srcData[groupFirstIdx[7]],
                                                                                             scalar,
                                                                                             static_cast< int32_t >( level ) );
    }
    else
    {
      vertexdof::macrocell::generated::add_3D_macrocell_vertexdof_1_rhs_function(
      dstData, srcData, scalar, static_cast< int32_t >( level ) );
    }
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
   if( isDummy() )
   {
      return;
   }
   this->startTiming( "Add" );
   // Collect all source IDs in a vector
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > > srcVertexIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >   srcEdgeIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >   srcFaceIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > >   srcCellIDs;

   for( const VertexDoFFunction< ValueType >& function : functions )
   {
      srcVertexIDs.push_back( function.vertexDataID_ );
      srcEdgeIDs.push_back( function.edgeDataID_ );
      srcFaceIDs.push_back( function.faceDataID_ );
      srcCellIDs.push_back( function.cellDataID_ );
   }

   for( const auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrovertex::add( vertex, scalars, srcVertexIDs, vertexDataID_, level );
      }
   }

   for( const auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroedge::add< ValueType >( level, edge, scalars, srcEdgeIDs, edgeDataID_ );
      }
   }

   for( const auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         macroFaceAdd< ValueType >( level, face, scalars, srcFaceIDs, faceDataID_, *this->getStorage() );
      }
   }

   for( const auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      if( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
        macroCellAdd< ValueType >( level, cell, scalars, srcCellIDs, cellDataID_ );
      }
   }
   this->stopTiming( "Add" );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::multElementwise( const std::vector< std::reference_wrapper< const VertexDoFFunction< ValueType > > >& functions,
                                                      const uint_t                                               level,
                                                      const DoFType                                              flag ) const
{
   if( isDummy() )
   {
      return;
   }
   this->startTiming( "Multiply elementwise" );
   // Collect all source IDs in a vector
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > > srcVertexIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >   srcEdgeIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >   srcFaceIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > >   srcCellIDs;

   for( const VertexDoFFunction& function : functions )
   {
      srcVertexIDs.push_back( function.vertexDataID_ );
      srcEdgeIDs.push_back( function.edgeDataID_ );
      srcFaceIDs.push_back( function.faceDataID_ );
      srcCellIDs.push_back( function.cellDataID_ );
   }

   for( const auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrovertex::multElementwise< ValueType >( vertex, srcVertexIDs, vertexDataID_, level );
      }
   }

   for( const auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroedge::multElementwise< ValueType >( level, edge, srcEdgeIDs, edgeDataID_ );
      }
   }

   for( const auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroface::multElementwise< ValueType >( level, face, srcFaceIDs, faceDataID_ );
      }
   }

   for( const auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      if( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrocell::multElementwise< ValueType >( level, cell, srcCellIDs, cellDataID_ );
      }
   }
   this->stopTiming( "Multiply elementwise" );
}

template < typename ValueType >
ValueType VertexDoFFunction< ValueType >::dotGlobal(const VertexDoFFunction< ValueType >& rhs, size_t level, DoFType flag ) const
{
   ValueType scalarProduct = dotLocal( rhs, level, flag );
   this->startTiming( "Dot (reduce)" );
   walberla::mpi::allReduceInplace( scalarProduct, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
   this->stopTiming( "Dot (reduce)" );
   return scalarProduct;
}

template < typename ValueType >
ValueType VertexDoFFunction< ValueType >::dotLocal(const VertexDoFFunction< ValueType >& rhs, size_t level, DoFType flag ) const
{
   if( isDummy() )
   {
      return ValueType( 0 );
   }
   this->startTiming( "Dot (local)" );
   auto scalarProduct = ValueType( 0 );

   for( const auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         scalarProduct += vertexdof::macrovertex::dot( vertex, vertexDataID_, rhs.vertexDataID_, level );
      }
   }

   for( const auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         scalarProduct += vertexdof::macroedge::dot< ValueType >( level, edge, edgeDataID_, rhs.edgeDataID_ );
      }
   }

   for( const auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         scalarProduct += vertexdof::macroface::dot< ValueType >( level, face, faceDataID_, rhs.faceDataID_ );
      }
   }

   for( const auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      if( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         scalarProduct += vertexdof::macrocell::dot< ValueType >( level, cell, cellDataID_, rhs.cellDataID_ );
      }
   }
   this->stopTiming( "Dot (local)" );
   return scalarProduct;
}

template < typename ValueType >
ValueType VertexDoFFunction< ValueType >::sumGlobal( const uint_t & level, const DoFType & flag, const bool & absolute ) const
{
   ValueType sum = sumLocal( level, flag, absolute );
  this->startTiming( "Sum (reduce)" );
  walberla::mpi::allReduceInplace( sum, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
  this->stopTiming( "Sum (reduce)" );
  return sum;
}

template < typename ValueType >
ValueType VertexDoFFunction< ValueType >::sumLocal( const uint_t & level, const DoFType & flag, const bool & absolute ) const
{
   if( isDummy() )
   {
      return ValueType( 0 );
   }
   this->startTiming( "Sum (local)" );
   auto sum = ValueType( 0 );

   for( const auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         sum += vertexdof::macrovertex::sum( level, vertex, vertexDataID_, absolute );
      }
   }

   for( const auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         sum += vertexdof::macroedge::sum< ValueType >( level, edge, edgeDataID_, absolute );
      }
   }

   for( const auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         sum += vertexdof::macroface::sum< ValueType >( level, face, faceDataID_, absolute );
      }
   }

   for( const auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      if( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
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
   if( isDummy() )
   {
      return;
   }

   this->startTiming( "Enumerate" );

   uint_t counter = hyteg::numberOfLocalDoFs< VertexDoFFunctionTag >( *( this->getStorage() ), level );

   std::vector< uint_t > dofs_per_rank = walberla::mpi::allGather( counter );

   auto startOnRank = ValueType( 0 );

   for( uint_t i = 0; i < uint_c( walberla::MPIManager::instance()->rank() ); ++i )
   {
      startOnRank += static_cast< ValueType >( dofs_per_rank[i] );
   }
   enumerate( level, startOnRank );
   this->stopTiming( "Enumerate" );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::enumerate( uint_t level, ValueType& offset ) const
{
   if( isDummy() )
   {
      return;
   }

   for( auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;
      vertexdof::macrovertex::enumerate( level, vertex, vertexDataID_, offset );
   }

   for( auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;
      vertexdof::macroedge::enumerate< ValueType >( level, edge, edgeDataID_, offset );
   }

   for( auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;
      vertexdof::macroface::enumerate< ValueType >( level, face, faceDataID_, offset );
   }

   for( auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      vertexdof::macrocell::enumerate< ValueType >( level, cell, cellDataID_, offset );
   }

   /// in contrast to other methods in the function class enumerate needs to communicate due to its usage in the PETSc solvers
   communication::syncFunctionBetweenPrimitives( *this, level );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::integrateDG( DGFunction< ValueType >&        rhs,
                                                  VertexDoFFunction< ValueType >& rhsP1,
                                                  uint_t                          level,
                                                  DoFType                         flag )
{
   if( isDummy() )
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

   for( auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrovertex::integrateDG< ValueType >(
             vertex, this->getStorage(), rhs.getVertexDataID(), rhsP1.getVertexDataID(), vertexDataID_, level );
      }
   }

   communicators_[level]->template startCommunication< Vertex, Edge >();
   rhsP1.endCommunication< Face, Edge >( level );

   for( auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroedge::integrateDG< ValueType >(
             level, edge, this->getStorage(), rhs.getEdgeDataID(), rhsP1.getEdgeDataID(), edgeDataID_ );
      }
   }

   communicators_[level]->template endCommunication< Vertex, Edge >();
   communicators_[level]->template startCommunication< Edge, Face >();

   for( auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      if( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
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
   if( isDummy() )
   {
      return ValueType( 0 );
   }
   auto localMax = -std::numeric_limits< ValueType >::max();

   for( auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      localMax = std::max( localMax, vertexdof::macrocell::getMaxValue< ValueType >( level, cell, cellDataID_ ) );
   }

   for( auto& it : this->getStorage()->getFaces() )
   {
      Face&         face   = *it.second;
      const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if( testFlag( faceBC, flag ) )
      {
         localMax = std::max( localMax, vertexdof::macroface::getMaxValue< ValueType >( level, face, faceDataID_ ) );
      }
   }

   for( auto& it : this->getStorage()->getEdges() )
   {
      Edge&         edge   = *it.second;
      const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if( testFlag( edgeBC, flag ) )
      {
         localMax = std::max( localMax, vertexdof::macroedge::getMaxValue< ValueType >( level, edge, edgeDataID_ ) );
      }
   }

   for( auto& it : this->getStorage()->getVertices() )
   {
      Vertex&       vertex   = *it.second;
      const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if( testFlag( vertexBC, flag ) )
      {
         localMax = std::max( localMax, vertexdof::macrovertex::getMaxValue< ValueType >( level, vertex, vertexDataID_ ) );
      }
   }

   ValueType globalMax = localMax;
   if( mpiReduce )
   {
      globalMax = walberla::mpi::allReduce( localMax, walberla::mpi::MAX );
   }

   return globalMax;
}

template < typename ValueType >
ValueType VertexDoFFunction< ValueType >::getMaxMagnitude( uint_t level, DoFType flag, bool mpiReduce ) const
{
   if( isDummy() )
   {
      return ValueType( 0 );
   }
   auto localMax = ValueType( 0.0 );

   for( auto& it : this->getStorage()->getCells() )
   {
     Cell& cell = *it.second;
     localMax = std::max( localMax, vertexdof::macrocell::getMaxMagnitude< ValueType >( level, cell, cellDataID_ ) );
   }

   for( auto& it : this->getStorage()->getFaces() )
   {
      Face&         face   = *it.second;
      const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if( testFlag( faceBC, flag ) )
      {
         localMax = std::max( localMax, vertexdof::macroface::getMaxMagnitude< ValueType >( level, face, faceDataID_ ) );
      }
   }

   for( auto& it : this->getStorage()->getEdges() )
   {
      Edge&         edge   = *it.second;
      const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if( testFlag( edgeBC, flag ) )
      {
         localMax = std::max( localMax, vertexdof::macroedge::getMaxMagnitude< ValueType >( level, edge, edgeDataID_ ) );
      }
   }

   for( auto& it : this->getStorage()->getVertices() )
   {
      Vertex&       vertex   = *it.second;
      const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if( testFlag( vertexBC, flag ) )
      {
         localMax = std::max( localMax, vertexdof::macrovertex::getMaxMagnitude< ValueType >( level, vertex, vertexDataID_ ) );
      }
   }

   ValueType globalMax = localMax;
   if( mpiReduce )
   {
      globalMax = walberla::mpi::allReduce( localMax, walberla::mpi::MAX );
   }

   return globalMax;
}

template < typename ValueType >
ValueType VertexDoFFunction< ValueType >::getMinValue( uint_t level, DoFType flag, bool mpiReduce ) const
{
   if( isDummy() )
   {
      return ValueType( 0 );
   }
   auto localMin = std::numeric_limits< ValueType >::max();

   for( auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      localMin = std::min( localMin, vertexdof::macrocell::getMinValue< ValueType >( level, cell, cellDataID_ ) );
   }

   for( auto& it : this->getStorage()->getFaces() )
   {
      Face&         face   = *it.second;
      const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if( testFlag( faceBC, flag ) )
      {
         localMin = std::min( localMin, vertexdof::macroface::getMinValue< ValueType >( level, face, faceDataID_ ) );
      }
   }

   for( auto& it : this->getStorage()->getEdges() )
   {
      Edge&         edge   = *it.second;
      const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if( testFlag( edgeBC, flag ) )
      {
         localMin = std::min( localMin, vertexdof::macroedge::getMinValue< ValueType >( level, edge, edgeDataID_ ) );
      }
   }

   for( auto& it : this->getStorage()->getVertices() )
   {
      Vertex&       vertex   = *it.second;
      const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if( testFlag( vertexBC, flag ) )
      {
         localMin = std::min( localMin, vertexdof::macrovertex::getMinValue< ValueType >( level, vertex, vertexDataID_ ) );
      }
   }

   ValueType globalMin = localMin;
   if( mpiReduce )
   {
      globalMin = -walberla::mpi::allReduce( -localMin, walberla::mpi::MAX );
   }

   return globalMin;
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::setLocalCommunicationMode(
    const communication::BufferedCommunicator::LocalCommunicationMode& localCommunicationMode )
{
   if( isDummy() )
   {
      return;
   }
   for( auto& communicator : communicators_ )
   {
      communicator.second->setLocalCommunicationMode( localCommunicationMode );
   }
   for( auto& communicator : additiveCommunicators_ )
   {
      communicator.second->setLocalCommunicationMode( localCommunicationMode );
   }
}

template < typename ValueType >
template < typename PrimitiveType >
void VertexDoFFunction< ValueType >::interpolateByPrimitiveType( const ValueType& constant, uint_t level, DoFType flag ) const
{
   if( isDummy() )
   {
      return;
   }
   this->startTiming( "Interpolate" );

   if( std::is_same< PrimitiveType, Vertex >::value )
   {
      for( const auto& it : this->getStorage()->getVertices() )
      {
         Vertex& vertex = *it.second;

         if( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
         {
            vertexdof::macrovertex::interpolate( level, vertex, vertexDataID_, constant );
         }
      }
   } else if( std::is_same< PrimitiveType, Edge >::value )
   {
      for( const auto& it : this->getStorage()->getEdges() )
      {
         Edge& edge = *it.second;

         if( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
         {
            vertexdof::macroedge::interpolate( level, edge, edgeDataID_, constant );
         }
      }
   } else if( std::is_same< PrimitiveType, Face >::value )
   {
      for( const auto& it : this->getStorage()->getFaces() )
      {
         Face& face = *it.second;

         if( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
         {
            vertexdof::macroface::interpolate( level, face, faceDataID_, constant );
         }
      }
   } else if( std::is_same< PrimitiveType, Cell >::value )
   {
      for( const auto& it : this->getStorage()->getCells() )
      {
         Cell& cell = *it.second;

         if( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
         {
            vertexdof::macrocell::interpolate( level, cell, cellDataID_, constant );
         }
      }
   }

   this->stopTiming( "Interpolate" );
}

template class VertexDoFFunction< double >;
template class VertexDoFFunction< int >;

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
