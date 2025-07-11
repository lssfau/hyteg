/*
 * Copyright (c) 2017-2023 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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
#include "core/math/KahanSummation.h"

#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/geometry/BlendingHelpers.hpp"
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
: VertexDoFFunction( name, storage, minLevel, maxLevel, boundaryCondition, false )
{}

template < typename ValueType >
VertexDoFFunction< ValueType >::VertexDoFFunction( const std::string&                         name,
                                                   const std::shared_ptr< PrimitiveStorage >& storage,
                                                   uint_t                                     minLevel,
                                                   uint_t                                     maxLevel,
                                                   BoundaryCondition                          boundaryCondition,
                                                   bool                                       addVolumeGhostLayer )
: Function< VertexDoFFunction< ValueType > >( name, storage, minLevel, maxLevel )
, boundaryCondition_( std::move( boundaryCondition ) )
, referenceCounter_( new internal::ReferenceCounter() )
{
   auto cellVertexDoFFunctionMemoryDataHandling = std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Cell > >();
   auto faceVertexDoFFunctionMemoryDataHandling = std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Face > >();
   auto edgeVertexDoFFunctionMemoryDataHandling = std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Edge > >();
   auto vertexVertexDoFFunctionMemoryDataHandling =
       std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Vertex > >();

   storage->addCellData( cellDataID_, cellVertexDoFFunctionMemoryDataHandling, name );
   storage->addFaceData( faceDataID_, faceVertexDoFFunctionMemoryDataHandling, name );
   storage->addEdgeData( edgeDataID_, edgeVertexDoFFunctionMemoryDataHandling, name );
   storage->addVertexData( vertexDataID_, vertexVertexDoFFunctionMemoryDataHandling, name );

   if ( addVolumeGhostLayer )
   {
      WALBERLA_CHECK_GREATER_EQUAL( storage->getAdditionalHaloDepth(),
                                    1,
                                    "You are trying to extend the ghost-layers to the neighbor volume interior for "
                                    "VertexDoFFunctions. This requires an additional halo depth of at least 1 in the "
                                    "PrimitiveStorage. This can be set in the PrimitiveStorage constructor. Bye." )

      if ( !storage->hasGlobalCells() )
      {
         // Create a data handling instance that handles the initialization, serialization, and deserialization of data.
         const auto dofDataHandlingGL = std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Face > >();

         // Create three data IDs for all faces.
         storage->addFaceData( faceGhostLayerDataIDs_[0], dofDataHandlingGL, "VolumeDoFMacroFaceGL0Data" );
         storage->addFaceData( faceGhostLayerDataIDs_[1], dofDataHandlingGL, "VolumeDoFMacroFaceGL1Data" );
         storage->addFaceData( faceGhostLayerDataIDs_[2], dofDataHandlingGL, "VolumeDoFMacroFaceGL2Data" );
      }
      else
      {
         // Create a data handling instance that handles the initialization, serialization, and deserialization of data.
         const auto dofDataHandlingGL = std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Cell > >();

         // Create three data IDs for all faces.
         storage->addCellData( cellGhostLayerDataIDs_[0], dofDataHandlingGL, "VolumeDoFMacroCellGL0Data" );
         storage->addCellData( cellGhostLayerDataIDs_[1], dofDataHandlingGL, "VolumeDoFMacroCellGL1Data" );
         storage->addCellData( cellGhostLayerDataIDs_[2], dofDataHandlingGL, "VolumeDoFMacroCellGL2Data" );
         storage->addCellData( cellGhostLayerDataIDs_[3], dofDataHandlingGL, "VolumeDoFMacroCellGL3Data" );
      }
   }

   for ( uint_t level = minLevel; level <= maxLevel; ++level )
   {
      for ( const auto& it : storage->getVertices() )
      {
         allocateMemory( level, *it.second );
      }
      for ( const auto& it : storage->getEdges() )
      {
         allocateMemory( level, *it.second );
      }
      for ( const auto& it : storage->getFaces() )
      {
         allocateMemory( level, *it.second );
      }
      for ( const auto& it : storage->getCells() )
      {
         allocateMemory( level, *it.second );
      }

      communicators_[level]->addPackInfo( std::make_shared< VertexDoFPackInfo< ValueType > >( level,
                                                                                              vertexDataID_,
                                                                                              edgeDataID_,
                                                                                              faceDataID_,
                                                                                              cellDataID_,
                                                                                              faceGhostLayerDataIDs_,
                                                                                              cellGhostLayerDataIDs_,
                                                                                              this->getStorage() ) );
      additiveCommunicators_[level]->addPackInfo( std::make_shared< VertexDoFAdditivePackInfo< ValueType > >(
          level, vertexDataID_, edgeDataID_, faceDataID_, cellDataID_, this->getStorage() ) );
   }

   referenceCounter_->increaseRefs();
}

template < typename ValueType >
VertexDoFFunction< ValueType >::~VertexDoFFunction()
{
   referenceCounter_->decreaseRefs();
   if ( referenceCounter_->refs() <= 0 )
   {
      // There are no copies of this handle left. We can delete the allocated DoFs.
      deleteFunctionMemory();
   }
}

template < typename ValueType >
VertexDoFFunction< ValueType >::VertexDoFFunction( const VertexDoFFunction< ValueType >& other )
: Function< VertexDoFFunction< ValueType > >( other )
, vertexDataID_( other.vertexDataID_ )
, edgeDataID_( other.edgeDataID_ )
, faceDataID_( other.faceDataID_ )
, cellDataID_( other.cellDataID_ )
, boundaryCondition_( other.boundaryCondition_ )
, referenceCounter_( other.referenceCounter_ )
{
   referenceCounter_->increaseRefs();
}

template < typename ValueType >
VertexDoFFunction< ValueType >& VertexDoFFunction< ValueType >::operator=( const VertexDoFFunction< ValueType >& other )
{
   if ( this == &other )
   {
      return *this;
   }
   else if ( other.referenceCounter_ == referenceCounter_ )
   {
      WALBERLA_CHECK_EQUAL( vertexDataID_, other.vertexDataID_ )
      WALBERLA_CHECK_EQUAL( edgeDataID_, other.edgeDataID_ )
      WALBERLA_CHECK_EQUAL( faceDataID_, other.faceDataID_ )
      WALBERLA_CHECK_EQUAL( cellDataID_, other.cellDataID_ )
      WALBERLA_CHECK_EQUAL( boundaryCondition_, other.boundaryCondition_ )

      WALBERLA_CHECK_EQUAL( this->functionName_, other.functionName_ )
      WALBERLA_CHECK_EQUAL( this->minLevel_, other.minLevel_ )
      WALBERLA_CHECK_EQUAL( this->maxLevel_, other.maxLevel_ )
      WALBERLA_CHECK_EQUAL( this->storage_, other.storage_ )
   }
   else
   {
      if ( this->storage_ == other.storage_ )
      {
         WALBERLA_CHECK_UNEQUAL( vertexDataID_, other.vertexDataID_ )
         WALBERLA_CHECK_UNEQUAL( edgeDataID_, other.edgeDataID_ )
         WALBERLA_CHECK_UNEQUAL( faceDataID_, other.faceDataID_ )
         WALBERLA_CHECK_UNEQUAL( cellDataID_, other.cellDataID_ )
      }

      referenceCounter_->decreaseRefs();

      if ( referenceCounter_->refs() == 0 )
      {
         // There are no copies of this handle left. We can delete the allocated DoFs.
         deleteFunctionMemory();
      }

      vertexDataID_      = other.vertexDataID_;
      edgeDataID_        = other.edgeDataID_;
      faceDataID_        = other.faceDataID_;
      cellDataID_        = other.cellDataID_;
      boundaryCondition_ = other.boundaryCondition_;
      referenceCounter_  = other.referenceCounter_;
      referenceCounter_->increaseRefs();

      this->functionName_ = other.functionName_;
      this->storage_      = other.storage_;
      this->minLevel_     = other.minLevel_;
      this->maxLevel_     = other.maxLevel_;
   }
   return *this;
}

template < typename ValueType >
bool VertexDoFFunction< ValueType >::hasMemoryAllocated( const uint_t& level, const Vertex& vertex ) const
{
   WALBERLA_CHECK( this->getStorage()->vertexExistsLocally( vertex.getID() ) );
   return vertex.hasData( getVertexDataID() ) && vertex.getData( getVertexDataID() )->hasLevel( level );
}

template < typename ValueType >
bool VertexDoFFunction< ValueType >::hasMemoryAllocated( const uint_t& level, const Edge& edge ) const
{
   WALBERLA_CHECK( this->getStorage()->edgeExistsLocally( edge.getID() ) );
   return edge.hasData( getEdgeDataID() ) && edge.getData( getEdgeDataID() )->hasLevel( level );
}

template < typename ValueType >
bool VertexDoFFunction< ValueType >::hasMemoryAllocated( const uint_t& level, const Face& face ) const
{
   WALBERLA_CHECK( this->getStorage()->faceExistsLocally( face.getID() ) );
   return face.hasData( getFaceDataID() ) && face.getData( getFaceDataID() )->hasLevel( level );
}

template < typename ValueType >
bool VertexDoFFunction< ValueType >::hasMemoryAllocated( const uint_t& level, const Cell& cell ) const
{
   WALBERLA_CHECK( this->getStorage()->cellExistsLocally( cell.getID() ) );
   return cell.hasData( getCellDataID() ) && cell.getData( getCellDataID() )->hasLevel( level );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::allocateMemory( const uint_t& level, const Vertex& vertex )
{
   WALBERLA_CHECK( this->getStorage()->vertexExistsLocally( vertex.getID() ) );
   WALBERLA_CHECK( vertex.hasData( getVertexDataID() ) )
   if ( hasMemoryAllocated( level, vertex ) )
      return;
   vertex.getData( getVertexDataID() )->addData( level, vertexDoFMacroVertexFunctionMemorySize( level, vertex ), 0 );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::allocateMemory( const uint_t& level, const Edge& edge )
{
   WALBERLA_CHECK( this->getStorage()->edgeExistsLocally( edge.getID() ) );
   WALBERLA_CHECK( edge.hasData( getEdgeDataID() ) )
   if ( hasMemoryAllocated( level, edge ) )
      return;
   edge.getData( getEdgeDataID() )->addData( level, vertexDoFMacroEdgeFunctionMemorySize( level, edge ), 0 );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::allocateMemory( const uint_t& level, const Face& face )
{
   WALBERLA_CHECK( this->getStorage()->faceExistsLocally( face.getID() ) );
   WALBERLA_CHECK( face.hasData( getFaceDataID() ) )
   if ( hasMemoryAllocated( level, face ) )
      return;
   face.getData( getFaceDataID() )->addData( level, vertexDoFMacroFaceFunctionMemorySize( level, face ), 0 );

   if ( !this->getStorage()->hasGlobalCells() && hasVolumeGhostLayer() )
   {
      for ( uint_t glID = 0; glID < 3; glID++ )
      {
         face.getData( getFaceGLDataID( glID ) )->addData( level, levelinfo::num_microedges_per_edge( level ), 0 );
      }
   }
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::allocateMemory( const uint_t& level, const Cell& cell )
{
   WALBERLA_CHECK( this->getStorage()->cellExistsLocally( cell.getID() ) );
   WALBERLA_CHECK( cell.hasData( getCellDataID() ) )
   if ( hasMemoryAllocated( level, cell ) )
      return;
   cell.getData( getCellDataID() )->addData( level, vertexDoFMacroCellFunctionMemorySize( level, cell ), 0 );

   if ( this->getStorage()->hasGlobalCells() && hasVolumeGhostLayer() )
   {
      for ( uint_t glID = 0; glID < 4; glID++ )
      {
         cell.getData( getCellGLDataID( glID ) )
             ->addData( level, facedof::macroface::numMicroFacesPerMacroFace( level, facedof::FaceType::GRAY ), 0 );
      }
   }
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::deleteMemory( const uint_t& level, const Vertex& vertex )
{
   WALBERLA_CHECK( this->getStorage()->vertexExistsLocally( vertex.getID() ) );
   if ( !hasMemoryAllocated( level, vertex ) )
      return;
   vertex.getData( getVertexDataID() )->deleteData( level );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::deleteMemory( const uint_t& level, const Edge& edge )
{
   WALBERLA_CHECK( this->getStorage()->edgeExistsLocally( edge.getID() ) );
   if ( !hasMemoryAllocated( level, edge ) )
      return;
   edge.getData( getEdgeDataID() )->deleteData( level );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::deleteMemory( const uint_t& level, const Face& face )
{
   WALBERLA_CHECK( this->getStorage()->faceExistsLocally( face.getID() ) );
   if ( !hasMemoryAllocated( level, face ) )
      return;
   face.getData( getFaceDataID() )->deleteData( level );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::deleteMemory( const uint_t& level, const Cell& cell )
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
   interpolate( [constant]( const Point3D& ) { return constant; }, level, boundaryUID );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::interpolate( const std::function< ValueType( const Point3D& ) >& expr,
                                                  uint_t                                              level,
                                                  DoFType                                             flag ) const
{
   std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) > exprExtended =
       [&expr]( const hyteg::Point3D& x, const std::vector< ValueType >& ) { return expr( x ); };
   interpolate( exprExtended, {}, level, flag );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::interpolate( const std::function< ValueType( const Point3D& ) >& expr,
                                                  uint_t                                              level,
                                                  BoundaryUID                                         boundaryUID ) const
{
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
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( vertexIDs.size() ); i++ )
   {
      Vertex& vertex = *this->getStorage()->getVertex( vertexIDs[uint_c( i )] );

      if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrovertex::interpolate( this->getStorage(), vertex, vertexDataID_, srcVertexIDs, expr, level );
      }
   }

   std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
   {
      Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c( i )] );

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroedge::interpolate< ValueType >( this->getStorage(), level, edge, edgeDataID_, srcEdgeIDs, expr );
      }
   }

   std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
   {
      Face& face = *this->getStorage()->getFace( faceIDs[uint_c( i )] );

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroface::interpolate< ValueType >( this->getStorage(), level, face, faceDataID_, srcFaceIDs, expr );
      }
   }

   std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
   {
      Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c( i )] );

      if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrocell::interpolate< ValueType >( this->getStorage(), level, cell, cellDataID_, srcCellIDs, expr );
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
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( vertexIDs.size() ); i++ )
   {
      Vertex& vertex = *this->getStorage()->getVertex( vertexIDs[uint_c( i )] );

      if ( boundaryCondition_.getBoundaryUIDFromMeshFlag( vertex.getMeshBoundaryFlag() ) == boundaryUID )
      {
         vertexdof::macrovertex::interpolate( this->getStorage(), vertex, vertexDataID_, srcVertexIDs, expr, level );
      }
   }

   std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
   {
      Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c( i )] );

      if ( boundaryCondition_.getBoundaryUIDFromMeshFlag( edge.getMeshBoundaryFlag() ) == boundaryUID )
      {
         vertexdof::macroedge::interpolate< ValueType >( this->getStorage(), level, edge, edgeDataID_, srcEdgeIDs, expr );
      }
   }

   std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
   {
      Face& face = *this->getStorage()->getFace( faceIDs[uint_c( i )] );

      if ( boundaryCondition_.getBoundaryUIDFromMeshFlag( face.getMeshBoundaryFlag() ) == boundaryUID )
      {
         vertexdof::macroface::interpolate< ValueType >( this->getStorage(), level, face, faceDataID_, srcFaceIDs, expr );
      }
   }

   std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
   {
      Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c( i )] );

      if ( boundaryCondition_.getBoundaryUIDFromMeshFlag( cell.getMeshBoundaryFlag() ) == boundaryUID )
      {
         vertexdof::macrocell::interpolate< ValueType >( this->getStorage(), level, cell, cellDataID_, srcCellIDs, expr );
      }
   }
   this->stopTiming( "Interpolate" );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::interpolate( const VertexDoFFunction< ValueType >& functionOnParentGrid,
                                                  uint_t                                level,
                                                  uint_t                                srcLevel )
{
   this->startTiming( "Interpolate" );

   auto storage       = this->getStorage();
   auto parentStorage = functionOnParentGrid.getStorage();
   auto threeD        = storage->hasGlobalCells();

   auto getParentID = [&]( const PrimitiveID& id ) -> PrimitiveID {
      // the element also exists on the coarse grid
      if ( functionOnParentGrid.getStorage()->primitiveExistsLocally( id ) )
         return id;

      // the parent element exits on the coarse grid
      if ( functionOnParentGrid.getStorage()->primitiveExistsLocally( id.getParent() ) )
         return id.getParent();

      // the grand-parent element exits on the coarse grid (may happen when element is subject to both red and green refinement)
      if ( functionOnParentGrid.getStorage()->primitiveExistsLocally( id.getParent().getParent() ) )
         return id.getParent().getParent();

      // "pseudo-siblings" from a green refinement step exist on the coarse grid
      for ( auto& localID : functionOnParentGrid.getStorage()->getPrimitiveIDs() )
      {
         if ( localID.getParent() == id.getParent() )
            return PrimitiveID();
      }

      // parent element has "pseudo-siblings" on the coarse grid
      for ( auto& localID : functionOnParentGrid.getStorage()->getPrimitiveIDs() )
      {
         if ( localID.getParent() == id.getParent().getParent() )
            return PrimitiveID();
      }

      WALBERLA_ABORT( "Interpolation between parent and new grid failed: Primitive doesn't exist locally on parentStorage!" );
      return PrimitiveID();
   };

   // update ghostlayers on source
   functionOnParentGrid.communicate< Vertex, Edge >( srcLevel );
   functionOnParentGrid.communicate< Edge, Face >( srcLevel );
   functionOnParentGrid.communicate< Face, Cell >( srcLevel );

   // interpolate function on volume elements, i.e., cells or faces
   auto primitiveIDs = threeD ? storage->getCellIDs() : storage->getFaceIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( primitiveIDs.size() ); i++ )
   {
      const auto& id       = primitiveIDs[uint_c( i )];
      const auto  parentID = getParentID( id );

      // if the element also exists in the coarser mesh, we simply copy the data
      if ( parentID == id && srcLevel == level )
      {
         if ( threeD )
         {
            storage->getCell( id )
                ->getData( cellDataID_ )
                ->copyFrom( *parentStorage->getCell( id )->getData( functionOnParentGrid.getCellDataID() ), level );
         }
         else
         {
            storage->getFace( id )
                ->getData( faceDataID_ )
                ->copyFrom( *parentStorage->getFace( id )->getData( functionOnParentGrid.getFaceDataID() ), level );
         }
      }
      // otherwise, we evaluate the coarse function at each point
      else
      {
         auto expr = [&]( const Point3D& x, const std::vector< ValueType >& ) {
            ValueType value;
            functionOnParentGrid.evaluate( x, srcLevel, value, 0, parentID );
            return value;
         };
         if ( threeD )
            vertexdof::macrocell::interpolate< ValueType >( storage, level, *storage->getCell( id ), cellDataID_, {}, expr, 0 );
         else
            vertexdof::macroface::interpolate< ValueType >( storage, level, *storage->getFace( id ), faceDataID_, {}, expr, 0 );
      }
   }

   // update interface primitives
   this->communicate< Cell, Face >( level );
   this->communicate< Face, Edge >( level );
   this->communicate< Edge, Vertex >( level );

   this->stopTiming( "Interpolate" );
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::setToZero( uint_t level ) const
{
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
bool VertexDoFFunction< ValueType >::evaluate( const Point3D& physicalCoords,
                                               uint_t         level,
                                               ValueType&     value,
                                               real_t         searchToleranceRadius,
                                               PrimitiveID    id,
                                               real_t         distanceTolerance,
                                               bool           useBestGuess ) const
{
   if constexpr ( !std::is_same< ValueType, real_t >::value )
   {
      WALBERLA_UNUSED( physicalCoords );
      WALBERLA_UNUSED( level );
      WALBERLA_UNUSED( value );
      WALBERLA_UNUSED( searchToleranceRadius );
      WALBERLA_UNUSED( id );
      WALBERLA_ABORT( "VertexDoFFunction< ValueType >::evaluate not implemented for requested template parameter" );
      return false;
   }
   else
   {
      auto storage = this->getStorage();
      auto threeD  = storage->hasGlobalCells();

      auto    coordExists = threeD ? storage->cellExistsLocally( id ) : storage->faceExistsLocally( id );
      Point3D computationalCoords;

      if ( coordExists )
      {
         if ( threeD )
            storage->getCell( id )->getGeometryMap()->evalFinv( physicalCoords, computationalCoords );
         else
            storage->getFace( id )->getGeometryMap()->evalFinv( physicalCoords, computationalCoords );
      }
      else
      {
         auto target         = threeD ? mapFromPhysicalToComputationalDomain3D(
                                    storage, physicalCoords, searchToleranceRadius, distanceTolerance, useBestGuess ) :
                                        mapFromPhysicalToComputationalDomain2D(
                                    storage, physicalCoords, searchToleranceRadius, distanceTolerance, useBestGuess );
         coordExists         = std::get< 0 >( target );
         id                  = std::get< 1 >( target );
         computationalCoords = std::get< 2 >( target );
      }

      if ( coordExists )
      {
         value = threeD ? vertexdof::macrocell::evaluate( level, *( storage->getCell( id ) ), computationalCoords, cellDataID_ ) :
                          vertexdof::macroface::evaluate( level, *( storage->getFace( id ) ), computationalCoords, faceDataID_ );
      }

      return coordExists;
   }

   // will not be reached, but some compilers complain otherwise
   return false;
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::evaluateGradient( const Point3D& physicalCoords, uint_t level, Point3D& gradient ) const
{
   if constexpr ( !std::is_same< ValueType, real_t >::value )
   {
      WALBERLA_UNUSED( physicalCoords );
      WALBERLA_UNUSED( level );
      WALBERLA_UNUSED( gradient );
      WALBERLA_ABORT( "P1Function< ValueType >::evaluateGradient not implemented for requested template parameter" );
   }
   else
   {
      // negative value would exclude this alternative feature in finding primitive ID
      real_t searchToleranceRadius = real_c( 1e-12 );

      // Check if 2D or 3D function
      if ( !this->getStorage()->hasGlobalCells() )
      {
         auto [found, faceID, computationalCoords] =
             mapFromPhysicalToComputationalDomain2D( this->getStorage(), physicalCoords, searchToleranceRadius );
         if ( found )
         {
            Face& face = *( this->getStorage()->getFace( faceID ) );

            // evaluate gradient on computational domain
            vertexdof::macroface::evaluateGradient< ValueType >( level, face, computationalCoords, faceDataID_, gradient );

            // transform gradient to physical coordinates
            Matrix2r DFinv;
            face.getGeometryMap()->evalDFinv( physicalCoords, DFinv );
            real_t aux0 = gradient[0];
            real_t aux1 = gradient[1];
            gradient[0] = DFinv( 0, 0 ) * aux0 + DFinv( 0, 1 ) * aux1;
            gradient[1] = DFinv( 1, 0 ) * aux0 + DFinv( 1, 1 ) * aux1;

            return;
         }
      }
      else
      {
         WALBERLA_ABORT( "VertexDoFFunction< real_t >::evaluateGradient() not implemented for 3D" )
      }

      WALBERLA_ABORT( "There is no local macro element including a point at the given mapped back coordinates for "
                      << physicalCoords )
   }
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::swap( const VertexDoFFunction< ValueType >& other,
                                           const uint_t&                         level,
                                           const DoFType&                        flag ) const
{
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
void VertexDoFFunction< ValueType >::copyFrom( const VertexDoFFunction< ValueType >&  other,
                                               const uint_t&                          level,
                                               const std::map< PrimitiveID, uint_t >& localPrimitiveIDsToRank,
                                               const std::map< PrimitiveID, uint_t >& otherPrimitiveIDsToRank ) const
{
   this->startTiming( "Copy" );

   walberla::mpi::BufferSystem        bufferSystem( walberla::mpi::MPIManager::instance()->comm(), 9563 );
   std::set< walberla::mpi::MPIRank > receiverRanks;
   for ( auto it : localPrimitiveIDsToRank )
   {
      receiverRanks.insert( walberla::mpi::MPIRank( it.second ) );
   }
   bufferSystem.setReceiverInfo( receiverRanks, true );

   for ( auto& it : other.getStorage()->getVertices() )
   {
      PrimitiveID otherPrimitiveID = it.first;
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
      PrimitiveID otherPrimitiveID = it.first;
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
      PrimitiveID otherPrimitiveID = it.first;
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
      PrimitiveID otherPrimitiveID = it.first;
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
         PrimitiveID otherID;
         uint_t      primitiveType = 4;
         uint_t      dataSize      = 0;
         ValueType   value;
         ValueType*  dstPointer;

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
                      const PrimitiveStorage&                                                    storage )
{
   if ( hyteg::globalDefines::useGeneratedKernels && scalars.size() == 1 )
   {
      WALBERLA_NON_OPENMP_SECTION()
      {
         storage.getTimingTree()->start( "1 RHS function" );
      }
      auto dstData = face.getData( dstFaceID )->getPointer( level );
      auto srcData = face.getData( srcFaceIDs.at( 0 ) )->getPointer( level );
      auto scalar  = scalars.at( 0 );
      vertexdof::macroface::generated::assign_2D_macroface_vertexdof_1_rhs_function(
          dstData, srcData, scalar, static_cast< int32_t >( level ) );
      WALBERLA_NON_OPENMP_SECTION()
      {
         storage.getTimingTree()->stop( "1 RHS function" );
      }
   }
   else if ( hyteg::globalDefines::useGeneratedKernels && scalars.size() == 2 )
   {
      WALBERLA_NON_OPENMP_SECTION()
      {
         storage.getTimingTree()->start( "2 RHS functions" );
      }
      auto dstData  = face.getData( dstFaceID )->getPointer( level );
      auto srcData0 = face.getData( srcFaceIDs.at( 0 ) )->getPointer( level );
      auto srcData1 = face.getData( srcFaceIDs.at( 1 ) )->getPointer( level );
      auto scalar0  = scalars.at( 0 );
      auto scalar1  = scalars.at( 1 );
      vertexdof::macroface::generated::assign_2D_macroface_vertexdof_2_rhs_functions(
          dstData, srcData0, srcData1, scalar0, scalar1, static_cast< int32_t >( level ) );

      WALBERLA_NON_OPENMP_SECTION()
      {
         storage.getTimingTree()->stop( "2 RHS functions" );
      }
   }
   else if ( hyteg::globalDefines::useGeneratedKernels && scalars.size() == 3 )
   {
      WALBERLA_NON_OPENMP_SECTION()
      {
         storage.getTimingTree()->start( "3 RHS functions" );
      }
      auto dstData  = face.getData( dstFaceID )->getPointer( level );
      auto srcData0 = face.getData( srcFaceIDs.at( 0 ) )->getPointer( level );
      auto srcData1 = face.getData( srcFaceIDs.at( 1 ) )->getPointer( level );
      auto srcData2 = face.getData( srcFaceIDs.at( 2 ) )->getPointer( level );
      auto scalar0  = scalars.at( 0 );
      auto scalar1  = scalars.at( 1 );
      auto scalar2  = scalars.at( 2 );
      vertexdof::macroface::generated::assign_2D_macroface_vertexdof_3_rhs_functions(
          dstData, srcData0, srcData1, srcData2, scalar0, scalar1, scalar2, static_cast< int32_t >( level ) );

      WALBERLA_NON_OPENMP_SECTION()
      {
         storage.getTimingTree()->stop( "3 RHS functions" );
      }
   }
   else
   {
      vertexdof::macroface::assign( level, face, scalars, srcFaceIDs, dstFaceID );
   }
}

template <>
void macroFaceAssign< int32_t >( const uint_t&                                                            level,
                                 Face&                                                                    face,
                                 const std::vector< int32_t >&                                            scalars,
                                 const std::vector< PrimitiveDataID< FunctionMemory< int32_t >, Face > >& srcFaceIDs,
                                 const PrimitiveDataID< FunctionMemory< int32_t >, Face >&                dstFaceID,
                                 const PrimitiveStorage& )
{
   vertexdof::macroface::assign< int32_t >( level, face, scalars, srcFaceIDs, dstFaceID );
}

template <>
void macroFaceAssign< int64_t >( const uint_t&                                                            level,
                                 Face&                                                                    face,
                                 const std::vector< int64_t >&                                            scalars,
                                 const std::vector< PrimitiveDataID< FunctionMemory< int64_t >, Face > >& srcFaceIDs,
                                 const PrimitiveDataID< FunctionMemory< int64_t >, Face >&                dstFaceID,
                                 const PrimitiveStorage& )
{
   vertexdof::macroface::assign< int64_t >( level, face, scalars, srcFaceIDs, dstFaceID );
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
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( vertexIDs.size() ); i++ )
   {
      Vertex& vertex = *this->getStorage()->getVertex( vertexIDs[uint_c( i )] );

      if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrovertex::assign< ValueType >( vertex, scalars, srcVertexIDs, vertexDataID_, level );
      }
   }

   this->getStorage()->getTimingTree()->stop( "Vertex" );
   this->getStorage()->getTimingTree()->start( "Edge" );

   std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
   {
      Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c( i )] );

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroedge::assign< ValueType >( level, edge, scalars, srcEdgeIDs, edgeDataID_ );
      }
   }

   this->getStorage()->getTimingTree()->stop( "Edge" );
   this->getStorage()->getTimingTree()->start( "Face" );

   std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
   {
      Face& face = *this->getStorage()->getFace( faceIDs[uint_c( i )] );

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         macroFaceAssign< ValueType >( level, face, scalars, srcFaceIDs, faceDataID_, *this->getStorage() );
      }
   }

   this->getStorage()->getTimingTree()->stop( "Face" );
   this->getStorage()->getTimingTree()->start( "Cell" );

   std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
   {
      Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c( i )] );
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
   this->startTiming( "Add" );

   std::vector< PrimitiveID > vertexIDs = this->getStorage()->getVertexIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( vertexIDs.size() ); i++ )
   {
      Vertex& vertex = *this->getStorage()->getVertex( vertexIDs[uint_c( i )] );

      if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrovertex::add< ValueType >( vertex, scalar, vertexDataID_, level );
      }
   }

   std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
   {
      Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c( i )] );

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroedge::add< ValueType >( level, edge, scalar, edgeDataID_ );
      }
   }

   std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
   {
      Face& face = *this->getStorage()->getFace( faceIDs[uint_c( i )] );

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroface::add< ValueType >( level, face, scalar, faceDataID_ );
      }
   }

   std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
   {
      Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c( i )] );
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
                   const PrimitiveStorage&                                                    storage )
{
   if ( hyteg::globalDefines::useGeneratedKernels && scalars.size() == 1 )
   {
      WALBERLA_NON_OPENMP_SECTION()
      {
         storage.getTimingTree()->start( "1 RHS function" );
      }
      auto dstData = face.getData( dstFaceID )->getPointer( level );
      auto srcData = face.getData( srcFaceIDs.at( 0 ) )->getPointer( level );
      auto scalar  = scalars.at( 0 );
      vertexdof::macroface::generated::add_2D_macroface_vertexdof_1_rhs_function(
          dstData, srcData, scalar, static_cast< int32_t >( level ) );
      WALBERLA_NON_OPENMP_SECTION()
      {
         storage.getTimingTree()->stop( "1 RHS function" );
      }
   }
   else if ( hyteg::globalDefines::useGeneratedKernels && scalars.size() == 2 )
   {
      WALBERLA_NON_OPENMP_SECTION()
      {
         storage.getTimingTree()->start( "2 RHS functions" );
      }
      auto dstData  = face.getData( dstFaceID )->getPointer( level );
      auto srcData0 = face.getData( srcFaceIDs.at( 0 ) )->getPointer( level );
      auto srcData1 = face.getData( srcFaceIDs.at( 1 ) )->getPointer( level );
      auto scalar0  = scalars.at( 0 );
      auto scalar1  = scalars.at( 1 );
      vertexdof::macroface::generated::add_2D_macroface_vertexdof_2_rhs_functions(
          dstData, srcData0, srcData1, scalar0, scalar1, static_cast< int32_t >( level ) );
      WALBERLA_NON_OPENMP_SECTION()
      {
         storage.getTimingTree()->stop( "2 RHS functions" );
      }
   }
   else if ( hyteg::globalDefines::useGeneratedKernels && scalars.size() == 3 )
   {
      WALBERLA_NON_OPENMP_SECTION()
      {
         storage.getTimingTree()->start( "3 RHS functions" );
      }
      auto dstData  = face.getData( dstFaceID )->getPointer( level );
      auto srcData0 = face.getData( srcFaceIDs.at( 0 ) )->getPointer( level );
      auto srcData1 = face.getData( srcFaceIDs.at( 1 ) )->getPointer( level );
      auto srcData2 = face.getData( srcFaceIDs.at( 2 ) )->getPointer( level );
      auto scalar0  = scalars.at( 0 );
      auto scalar1  = scalars.at( 1 );
      auto scalar2  = scalars.at( 2 );
      vertexdof::macroface::generated::add_2D_macroface_vertexdof_3_rhs_functions(
          dstData, srcData0, srcData1, srcData2, scalar0, scalar1, scalar2, static_cast< int32_t >( level ) );
      WALBERLA_NON_OPENMP_SECTION()
      {
         storage.getTimingTree()->stop( "3 RHS functions" );
      }
   }
   else
   {
      vertexdof::macroface::add( level, face, scalars, srcFaceIDs, dstFaceID );
   }
}

template <>
void macroFaceAdd< int32_t >( const uint_t&                                                            level,
                              Face&                                                                    face,
                              const std::vector< int32_t >&                                            scalars,
                              const std::vector< PrimitiveDataID< FunctionMemory< int32_t >, Face > >& srcFaceIDs,
                              const PrimitiveDataID< FunctionMemory< int32_t >, Face >&                dstFaceID,
                              const PrimitiveStorage& )
{
   vertexdof::macroface::add< int32_t >( level, face, scalars, srcFaceIDs, dstFaceID );
}

template <>
void macroFaceAdd< int64_t >( const uint_t&                                                            level,
                              Face&                                                                    face,
                              const std::vector< int64_t >&                                            scalars,
                              const std::vector< PrimitiveDataID< FunctionMemory< int64_t >, Face > >& srcFaceIDs,
                              const PrimitiveDataID< FunctionMemory< int64_t >, Face >&                dstFaceID,
                              const PrimitiveStorage& )
{
   vertexdof::macroface::add< int64_t >( level, face, scalars, srcFaceIDs, dstFaceID );
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
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( vertexIDs.size() ); i++ )
   {
      Vertex& vertex = *this->getStorage()->getVertex( vertexIDs[uint_c( i )] );

      if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrovertex::add( vertex, scalars, srcVertexIDs, vertexDataID_, level );
      }
   }

   std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
   {
      Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c( i )] );

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroedge::add< ValueType >( level, edge, scalars, srcEdgeIDs, edgeDataID_ );
      }
   }

   std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
   {
      Face& face = *this->getStorage()->getFace( faceIDs[uint_c( i )] );

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         macroFaceAdd< ValueType >( level, face, scalars, srcFaceIDs, faceDataID_, *this->getStorage() );
      }
   }

   std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
   {
      Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c( i )] );
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
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( vertexIDs.size() ); i++ )
   {
      Vertex& vertex = *this->getStorage()->getVertex( vertexIDs[uint_c( i )] );

      if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrovertex::multElementwise< ValueType >( level, vertex, srcVertexIDs, vertexDataID_ );
      }
   }

   std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
   {
      Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c( i )] );

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroedge::multElementwise< ValueType >( level, edge, srcEdgeIDs, edgeDataID_ );
      }
   }

   std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
   {
      Face& face = *this->getStorage()->getFace( faceIDs[uint_c( i )] );

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroface::multElementwise< ValueType >( level, face, srcFaceIDs, faceDataID_ );
      }
   }

   std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
   {
      Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c( i )] );

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
   if constexpr ( !std::is_floating_point< ValueType >::value )
   {
      WALBERLA_UNUSED( level );
      WALBERLA_UNUSED( flag );
      WALBERLA_UNUSED( workOnHalos );
      WALBERLA_ABORT( "VertexDoFFunction< ValueType >::invertElementwise not available for requested ValueType" );
   }
   else
   {
      this->startTiming( "Invert elementwise" );

      if ( workOnHalos )
      {
         for ( const auto& it : this->getStorage()->getVertices() )
         {
            Vertex& vertex = *it.second;

            if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
            {
               ValueType* data = vertex.getData( vertexDataID_ )->getPointer( level );
               uint_t     size = vertex.getData( vertexDataID_ )->getSize( level );
               for ( uint_t k = 0; k < size; ++k )
               {
                  data[k] = walberla::numeric_cast< ValueType >( 1.0 ) / data[k];
                  // data[0]      = walberla::numeric_cast<ValueType>( 1.0 ) / data[0];
               }
            }
         }

         for ( const auto& it : this->getStorage()->getEdges() )
         {
            Edge& edge = *it.second;

            if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
            {
               ValueType* data = edge.getData( edgeDataID_ )->getPointer( level );
               uint_t     size = edge.getData( edgeDataID_ )->getSize( level );
               for ( uint_t k = 0; k < size; ++k )
               {
                  data[k] = walberla::numeric_cast< ValueType >( 1.0 ) / data[k];
               }
            }
         }

         for ( const auto& it : this->getStorage()->getFaces() )
         {
            Face& face = *it.second;

            if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
            {
               ValueType* data = face.getData( faceDataID_ )->getPointer( level );
               uint_t     size = face.getData( faceDataID_ )->getSize( level );
               for ( uint_t k = 0; k < size; ++k )
               {
                  data[k] = walberla::numeric_cast< ValueType >( 1.0 ) / data[k];
               }
            }
         }

         for ( const auto& it : this->getStorage()->getCells() )
         {
            Cell& cell = *it.second;

            if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
            {
               ValueType* data = cell.getData( cellDataID_ )->getPointer( level );
               uint_t     size = cell.getData( cellDataID_ )->getSize( level );
               for ( uint_t k = 0; k < size; ++k )
               {
                  data[k] = walberla::numeric_cast< ValueType >( 1.0 ) / data[k];
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
               ValueType* data = vertex.getData( vertexDataID_ )->getPointer( level );
               data[0]         = walberla::numeric_cast< ValueType >( 1.0 ) / data[0];
            }
         }

         for ( const auto& it : this->getStorage()->getEdges() )
         {
            Edge& edge = *it.second;

            if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
            {
               ValueType* data = edge.getData( edgeDataID_ )->getPointer( level );
               for ( const auto& iter : vertexdof::macroedge::Iterator( level, 1 ) )
               {
                  const uint_t idx = vertexdof::macroedge::indexFromVertex( level, iter.x(), stencilDirection::VERTEX_C );
                  data[idx]        = walberla::numeric_cast< ValueType >( 1.0 ) / data[idx];
               }
            }
         }

         for ( const auto& it : this->getStorage()->getFaces() )
         {
            Face& face = *it.second;

            if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
            {
               ValueType* data = face.getData( faceDataID_ )->getPointer( level );
               for ( const auto& iter : vertexdof::macroface::Iterator( level, 1 ) )
               {
                  const uint_t idx =
                      vertexdof::macroface::indexFromVertex( level, iter.x(), iter.y(), stencilDirection::VERTEX_C );
                  data[idx] = walberla::numeric_cast< ValueType >( 1.0 ) / data[idx];
               }
            }
         }

         for ( const auto& it : this->getStorage()->getCells() )
         {
            Cell& cell = *it.second;

            if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
            {
               ValueType* data = cell.getData( cellDataID_ )->getPointer( level );
               for ( const auto& iter : vertexdof::macrocell::Iterator( level, 1 ) )
               {
                  const uint_t idx =
                      vertexdof::macrocell::indexFromVertex( level, iter.x(), iter.y(), iter.z(), stencilDirection::VERTEX_C );
                  data[idx] = walberla::numeric_cast< ValueType >( 1.0 ) / data[idx];
               }
            }
         }
      }

      this->stopTiming( "Invert elementwise" );
   }
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
   this->startTiming( "Dot (local)" );
   auto scalarProduct = ValueType( 0 );

   ValueType                  scalarProductVertices = 0;
   std::vector< PrimitiveID > vertexIDs             = this->getStorage()->getVertexIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for reduction( + : scalarProductVertices )
#endif
   for ( int i = 0; i < int_c( vertexIDs.size() ); i++ )
   {
      Vertex& vertex = *this->getStorage()->getVertex( vertexIDs[uint_c( i )] );

      if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         scalarProductVertices += vertexdof::macrovertex::dot( vertex, vertexDataID_, rhs.vertexDataID_, level );
      }
   }
   scalarProduct += scalarProductVertices;

   if ( level >= 1 )
   {
      ValueType                  scalarProductEdges = 0;
      std::vector< PrimitiveID > edgeIDs            = this->getStorage()->getEdgeIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for reduction( + : scalarProductEdges )
#endif
      for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
      {
         Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c( i )] );

         if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
         {
            scalarProductEdges += vertexdof::macroedge::dot< ValueType >( level, edge, edgeDataID_, rhs.edgeDataID_ );
         }
      }
      scalarProduct += scalarProductEdges;

      ValueType                  scalarProductFaces = 0;
      std::vector< PrimitiveID > faceIDs            = this->getStorage()->getFaceIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for reduction( + : scalarProductFaces )
#endif
      for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
      {
         Face& face = *this->getStorage()->getFace( faceIDs[uint_c( i )] );

         if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
         {
            scalarProductFaces += vertexdof::macroface::dot< ValueType >( level, face, faceDataID_, rhs.faceDataID_ );
         }
      }
      scalarProduct += scalarProductFaces;

      ValueType                  scalarProductCells = 0;
      std::vector< PrimitiveID > cellIDs            = this->getStorage()->getCellIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for reduction( + : scalarProductCells )
#endif
      for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
      {
         Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c( i )] );

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
   this->startTiming( "Sum (local)" );
   walberla::math::KahanAccumulator< ValueType > sum;

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
   return sum.get();
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::enumerate( uint_t level ) const
{
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

   // in contrast to other methods in the function class enumerate needs to communicate due to its usage in the PETSc solvers
   communication::syncFunctionBetweenPrimitives( *this, level );
}

template < typename ValueType >
ValueType VertexDoFFunction< ValueType >::getMaxDoFValue( uint_t level, DoFType flag, bool mpiReduce ) const
{
   auto localMax = -std::numeric_limits< ValueType >::max();

   for ( auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      localMax   = std::max( localMax, vertexdof::macrocell::getMaxDoFValue< ValueType >( level, cell, cellDataID_ ) );
   }

   for ( auto& it : this->getStorage()->getFaces() )
   {
      Face&         face   = *it.second;
      const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         localMax = std::max( localMax, vertexdof::macroface::getMaxDoFValue< ValueType >( level, face, faceDataID_ ) );
      }
   }

   for ( auto& it : this->getStorage()->getEdges() )
   {
      Edge&         edge   = *it.second;
      const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         localMax = std::max( localMax, vertexdof::macroedge::getMaxDoFValue< ValueType >( level, edge, edgeDataID_ ) );
      }
   }

   for ( auto& it : this->getStorage()->getVertices() )
   {
      Vertex&       vertex   = *it.second;
      const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         localMax = std::max( localMax, vertexdof::macrovertex::getMaxDoFValue< ValueType >( level, vertex, vertexDataID_ ) );
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
ValueType VertexDoFFunction< ValueType >::getMaxDoFMagnitude( uint_t level, DoFType flag, bool mpiReduce ) const
{
   auto localMax = ValueType( 0.0 );

   for ( auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      localMax   = std::max( localMax, vertexdof::macrocell::getMaxDoFMagnitude< ValueType >( level, cell, cellDataID_ ) );
   }

   for ( auto& it : this->getStorage()->getFaces() )
   {
      Face&         face   = *it.second;
      const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         localMax = std::max( localMax, vertexdof::macroface::getMaxDoFMagnitude< ValueType >( level, face, faceDataID_ ) );
      }
   }

   for ( auto& it : this->getStorage()->getEdges() )
   {
      Edge&         edge   = *it.second;
      const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         localMax = std::max( localMax, vertexdof::macroedge::getMaxDoFMagnitude< ValueType >( level, edge, edgeDataID_ ) );
      }
   }

   for ( auto& it : this->getStorage()->getVertices() )
   {
      Vertex&       vertex   = *it.second;
      const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         localMax = std::max( localMax, vertexdof::macrovertex::getMaxDoFMagnitude< ValueType >( level, vertex, vertexDataID_ ) );
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
ValueType VertexDoFFunction< ValueType >::getMinDoFValue( uint_t level, DoFType flag, bool mpiReduce ) const
{
   auto localMin = std::numeric_limits< ValueType >::max();

   for ( auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      localMin   = std::min( localMin, vertexdof::macrocell::getMinDoFValue< ValueType >( level, cell, cellDataID_ ) );
   }

   for ( auto& it : this->getStorage()->getFaces() )
   {
      Face&         face   = *it.second;
      const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         localMin = std::min( localMin, vertexdof::macroface::getMinDoFValue< ValueType >( level, face, faceDataID_ ) );
      }
   }

   for ( auto& it : this->getStorage()->getEdges() )
   {
      Edge&         edge   = *it.second;
      const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         localMin = std::min( localMin, vertexdof::macroedge::getMinDoFValue< ValueType >( level, edge, edgeDataID_ ) );
      }
   }

   for ( auto& it : this->getStorage()->getVertices() )
   {
      Vertex&       vertex   = *it.second;
      const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         localMin = std::min( localMin, vertexdof::macrovertex::getMinDoFValue< ValueType >( level, vertex, vertexDataID_ ) );
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
ValueType VertexDoFFunction< ValueType >::reduceGlobal( uint_t                                              level,
                                                        std::function< ValueType( ValueType, ValueType ) >& reduceOperation,
                                                        ValueType                                           initialValue,
                                                        walberla::mpi::Operation                            mpiReduceOperation,
                                                        DoFType                                             flag ) const
{
   initialValue = reduceLocal( level, reduceOperation, initialValue, flag );
   walberla::mpi::allReduceInplace( initialValue, mpiReduceOperation );
   return initialValue;
}

template < typename ValueType >
ValueType VertexDoFFunction< ValueType >::reduceGlobal( uint_t                                              level,
                                                        std::function< ValueType( ValueType, ValueType ) >& reduceOperation,
                                                        ValueType                                           initialValue,
                                                        DoFType                                             flag ) const
{
   ValueType localReduce  = reduceLocal( level, reduceOperation, initialValue, flag );
   auto      gatherVector = walberla::mpi::allGather( localReduce );
   for ( auto localValue : gatherVector )
   {
      localReduce = reduceOperation( initialValue, localValue );
   }
   return localReduce;
}

template < typename ValueType >
ValueType VertexDoFFunction< ValueType >::reduceLocal( uint_t                                              level,
                                                       std::function< ValueType( ValueType, ValueType ) >& reduceOperation,
                                                       ValueType                                           initialValue,
                                                       DoFType                                             flag ) const
{
   for ( auto& it : this->getStorage()->getCells() )
   {
      //Cell&     cell = *it.second;
      auto cellReduced =
          vertexdof::macrocell::reduce< ValueType >( level, reduceOperation, initialValue, *it.second, cellDataID_ );
      initialValue = reduceOperation( initialValue, cellReduced );
   }

   for ( auto& it : this->getStorage()->getFaces() )
   {
      const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( it.second->getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         auto faceReduced =
             vertexdof::macroface::reduce< ValueType >( level, reduceOperation, initialValue, *it.second, faceDataID_ );
         initialValue = reduceOperation( initialValue, faceReduced );
      }
   }

   for ( auto& it : this->getStorage()->getEdges() )
   {
      const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( it.second->getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         auto edgeReduced =
             vertexdof::macroedge::reduce< ValueType >( level, reduceOperation, initialValue, *it.second, edgeDataID_ );
         initialValue = reduceOperation( initialValue, edgeReduced );
      }
   }

   for ( auto& it : this->getStorage()->getVertices() )
   {
      Vertex&       vertex   = *it.second;
      const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         initialValue = reduceOperation( initialValue, vertex.getData( vertexDataID_ )->getPointer( level )[0] );
      }
   }
   return initialValue;
}

template < typename ValueType >
void VertexDoFFunction< ValueType >::setLocalCommunicationMode(
    const communication::BufferedCommunicator::LocalCommunicationMode& localCommunicationMode )
{
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
   this->startTiming( "Interpolate" );

   if ( std::is_same< PrimitiveType, Vertex >::value )
   {
      std::vector< PrimitiveID > vertexIDs = this->getStorage()->getVertexIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
      for ( int i = 0; i < int_c( vertexIDs.size() ); i++ )
      {
         Vertex& vertex = *this->getStorage()->getVertex( vertexIDs[uint_c( i )] );

         if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
         {
            vertexdof::macrovertex::interpolate( level, vertex, vertexDataID_, constant );
         }
      }
   }
   else if ( std::is_same< PrimitiveType, Edge >::value )
   {
      std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
      for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
      {
         Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c( i )] );

         if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
         {
            vertexdof::macroedge::interpolate( level, edge, edgeDataID_, constant );
         }
      }
   }
   else if ( std::is_same< PrimitiveType, Face >::value )
   {
      std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
      for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
      {
         Face& face = *this->getStorage()->getFace( faceIDs[uint_c( i )] );

         if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
         {
            vertexdof::macroface::interpolate( level, face, faceDataID_, constant );
         }
      }
   }
   else if ( std::is_same< PrimitiveType, Cell >::value )
   {
      std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
      for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
      {
         Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c( i )] );

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

// ========================
//  explicit instantiation
// ========================
template class VertexDoFFunction< double >;
template class VertexDoFFunction< float >;
#ifdef WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT
template class VertexDoFFunction< walberla::float16 >;
#endif
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

template void VertexDoFFunction< float >::interpolateByPrimitiveType< hyteg::Vertex >( const float& constant,
                                                                                       uint_t       level,
                                                                                       DoFType      flag ) const;

template void VertexDoFFunction< float >::interpolateByPrimitiveType< hyteg::Edge >( const float& constant,
                                                                                     uint_t       level,
                                                                                     DoFType      flag ) const;

template void VertexDoFFunction< float >::interpolateByPrimitiveType< hyteg::Face >( const float& constant,
                                                                                     uint_t       level,
                                                                                     DoFType      flag ) const;

template void VertexDoFFunction< float >::interpolateByPrimitiveType< hyteg::Cell >( const float& constant,
                                                                                     uint_t       level,
                                                                                     DoFType      flag ) const;

template void VertexDoFFunction< int64_t >::interpolateByPrimitiveType< hyteg::Vertex >( const int64_t& constant,
                                                                                         uint_t         level,
                                                                                         DoFType        flag ) const;

template void VertexDoFFunction< int64_t >::interpolateByPrimitiveType< hyteg::Edge >( const int64_t& constant,
                                                                                       uint_t         level,
                                                                                       DoFType        flag ) const;

template void VertexDoFFunction< int64_t >::interpolateByPrimitiveType< hyteg::Face >( const int64_t& constant,
                                                                                       uint_t         level,
                                                                                       DoFType        flag ) const;

template void VertexDoFFunction< int64_t >::interpolateByPrimitiveType< hyteg::Cell >( const int64_t& constant,
                                                                                       uint_t         level,
                                                                                       DoFType        flag ) const;

} // namespace vertexdof
} // namespace hyteg
