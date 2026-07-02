/*
* Copyright (c) 2017-2026 Nils Kohl, Marcus Mohr.
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

#include "hyteg/mesh/micro/MicroMesh.hpp"

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroCell.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroEdge.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroFace.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/p3functionspace/P3VectorFunction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

namespace hyteg::micromesh {

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Local helper functions for the computation of micro-node positions in the absence of a micro-mesh ///
/////////////////////////////////////////////////////////////////////////////////////////////////////////

static Point3D applyBlending( const Point3D& pos, const Primitive& primitive )
{
   Point3D posBlending;
   primitive.getGeometryMap()->evalF( pos, posBlending );
   return posBlending;
}

//////////////////////
/// Micro-vertices ///
//////////////////////

static Point3D
    microVertexPositionNoMesh( uint_t level, const Vertex& vertex, const indexing::Index& microVertexIndex, bool withBlending )
{
   WALBERLA_UNUSED( level );
   WALBERLA_UNUSED( microVertexIndex );

   Point3D pos = vertex.getCoordinates();

   if ( !withBlending )
   {
      return pos;
   }

   return applyBlending( pos, vertex );
}

static Point3D
    microVertexPositionNoMesh( uint_t level, const Edge& edge, const indexing::Index& microVertexIndex, bool withBlending )
{
   Point3D pos = vertexdof::macroedge::coordinateFromIndex( level, edge, microVertexIndex );

   if ( !withBlending )
   {
      return pos;
   }

   return applyBlending( pos, edge );
}

static Point3D
    microVertexPositionNoMesh( uint_t level, const Face& face, const indexing::Index& microVertexIndex, bool withBlending )
{
   Point3D pos = vertexdof::macroface::coordinateFromIndex( level, face, microVertexIndex );

   if ( !withBlending )
   {
      return pos;
   }

   return applyBlending( pos, face );
}

static Point3D
    microVertexPositionNoMesh( uint_t level, const Cell& cell, const indexing::Index& microVertexIndex, bool withBlending )
{
   Point3D pos = vertexdof::macrocell::coordinateFromIndex( level, cell, microVertexIndex );

   if ( !withBlending )
   {
      return pos;
   }

   return applyBlending( pos, cell );
}

static Point3D microVertexPositionNoMesh( uint_t                                     level,
                                          const std::shared_ptr< PrimitiveStorage >& storage,
                                          PrimitiveID                                primitiveId,
                                          const indexing::Index&                     microVertexIndex,
                                          bool                                       withBlending )
{
   if ( storage->cellExistsLocally( primitiveId ) )
   {
      return microVertexPositionNoMesh( level, *storage->getCell( primitiveId ), microVertexIndex, withBlending );
   }
   else if ( storage->faceExistsLocally( primitiveId ) )
   {
      return microVertexPositionNoMesh( level, *storage->getFace( primitiveId ), microVertexIndex, withBlending );
   }
   else if ( storage->edgeExistsLocally( primitiveId ) )
   {
      return microVertexPositionNoMesh( level, *storage->getEdge( primitiveId ), microVertexIndex, withBlending );
   }
   else if ( storage->vertexExistsLocally( primitiveId ) )
   {
      return microVertexPositionNoMesh( level, *storage->getVertex( primitiveId ), microVertexIndex, withBlending );
   }
   else
   {
      WALBERLA_ABORT( "MicroMesh: Primitive does not exist locally!" );
   }
}

///////////////////
/// Micro-edges ///
///////////////////

static Point3D microEdgeCenterPositionNoMesh( uint_t                      level,
                                              const Edge&                 edge,
                                              const indexing::Index&      microEdgeIndex,
                                              edgedof::EdgeDoFOrientation microEdgeOrientation,
                                              bool                        withBlending )
{
   WALBERLA_UNUSED( microEdgeOrientation );

   Point3D pos = edgedof::macroedge::coordinateFromIndex( level, edge, microEdgeIndex );

   if ( !withBlending )
   {
      return pos;
   }

   return applyBlending( pos, edge );
}

static Point3D microEdgeCenterPositionNoMesh( uint_t                      level,
                                              const Face&                 face,
                                              const indexing::Index&      microEdgeIndex,
                                              edgedof::EdgeDoFOrientation microEdgeOrientation,
                                              bool                        withBlending )
{
   Point3D pos = edgedof::macroface::coordinateFromIndex( level, face, microEdgeIndex, microEdgeOrientation );

   if ( !withBlending )
   {
      return pos;
   }

   return applyBlending( pos, face );
}

static Point3D microEdgeCenterPositionNoMesh( uint_t                      level,
                                              const Cell&                 cell,
                                              const indexing::Index&      microEdgeIndex,
                                              edgedof::EdgeDoFOrientation microEdgeOrientation,
                                              bool                        withBlending )
{
   Point3D pos = edgedof::macrocell::coordinateFromIndex( level, cell, microEdgeIndex, microEdgeOrientation );

   if ( !withBlending )
   {
      return pos;
   }

   return applyBlending( pos, cell );
}

static Point3D microEdgeCenterPositionNoMesh( uint_t                                     level,
                                              const std::shared_ptr< PrimitiveStorage >& storage,
                                              PrimitiveID                                primitiveId,
                                              const indexing::Index&                     microEdgeIndex,
                                              edgedof::EdgeDoFOrientation                microEdgeOrientation,
                                              bool                                       withBlending )
{
   auto cellPtr = storage->getLocalCell( primitiveId );
   if ( cellPtr != nullptr )
   {
      return microEdgeCenterPositionNoMesh( level, *cellPtr, microEdgeIndex, microEdgeOrientation, withBlending );
   }
   auto facePtr = storage->getLocalFace( primitiveId );
   if ( facePtr != nullptr )
   {
      return microEdgeCenterPositionNoMesh( level, *facePtr, microEdgeIndex, microEdgeOrientation, withBlending );
   }
   auto edgePtr = storage->getLocalEdge( primitiveId );
   if ( edgePtr != nullptr )
   {
      return microEdgeCenterPositionNoMesh( level, *edgePtr, microEdgeIndex, microEdgeOrientation, withBlending );
   }

   WALBERLA_ABORT( "MicroMesh: Primitive does not exist locally, or you passed a macro-vertex!" );
}

/////////////////////////////////
/// Micro-mesh initialization ///
/////////////////////////////////

static void initMicroMeshFromMacroMesh( P1VectorFunction< real_t >& p1Mesh, uint_t level, bool withBlending )
{
   auto         storage   = p1Mesh.getStorage();
   const uint_t dimension = p1Mesh.getDimension();

   for ( const auto& [pid, vertex] : storage->getVertices() )
   {
      auto pos = microVertexPositionNoMesh( level, *vertex, indexing::Index( 0, 0, 0 ), withBlending );
      for ( uint_t i = 0; i < dimension; i++ )
      {
         vertex->getData( p1Mesh.component( i ).getVertexDataID() )->getPointer( level )[0] = pos( Eigen::Index( i ) );
      }
   }

   for ( const auto& [pid, edge] : storage->getEdges() )
   {
      for ( auto idx : vertexdof::macroedge::Iterator( level ) )
      {
         auto pos = microVertexPositionNoMesh( level, *edge, idx, withBlending );
         for ( uint_t i = 0; i < dimension; i++ )
         {
            auto data = edge->getData( p1Mesh.component( i ).getEdgeDataID() )->getPointer( level );
            data[vertexdof::macroedge::index( level, idx.x() )] = pos( Eigen::Index( i ) );
         }
      }
   }

   for ( const auto& [pid, face] : storage->getFaces() )
   {
      for ( auto idx : vertexdof::macroface::Iterator( level ) )
      {
         auto pos = microVertexPositionNoMesh( level, *face, idx, withBlending );
         for ( uint_t i = 0; i < dimension; i++ )
         {
            auto data = face->getData( p1Mesh.component( i ).getFaceDataID() )->getPointer( level );
            data[vertexdof::macroface::index( level, idx.x(), idx.y() )] = pos( Eigen::Index( i ) );
         }
      }
   }

   for ( const auto& [pid, cell] : storage->getCells() )
   {
      for ( auto idx : vertexdof::macrocell::Iterator( level ) )
      {
         auto pos = microVertexPositionNoMesh( level, *cell, idx, withBlending );
         for ( uint_t i = 0; i < dimension; i++ )
         {
            auto data = cell->getData( p1Mesh.component( i ).getCellDataID() )->getPointer( level );
            data[vertexdof::macrocell::index( level, idx.x(), idx.y(), idx.z() )] = pos( Eigen::Index( i ) );
         }
      }
   }
}

static void initMicroMeshFromMacroMesh( P2VectorFunction< real_t >& p2Mesh, uint_t level, bool withBlending )
{
   auto         storage   = p2Mesh.getStorage();
   const uint_t dimension = p2Mesh.getDimension();

   for ( const auto& [pid, vertex] : storage->getVertices() )
   {
      auto pos = microVertexPositionNoMesh( level, *vertex, indexing::Index( 0, 0, 0 ), withBlending );
      for ( uint_t i = 0; i < dimension; i++ )
      {
         vertex->getData( p2Mesh.component( i ).getVertexDoFFunction().getVertexDataID() )->getPointer( level )[0] =
             pos( Eigen::Index( i ) );
      }
   }

   for ( const auto& [pid, edge] : storage->getEdges() )
   {
      for ( auto idx : vertexdof::macroedge::Iterator( level ) )
      {
         auto pos = microVertexPositionNoMesh( level, *edge, idx, withBlending );
         for ( uint_t i = 0; i < dimension; i++ )
         {
            auto vdata = edge->getData( p2Mesh.component( i ).getVertexDoFFunction().getEdgeDataID() )->getPointer( level );
            vdata[vertexdof::macroedge::index( level, idx.x() )] = pos( Eigen::Index( i ) );
         }
      }

      for ( auto idx : edgedof::macroedge::Iterator( level ) )
      {
         auto pos = microEdgeCenterPositionNoMesh( level, *edge, idx, edgedof::EdgeDoFOrientation::X, withBlending );
         for ( uint_t i = 0; i < dimension; i++ )
         {
            auto edata = edge->getData( p2Mesh.component( i ).getEdgeDoFFunction().getEdgeDataID() )->getPointer( level );
            edata[edgedof::macroedge::index( level, idx.x() )] = pos( Eigen::Index( i ) );
         }
      }
   }

   for ( const auto& [pid, face] : storage->getFaces() )
   {
      for ( auto idx : vertexdof::macroface::Iterator( level ) )
      {
         auto pos = microVertexPositionNoMesh( level, *face, idx, withBlending );
         for ( uint_t i = 0; i < dimension; i++ )
         {
            auto vdata = face->getData( p2Mesh.component( i ).getVertexDoFFunction().getFaceDataID() )->getPointer( level );
            vdata[vertexdof::macroface::index( level, idx.x(), idx.y() )] = pos( Eigen::Index( i ) );
         }
      }

      for ( const auto& orientation : edgedof::faceLocalEdgeDoFOrientations )
      {
         for ( auto idx : edgedof::macroface::Iterator( level ) )
         {
            auto pos = microEdgeCenterPositionNoMesh( level, *face, idx, orientation, withBlending );
            for ( uint_t i = 0; i < dimension; i++ )
            {
               auto edata = face->getData( p2Mesh.component( i ).getEdgeDoFFunction().getFaceDataID() )->getPointer( level );
               edata[edgedof::macroface::index( level, idx.x(), idx.y(), orientation )] = pos( Eigen::Index( i ) );
            }
         }
      }
   }

   for ( const auto& [pid, cell] : storage->getCells() )
   {
      for ( auto idx : vertexdof::macrocell::Iterator( level ) )
      {
         auto pos = microVertexPositionNoMesh( level, *cell, idx, withBlending );
         for ( uint_t i = 0; i < dimension; i++ )
         {
            auto vdata = cell->getData( p2Mesh.component( i ).getVertexDoFFunction().getCellDataID() )->getPointer( level );
            vdata[vertexdof::macrocell::index( level, idx.x(), idx.y(), idx.z() )] = pos( Eigen::Index( i ) );
         }
      }

      for ( const auto& orientation : edgedof::allEdgeDoFOrientationsWithoutXYZ )
      {
         for ( auto idx : edgedof::macrocell::Iterator( level ) )
         {
            auto pos = microEdgeCenterPositionNoMesh( level, *cell, idx, orientation, withBlending );
            for ( uint_t i = 0; i < dimension; i++ )
            {
               auto edata = cell->getData( p2Mesh.component( i ).getEdgeDoFFunction().getCellDataID() )->getPointer( level );
               edata[edgedof::macrocell::index( level, idx.x(), idx.y(), idx.z(), orientation )] = pos( Eigen::Index( i ) );
            }
         }
      }

      for ( auto idx : edgedof::macrocell::IteratorXYZ( level ) )
      {
         auto pos = microEdgeCenterPositionNoMesh( level, *cell, idx, edgedof::EdgeDoFOrientation::XYZ, withBlending );
         for ( uint_t i = 0; i < dimension; i++ )
         {
            auto edata = cell->getData( p2Mesh.component( i ).getEdgeDoFFunction().getCellDataID() )->getPointer( level );
            edata[edgedof::macrocell::index( level, idx.x(), idx.y(), idx.z(), edgedof::EdgeDoFOrientation::XYZ )] =
                pos( Eigen::Index( i ) );
         }
      }
   }
}

static void initMicroMeshFromMacroMesh( P3VectorFunction< real_t >& p3Mesh, uint_t level, bool withBlending )
{
   auto         storage   = p3Mesh.getStorage();
   const uint_t dimension = p3Mesh.getDimension();

   for ( const auto& [pid, vertex] : storage->getVertices() )
   {
      auto pos = microVertexPositionNoMesh( level, *vertex, indexing::Index( 0, 0, 0 ), withBlending );
      for ( uint_t i = 0; i < dimension; i++ )
      {
         vertex->getData( p3Mesh.component( i ).getVertexDoFFunction().getVertexDataID() )->getPointer( level )[0] =
             pos( Eigen::Index( i ) );
      }
   }

   for ( const auto& [pid, edge] : storage->getEdges() )
   {
      for ( auto idx : vertexdof::macroedge::Iterator( level ) )
      {
         auto pos = microVertexPositionNoMesh( level, *edge, idx, withBlending );
         for ( uint_t i = 0; i < dimension; i++ )
         {
            auto vdata = edge->getData( p3Mesh.component( i ).getVertexDoFFunction().getEdgeDataID() )->getPointer( level );
            vdata[vertexdof::macroedge::index( level, idx.x() )] = pos( Eigen::Index( i ) );
         }
      }

      for ( auto microEdgeIdx : edgedof::macroedge::Iterator( level ) )
      {
         indexing::Index v0 = microEdgeIdx;
         indexing::Index v1 = microEdgeIdx + indexing::Index( 1, 0, 0 );

         Point3D posBlue = microEdgeArbitraryPosition( storage, pid, level, v0, v1, real_c( 1.0 / 3.0 ) );
         Point3D posRed  = microEdgeArbitraryPosition( storage, pid, level, v0, v1, real_c( 2.0 / 3.0 ) );

         for ( uint_t i = 0; i < dimension; i++ )
         {
            auto eDataBlue = edge->getData( p3Mesh.component( i ).getEdgeDoFFunctionBlue().getEdgeDataID() )->getPointer( level );
            eDataBlue[edgedof::macroedge::index( level, microEdgeIdx.x() )] = posBlue( Eigen::Index( i ) );

            auto eDataRed = edge->getData( p3Mesh.component( i ).getEdgeDoFFunctionRed().getEdgeDataID() )->getPointer( level );
            eDataRed[edgedof::macroedge::index( level, microEdgeIdx.x() )] = posRed( Eigen::Index( i ) );
         }
      }
   }

   for ( const auto& [faceID, face] : storage->getFaces() )
   {
      for ( auto idx : vertexdof::macroface::Iterator( level ) )
      {
         auto pos = microVertexPositionNoMesh( level, *face, idx, withBlending );
         for ( uint_t i = 0; i < dimension; i++ )
         {
            auto vdata = face->getData( p3Mesh.component( i ).getVertexDoFFunction().getFaceDataID() )->getPointer( level );
            vdata[vertexdof::macroface::index( level, idx.x(), idx.y() )] = pos( Eigen::Index( i ) );
         }
      }

      for ( auto microEdgeIdx : edgedof::macroface::Iterator( level ) )
      {
         // horizontal DoFs
         {
            indexing::Index v0 = microEdgeIdx;
            indexing::Index v1 = microEdgeIdx + indexing::Index( 1, 0, 0 );

#ifndef NDEBUG
            std::array< indexing::Index, 2 > offset = edgedof::calcNeighboringVertexDoFIndices( edgedof::EdgeDoFOrientation::X );
            WALBERLA_ASSERT_EQUAL( offset[1], indexing::Index( 1, 0, 0 ) );
#endif

            const Point3D posBlue = microEdgeArbitraryPosition( storage, faceID, level, v0, v1, real_c( 1.0 / 3.0 ) );
            const Point3D posRed  = microEdgeArbitraryPosition( storage, faceID, level, v0, v1, real_c( 2.0 / 3.0 ) );

            for ( uint_t i = 0; i < dimension; i++ )
            {
               auto eDataBlue =
                   face->getData( p3Mesh.component( i ).getEdgeDoFFunctionBlue().getFaceDataID() )->getPointer( level );
               eDataBlue[edgedof::macroface::index( level, microEdgeIdx.x(), microEdgeIdx.y(), edgedof::EdgeDoFOrientation::X )] =
                   posBlue( Eigen::Index( i ) );

               auto eDataRed =
                   face->getData( p3Mesh.component( i ).getEdgeDoFFunctionRed().getFaceDataID() )->getPointer( level );
               eDataRed[edgedof::macroface::index( level, microEdgeIdx.x(), microEdgeIdx.y(), edgedof::EdgeDoFOrientation::X )] =
                   posRed( Eigen::Index( i ) );
            }
         }

         // Vertical DoFs
         {
            indexing::Index v0 = microEdgeIdx;
            indexing::Index v1 = microEdgeIdx + indexing::Index( 0, 1, 0 );

#ifndef NDEBUG
            std::array< indexing::Index, 2 > offset = edgedof::calcNeighboringVertexDoFIndices( edgedof::EdgeDoFOrientation::Y );
            WALBERLA_ASSERT_EQUAL( offset[1], indexing::Index( 0, 1, 0 ) );
#endif

            const Point3D posBlue = microEdgeArbitraryPosition( storage, faceID, level, v0, v1, real_c( 1.0 / 3.0 ) );
            const Point3D posRed  = microEdgeArbitraryPosition( storage, faceID, level, v0, v1, real_c( 2.0 / 3.0 ) );

            for ( uint_t i = 0; i < dimension; i++ )
            {
               auto eDataBlue =
                   face->getData( p3Mesh.component( i ).getEdgeDoFFunctionBlue().getFaceDataID() )->getPointer( level );
               eDataBlue[edgedof::macroface::index( level, microEdgeIdx.x(), microEdgeIdx.y(), edgedof::EdgeDoFOrientation::Y )] =
                   posBlue( Eigen::Index( i ) );

               auto eDataRed =
                   face->getData( p3Mesh.component( i ).getEdgeDoFFunctionRed().getFaceDataID() )->getPointer( level );
               eDataRed[edgedof::macroface::index( level, microEdgeIdx.x(), microEdgeIdx.y(), edgedof::EdgeDoFOrientation::Y )] =
                   posRed( Eigen::Index( i ) );
            }
         }

         // Diagonal DoFs
         {
            indexing::Index v0 = microEdgeIdx + indexing::Index( 1, 0, 0 );
            indexing::Index v1 = microEdgeIdx + indexing::Index( 0, 1, 0 );

#ifndef NDEBUG
            std::array< indexing::Index, 2 > offset = edgedof::calcNeighboringVertexDoFIndices( edgedof::EdgeDoFOrientation::XY );
            WALBERLA_ASSERT_EQUAL( offset[0], indexing::Index( 1, 0, 0 ) );
            WALBERLA_ASSERT_EQUAL( offset[1], indexing::Index( 0, 1, 0 ) );
#endif

            const Point3D posBlue = microEdgeArbitraryPosition( storage, faceID, level, v0, v1, real_c( 1.0 / 3.0 ) );
            const Point3D posRed  = microEdgeArbitraryPosition( storage, faceID, level, v0, v1, real_c( 2.0 / 3.0 ) );

            for ( uint_t i = 0; i < dimension; i++ )
            {
               auto eDataBlue =
                   face->getData( p3Mesh.component( i ).getEdgeDoFFunctionBlue().getFaceDataID() )->getPointer( level );
               eDataBlue[edgedof::macroface::index(
                   level, microEdgeIdx.x(), microEdgeIdx.y(), edgedof::EdgeDoFOrientation::XY )] = posBlue( Eigen::Index( i ) );

               auto eDataRed =
                   face->getData( p3Mesh.component( i ).getEdgeDoFFunctionRed().getFaceDataID() )->getPointer( level );
               eDataRed[edgedof::macroface::index( level, microEdgeIdx.x(), microEdgeIdx.y(), edgedof::EdgeDoFOrientation::XY )] =
                   posRed( Eigen::Index( i ) );
            }
         }
      }

      // FaceDoFs
      for ( auto faceType : facedof::allFaceTypes )
      {
         for ( const auto& idxIt : facedof::macroface::Iterator( level, faceType ) )
         {
            const Point3D pos = micromesh::microFaceCenterPosition( p3Mesh.getStorage(), faceID, level, idxIt, faceType );

            for ( uint_t i = 0; i < dimension; i++ )
            {
               const auto memLayout = p3Mesh.component( i ).getFaceDoFFunction().memoryLayout();
               auto       dofs      = p3Mesh.component( i ).getFaceDoFFunction().dofMemory( faceID, level );

               dofs[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), faceType, 0, 1, level, memLayout )] =
                   pos( Eigen::Index( i ) );
            }
         }
      }
   }

   for ( const auto& [cellID, cell] : storage->getCells() )
   {
      WALBERLA_ABORT( "initMicroMeshFromMacroMesh() not implemented for the p3Mesh and 3D combo, yet!" );
   }
}

//////////////////////////////////////
/// MicroMesh class implementation ///
//////////////////////////////////////

MicroMesh::MicroMesh( const std::shared_ptr< PrimitiveStorage >& storage,
                      uint_t                                     minLevel,
                      uint_t                                     maxLevel,
                      uint_t                                     polynomialDegree,
                      uint_t                                     dimension )
{
   WALBERLA_CHECK_GREATER_EQUAL( dimension, 2, "The only supported mesh dimensions are 2 and 3." )
   WALBERLA_CHECK_LESS_EQUAL( dimension, 3, "The only supported mesh dimensions are 2 and 3." )

   WALBERLA_CHECK_LESS_EQUAL( minLevel, maxLevel, "The MicroMesh's min level should be less or equal to its max level." )

   const bool withBlending = true;

   if ( polynomialDegree == 1 )
   {
      p1_ = std::make_shared< P1VectorFunction< real_t > >( "microMeshP1", storage, minLevel, maxLevel, dimension );
      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         initMicroMeshFromMacroMesh( *p1_, level, withBlending );
      }
   }
   else if ( polynomialDegree == 2 )
   {
      p2_ = std::make_shared< P2VectorFunction< real_t > >( "microMeshP2", storage, minLevel, maxLevel, dimension );
      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         initMicroMeshFromMacroMesh( *p2_, level, withBlending );
      }
   }
   else if ( polynomialDegree == 3 )
   {
      p3_ = std::make_shared< P3VectorFunction< real_t > >( "microMeshP3", storage, minLevel, maxLevel, dimension );
      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         initMicroMeshFromMacroMesh( *p3_, level, withBlending );
      }
   }
   else
   {
      WALBERLA_ABORT( "MicroMesh with polynomial degree " << polynomialDegree << " not supported." );
   }
}

MicroMesh::MicroMesh( const std::shared_ptr< P1VectorFunction< real_t > >& mesh )
: p1_( mesh )
{}

MicroMesh::MicroMesh( const std::shared_ptr< P2VectorFunction< real_t > >& mesh )
: p2_{ mesh }
{}

MicroMesh::MicroMesh( const std::shared_ptr< P3VectorFunction< real_t > >& mesh )
: p3_{ mesh }
{}

uint_t MicroMesh::polynomialDegree() const
{
   if ( p1Mesh() )
   {
      return 1;
   }

   if ( p2Mesh() )
   {
      return 2;
   }

   if ( p3Mesh() )
   {
      return 3;
   }

   WALBERLA_ABORT( "MicroMesh has no function allocated (this should never happen) :/" )
}

uint_t MicroMesh::dimension() const
{
   if ( p1Mesh() )
   {
      return p1Mesh()->getDimension();
   }

   if ( p2Mesh() )
   {
      return p2Mesh()->getDimension();
   }

   if ( p3Mesh() )
   {
      return p3Mesh()->getDimension();
   }

   WALBERLA_ABORT( "MicroMesh has no function allocated (this should never happen) :/" )
}

std::shared_ptr< P1VectorFunction< real_t > > MicroMesh::p1Mesh() const
{
   return p1_;
}

std::shared_ptr< P2VectorFunction< real_t > > MicroMesh::p2Mesh() const
{
   return p2_;
}

std::shared_ptr< P3VectorFunction< real_t > > MicroMesh::p3Mesh() const
{
   return p3_;
}

std::variant< std::shared_ptr< P1VectorFunction< real_t > >,
              std::shared_ptr< P2VectorFunction< real_t > >,
              std::shared_ptr< P3VectorFunction< real_t > > >
    MicroMesh::mesh() const
{
   std::variant< std::shared_ptr< P1VectorFunction< real_t > >,
                 std::shared_ptr< P2VectorFunction< real_t > >,
                 std::shared_ptr< P3VectorFunction< real_t > > >
       mesh;
   if ( p1Mesh() )
   {
      mesh = p1Mesh();
   }
   else if ( p2Mesh() )
   {
      mesh = p2Mesh();
   }
   else if ( p3Mesh() )
   {
      mesh = p3Mesh();
   }
   return mesh;
}

/////////////////////////////////////////////////
/// MicroMesh communication and interpolation ///
/////////////////////////////////////////////////

void communicate( MicroMesh& microMesh, uint_t level )
{
   if ( microMesh.p1Mesh() )
   {
      communication::syncVectorFunctionBetweenPrimitives( *microMesh.p1Mesh(), level );
   }
   else if ( microMesh.p2Mesh() )
   {
      communication::syncVectorFunctionBetweenPrimitives( *microMesh.p2Mesh(), level );
   }
   else if ( microMesh.p3Mesh() )
   {
      communication::syncVectorFunctionBetweenPrimitives( *microMesh.p3Mesh(), level );
   }
}

void communicate( const std::shared_ptr< PrimitiveStorage >& storage, uint_t level )
{
   auto microMesh = storage->getMicroMesh();

   if ( !microMesh )
   {
      return;
   }

   communicate( *microMesh, level );
}

void interpolate( MicroMesh&                                                      microMesh,
                  const std::vector< std::function< real_t( const Point3D& ) > >& map,
                  walberla::uint_t                                                level )
{
   WALBERLA_CHECK_GREATER_EQUAL( map.size(),
                                 microMesh.dimension(),
                                 "MicroMesh: map dimension is smaller than the dimension of the space the mesh is embedded in." )

   if ( microMesh.p1Mesh() )
   {
      microMesh.p1Mesh()->interpolate( map, level );
   }
   else if ( microMesh.p2Mesh() )
   {
      microMesh.p2Mesh()->interpolate( map, level );
   }
   else if ( microMesh.p3Mesh() )
   {
      microMesh.p3Mesh()->interpolate( map, level );
   }
   else
   {
      WALBERLA_ABORT( "MicroMesh has no function allocated :/" )
   }
}

void interpolate( const std::shared_ptr< PrimitiveStorage >&                      storage,
                  const std::vector< std::function< real_t( const Point3D& ) > >& map,
                  walberla::uint_t                                                level )
{
   auto microMesh = storage->getMicroMesh();

   WALBERLA_CHECK_NOT_NULLPTR( microMesh, "MicroMesh: Cannot interpolate, if no mesh was added to the PrimitiveStorage!" )

   interpolate( *microMesh, map, level );
}

void interpolateAndCommunicate( MicroMesh&                                                      microMesh,
                                const std::vector< std::function< real_t( const Point3D& ) > >& map,
                                walberla::uint_t                                                level )
{
   interpolate( microMesh, map, level );
   communicate( microMesh, level );
}

void interpolateAndCommunicate( const std::shared_ptr< PrimitiveStorage >&                      storage,
                                const std::vector< std::function< real_t( const Point3D& ) > >& map,
                                walberla::uint_t                                                level )
{
   interpolate( storage, map, level );
   communicate( storage, level );
}

void interpolateRefinedCoarseMesh( MicroMesh& microMesh, uint_t level, bool withBlending )
{
   if ( microMesh.p1Mesh() )
   {
      initMicroMeshFromMacroMesh( *microMesh.p1Mesh(), level, withBlending );
   }
   else if ( microMesh.p2Mesh() )
   {
      initMicroMeshFromMacroMesh( *microMesh.p2Mesh(), level, withBlending );
   }
   else if ( microMesh.p3Mesh() )
   {
      initMicroMeshFromMacroMesh( *microMesh.p3Mesh(), level, withBlending );
   }
}

void interpolateRefinedCoarseMesh( const std::shared_ptr< PrimitiveStorage >& storage, uint_t level, bool withBlending )
{
   auto microMesh = storage->getMicroMesh();

   WALBERLA_CHECK_NOT_NULLPTR( microMesh, "MicroMesh: Cannot interpolate if no mesh has been added to the PrimitiveStorage!" )

   interpolateRefinedCoarseMesh( *microMesh, level, withBlending );
}

////////////////////////////////////////////////////////////
/// Computation of micro-vertex and micro-edge positions ///
////////////////////////////////////////////////////////////

Point3D microVertexPosition( const std::shared_ptr< PrimitiveStorage >& storage,
                             PrimitiveID                                primitiveId,
                             uint_t                                     level,
                             const indexing::Index&                     microVertexIndex )
{
   WALBERLA_CHECK( storage->primitiveExistsLocally( primitiveId ), "Cannot compute micro-mesh index of non-local primitive." );

   auto microMesh = storage->getMicroMesh();

   if ( !microMesh )
   {
      return microVertexPositionNoMesh( level, storage, primitiveId, microVertexIndex, true );
   }

   WALBERLA_CHECK( microMesh->polynomialDegree() == 1 || microMesh->polynomialDegree() == 2 || microMesh->polynomialDegree() == 3,
                   "Invalid polynomial degree of MicroMesh." )

   const uint_t dim = microMesh->dimension();

   Point3D p;
   p.setZero();

   if ( storage->vertexExistsLocally( primitiveId ) )
   {
      auto vertex = storage->getVertex( primitiveId );
      for ( uint_t i = 0; i < dim; i++ )
      {
         real_t* vdata;
         if ( microMesh->polynomialDegree() == 1 )
         {
            vdata = vertex->getData( microMesh->p1Mesh()->component( i ).getVertexDataID() )->getPointer( level );
         }
         else if ( microMesh->polynomialDegree() == 2 )
         {
            vdata = vertex->getData( microMesh->p2Mesh()->component( i ).getVertexDoFFunction().getVertexDataID() )
                        ->getPointer( level );
         }
         else if ( microMesh->polynomialDegree() == 3 )
         {
            vdata = vertex->getData( microMesh->p3Mesh()->component( i ).getVertexDoFFunction().getVertexDataID() )
                        ->getPointer( level );
         }
         else
         {
            WALBERLA_ABORT( "MicroMesh: polynomial degree not supported." )
         }
         WALBERLA_CHECK_EQUAL( microVertexIndex, indexing::Index( 0, 0, 0 ) );
         p( Eigen::Index( i ) ) = vdata[0];
      }
   }
   else if ( storage->edgeExistsLocally( primitiveId ) )
   {
      auto edge = storage->getEdge( primitiveId );
      for ( uint_t i = 0; i < dim; i++ )
      {
         real_t* vdata;
         if ( microMesh->polynomialDegree() == 1 )
         {
            vdata = edge->getData( microMesh->p1Mesh()->component( i ).getEdgeDataID() )->getPointer( level );
         }
         else if ( microMesh->polynomialDegree() == 2 )
         {
            vdata =
                edge->getData( microMesh->p2Mesh()->component( i ).getVertexDoFFunction().getEdgeDataID() )->getPointer( level );
         }
         else if ( microMesh->polynomialDegree() == 3 )
         {
            vdata =
                edge->getData( microMesh->p3Mesh()->component( i ).getVertexDoFFunction().getEdgeDataID() )->getPointer( level );
         }
         else
         {
            WALBERLA_ABORT( "MicroMesh: polynomial degree not supported." )
         }
         p( Eigen::Index( i ) ) = vdata[vertexdof::macroedge::index( level, microVertexIndex.x() )];
      }
   }
   else if ( storage->faceExistsLocally( primitiveId ) )
   {
      auto face = storage->getFace( primitiveId );
      for ( uint_t i = 0; i < dim; i++ )
      {
         real_t* vdata;
         if ( microMesh->polynomialDegree() == 1 )
         {
            vdata = face->getData( microMesh->p1Mesh()->component( i ).getFaceDataID() )->getPointer( level );
         }
         else if ( microMesh->polynomialDegree() == 2 )
         {
            vdata =
                face->getData( microMesh->p2Mesh()->component( i ).getVertexDoFFunction().getFaceDataID() )->getPointer( level );
         }
         else if ( microMesh->polynomialDegree() == 3 )
         {
            vdata =
                face->getData( microMesh->p3Mesh()->component( i ).getVertexDoFFunction().getFaceDataID() )->getPointer( level );
         }
         else
         {
            WALBERLA_ABORT( "MicroMesh: this should not happen." );
         }

         p( Eigen::Index( i ) ) = vdata[vertexdof::macroface::index( level, microVertexIndex.x(), microVertexIndex.y() )];
      }
   }
   else if ( storage->cellExistsLocally( primitiveId ) )
   {
      auto cell = storage->getCell( primitiveId );
      for ( uint_t i = 0; i < dim; i++ )
      {
         real_t* vdata;
         if ( microMesh->polynomialDegree() == 1 )
         {
            vdata = cell->getData( microMesh->p1Mesh()->component( i ).getCellDataID() )->getPointer( level );
         }
         else if ( microMesh->polynomialDegree() == 2 )
         {
            vdata =
                cell->getData( microMesh->p2Mesh()->component( i ).getVertexDoFFunction().getCellDataID() )->getPointer( level );
         }
         else if ( microMesh->polynomialDegree() == 3 )
         {
            vdata =
                cell->getData( microMesh->p3Mesh()->component( i ).getVertexDoFFunction().getCellDataID() )->getPointer( level );
         }
         else
         {
            WALBERLA_ABORT( "MicroMesh: this should not happen." );
         }

         p( Eigen::Index( i ) ) =
             vdata[vertexdof::macrocell::index( level, microVertexIndex.x(), microVertexIndex.y(), microVertexIndex.z() )];
      }
   }
   else
   {
      WALBERLA_ABORT( "MicroMesh: PrimitiveID not existing locally. " )
   }

   return p;
}

Point3D microEdgeCenterPosition( const std::shared_ptr< PrimitiveStorage >& storage,
                                 PrimitiveID                                primitiveId,
                                 uint_t                                     level,
                                 const indexing::Index&                     microVertexIndexA,
                                 const indexing::Index&                     microVertexIndexB )
{
   return microEdgeCenterPosition( storage,
                                   primitiveId,
                                   level,
                                   edgedof::calcEdgeDoFIndex( microVertexIndexA, microVertexIndexB ),
                                   edgedof::calcEdgeDoFOrientation( microVertexIndexA, microVertexIndexB ) );
}

Point3D microEdgeCenterPosition( const std::shared_ptr< PrimitiveStorage >& storage,
                                 PrimitiveID                                primitiveId,
                                 uint_t                                     level,
                                 const indexing::Index&                     microEdgeIndex,
                                 const edgedof::EdgeDoFOrientation&         microEdgeOrientation )
{
   WALBERLA_CHECK( storage->primitiveExistsLocally( primitiveId ), "Cannot compute micro-mesh index of non-local primitive." );

   auto microMesh = storage->getMicroMesh();

   if ( !microMesh )
   {
      return microEdgeCenterPositionNoMesh( level, storage, primitiveId, microEdgeIndex, microEdgeOrientation, true );
   }

   WALBERLA_CHECK( microMesh->polynomialDegree() == 1 || microMesh->polynomialDegree() == 2 || microMesh->polynomialDegree() == 3,
                   "Invalid polynomial degree of MicroMesh." )

   const uint_t dim = microMesh->dimension();

   Point3D p;
   p.setZero();

   auto microVertexIndexOffsets = edgedof::calcNeighboringVertexDoFIndices( microEdgeOrientation );
   auto microVertexIndexA       = microEdgeIndex + microVertexIndexOffsets[0];
   auto microVertexIndexB       = microEdgeIndex + microVertexIndexOffsets[1];

   if ( storage->edgeExistsLocally( primitiveId ) )
   {
      auto edge = storage->getEdge( primitiveId );
      for ( uint_t i = 0; i < dim; i++ )
      {
         if ( microMesh->polynomialDegree() == 1 )
         {
            real_t* vdata          = edge->getData( microMesh->p1Mesh()->component( i ).getEdgeDataID() )->getPointer( level );
            p( Eigen::Index( i ) ) = real_c( 0.5 ) * ( vdata[vertexdof::macroedge::index( level, microVertexIndexA.x() )] +
                                                       vdata[vertexdof::macroedge::index( level, microVertexIndexB.x() )] );
         }
         else if ( microMesh->polynomialDegree() == 2 )
         {
            real_t* edata =
                edge->getData( microMesh->p2Mesh()->component( i ).getEdgeDoFFunction().getEdgeDataID() )->getPointer( level );
            p( Eigen::Index( i ) ) = edata[edgedof::macroedge::index( level, microEdgeIndex.x() )];
         }
         else if ( microMesh->polynomialDegree() == 3 )
         {
            const P3Function< real_t >& component = microMesh->p3Mesh()->component( i );
            real_t* vData     = edge->getData( component.getVertexDoFFunction().getEdgeDataID() )->getPointer( level );
            real_t* eDataBlue = edge->getData( component.getEdgeDoFFunctionBlue().getEdgeDataID() )->getPointer( level );
            real_t* eDataRed  = edge->getData( component.getEdgeDoFFunctionRed().getEdgeDataID() )->getPointer( level );

            // need to interpolate along edge
            p( Eigen::Index( i ) ) =
                -real_c( 1.0 / 16.0 ) * ( vData[vertexdof::macroedge::index( level, microVertexIndexA.x() )] +
                                          vData[vertexdof::macroedge::index( level, microVertexIndexB.x() )] ) +
                real_c( 9.0 / 16.0 ) * ( eDataBlue[edgedof::macroedge::index( level, microEdgeIndex.x() )] +
                                         eDataRed[edgedof::macroedge::index( level, microEdgeIndex.x() )] );
         }
      }
   }
   else if ( storage->faceExistsLocally( primitiveId ) )
   {
      auto face = storage->getFace( primitiveId );
      for ( uint_t i = 0; i < dim; i++ )
      {
         if ( microMesh->polynomialDegree() == 1 )
         {
            real_t* vdata = face->getData( microMesh->p1Mesh()->component( i ).getFaceDataID() )->getPointer( level );
            p( Eigen::Index( i ) ) =
                real_c( 0.5 ) * ( vdata[vertexdof::macroface::index( level, microVertexIndexA.x(), microVertexIndexA.y() )] +
                                  vdata[vertexdof::macroface::index( level, microVertexIndexB.x(), microVertexIndexB.y() )] );
         }
         else if ( microMesh->polynomialDegree() == 2 )
         {
            real_t* edata =
                face->getData( microMesh->p2Mesh()->component( i ).getEdgeDoFFunction().getFaceDataID() )->getPointer( level );
            p( Eigen::Index( i ) ) =
                edata[edgedof::macroface::index( level, microEdgeIndex.x(), microEdgeIndex.y(), microEdgeOrientation )];
         }
         else if ( microMesh->polynomialDegree() == 3 )
         {
            const P3Function< real_t >& component = microMesh->p3Mesh()->component( i );
            real_t* vData     = face->getData( component.getVertexDoFFunction().getFaceDataID() )->getPointer( level );
            real_t* eDataBlue = face->getData( component.getEdgeDoFFunctionBlue().getFaceDataID() )->getPointer( level );
            real_t* eDataRed  = face->getData( component.getEdgeDoFFunctionRed().getFaceDataID() )->getPointer( level );

            // need to interpolate along edge
            p( Eigen::Index( i ) ) =
                -real_c( 1.0 / 16.0 ) *
                    ( vData[vertexdof::macroface::index( level, microVertexIndexA.x(), microVertexIndexA.y() )] +
                      vData[vertexdof::macroface::index( level, microVertexIndexB.x(), microVertexIndexB.y() )] ) +
                real_c( 9.0 / 16.0 ) * ( eDataBlue[edgedof::macroface::index(
                                             level, microEdgeIndex.x(), microEdgeIndex.y(), microEdgeOrientation )] +
                                         eDataRed[edgedof::macroface::index(
                                             level, microEdgeIndex.x(), microEdgeIndex.y(), microEdgeOrientation )] );
         }
         else
         {
            WALBERLA_ABORT( "MicroMesh: this should not happen." );
         }
      }
   }
   else if ( storage->cellExistsLocally( primitiveId ) )
   {
      auto cell = storage->getCell( primitiveId );
      for ( uint_t i = 0; i < dim; i++ )
      {
         if ( microMesh->polynomialDegree() == 1 )
         {
            real_t* vdata = cell->getData( microMesh->p1Mesh()->component( i ).getCellDataID() )->getPointer( level );
            p( Eigen::Index( i ) ) =
                real_c( 0.5 ) * ( vdata[vertexdof::macrocell::index(
                                      level, microVertexIndexA.x(), microVertexIndexA.y(), microVertexIndexA.z() )] +
                                  vdata[vertexdof::macrocell::index(
                                      level, microVertexIndexB.x(), microVertexIndexB.y(), microVertexIndexB.z() )] );
         }
         else if ( microMesh->polynomialDegree() == 2 )
         {
            real_t* edata =
                cell->getData( microMesh->p2Mesh()->component( i ).getEdgeDoFFunction().getCellDataID() )->getPointer( level );
            p( Eigen::Index( i ) ) = edata[edgedof::macrocell::index(
                level, microEdgeIndex.x(), microEdgeIndex.y(), microEdgeIndex.z(), microEdgeOrientation )];
         }
         else if ( microMesh->polynomialDegree() == 3 )
         {
            const P3Function< real_t >& component = microMesh->p3Mesh()->component( i );
            real_t* vData     = cell->getData( component.getVertexDoFFunction().getCellDataID() )->getPointer( level );
            real_t* eDataBlue = cell->getData( component.getEdgeDoFFunctionBlue().getCellDataID() )->getPointer( level );
            real_t* eDataRed  = cell->getData( component.getEdgeDoFFunctionRed().getCellDataID() )->getPointer( level );

            // need to interpolate along edge
            p( Eigen::Index( i ) ) =
                -real_c( 1.0 / 16.0 ) * ( vData[vertexdof::macrocell::index(
                                              level, microVertexIndexA.x(), microVertexIndexA.y(), microVertexIndexA.z() )] +
                                          vData[vertexdof::macrocell::index(
                                              level, microVertexIndexB.x(), microVertexIndexB.y(), microVertexIndexB.z() )] ) +
                real_c( 9.0 / 16.0 ) *
                    ( eDataBlue[edgedof::macrocell::index(
                          level, microEdgeIndex.x(), microEdgeIndex.y(), microEdgeIndex.z(), microEdgeOrientation )] +
                      eDataRed[edgedof::macrocell::index(
                          level, microEdgeIndex.x(), microEdgeIndex.y(), microEdgeIndex.z(), microEdgeOrientation )] );
         }
         else
         {
            WALBERLA_ABORT( "MicroMesh: this should not happen." );
         }
      }
   }
   else
   {
      WALBERLA_ABORT( "MicroMesh: PrimitiveID not existing locally. " )
   }

   return p;
}

Point3D microEdgeArbitraryPosition( const std::shared_ptr< PrimitiveStorage >& storage,
                                    PrimitiveID                                primitiveId,
                                    uint_t                                     level,
                                    const indexing::Index&                     microVertexIndexA,
                                    const indexing::Index&                     microVertexIndexB,
                                    real_t                                     positionFactor )
{
   WALBERLA_CHECK( storage->primitiveExistsLocally( primitiveId ), "Cannot compute micro-mesh index of non-local primitive." );

   WALBERLA_ASSERT_GREATER_EQUAL( positionFactor, real_c( 0 ) );
   WALBERLA_ASSERT_LESS_EQUAL( positionFactor, real_c( 1 ) );

   auto microMesh = storage->getMicroMesh();

   // for a cubic micromesh we expect that the given vertex indices are consistent in the sense that when
   // the micro edge points from A -> B this is the same direction as the macro edge of the face/cell
   WALBERLA_DEBUG_SECTION()
   {
      if ( microMesh && microMesh->polynomialDegree() == 3 )
      {
         const edgedof::EdgeDoFOrientation microEdgeOrientation =
             edgedof::calcEdgeDoFOrientation( microVertexIndexA, microVertexIndexB );
         const indexing::Index offset = microVertexIndexB - microVertexIndexA;

         // WALBERLA_LOG_INFO( "edge orientation is: " << edgedof::edgeDoFOrientationToString.at( microEdgeOrientation ) );
         // WALBERLA_LOG_INFO( "offsets: x = " << offset.x() << ", y = " << offset.y() );

         if ( ( microEdgeOrientation == edgedof::EdgeDoFOrientation::X && offset.x() == -1 ) ||
              ( microEdgeOrientation == edgedof::EdgeDoFOrientation::Y && offset.y() == -1 ) ||
              ( microEdgeOrientation == edgedof::EdgeDoFOrientation::Z && offset.z() == -1 ) ||
              ( microEdgeOrientation == edgedof::EdgeDoFOrientation::XY && offset.x() == +1 ) ||
              ( microEdgeOrientation == edgedof::EdgeDoFOrientation::XZ && offset.x() == +1 ) ||
              ( microEdgeOrientation == edgedof::EdgeDoFOrientation::YZ && offset.y() == +1 ) ||
              ( microEdgeOrientation == edgedof::EdgeDoFOrientation::XYZ && offset.y() == -1 ) )
         {
            WALBERLA_ABORT( "microEdgeArbitraryPosition() requires a 'consistent' edge orientation when using a cubic micromesh." );
         }
      }
   }

   Point3D position;
   position.setZero();

   if ( storage->edgeExistsLocally( primitiveId ) )
   {
      auto edge = storage->getEdge( primitiveId );
      if ( !microMesh )
      {
         using vertexdof::macroedge::coordinateFromIndex;
         Point3D vertexA = coordinateFromIndex( level, *edge, microVertexIndexA );
         Point3D vertexB = coordinateFromIndex( level, *edge, microVertexIndexB );
         position        = vertexA + positionFactor * ( vertexB - vertexA );
         return applyBlending( position, *edge );
      }
      else
      {
         for ( uint_t i = 0; i < microMesh->dimension(); i++ )
         {
            if ( microMesh->polynomialDegree() == 1 )
            {
               real_t* vdata  = edge->getData( microMesh->p1Mesh()->component( i ).getEdgeDataID() )->getPointer( level );
               real_t  coordA = vdata[vertexdof::macroedge::index( level, microVertexIndexA.x() )];
               real_t  coordB = vdata[vertexdof::macroedge::index( level, microVertexIndexB.x() )];
               position( Eigen::Index( i ) ) = coordA + positionFactor * ( coordB - coordA );
            }
            else if ( microMesh->polynomialDegree() == 2 )
            {
               real_t* vdata = edge->getData( microMesh->p2Mesh()->component( i ).getVertexDoFFunction().getEdgeDataID() )
                                   ->getPointer( level );
               real_t* edata =
                   edge->getData( microMesh->p2Mesh()->component( i ).getEdgeDoFFunction().getEdgeDataID() )->getPointer( level );
               const auto edgeIndex = edgedof::calcEdgeDoFIndex( microVertexIndexA, microVertexIndexB );

               real_t weightA = vdata[vertexdof::macroedge::index( level, microVertexIndexA.x() )];
               real_t weightB = edata[edgedof::macroedge::index( level, edgeIndex.x() )];
               real_t weightC = vdata[vertexdof::macroedge::index( level, microVertexIndexB.x() )];

               const real_t xi = positionFactor;
               position( Eigen::Index( i ) ) =
                   ( real_c( 2 ) * weightA - real_c( 4 ) * weightB + real_c( 2 ) * weightC ) * xi * xi -
                   ( real_c( 3 ) * weightA - real_c( 4 ) * weightB + weightC ) * xi + weightA;
            }
            else if ( microMesh->polynomialDegree() == 3 )
            {
               const auto edgeIndex = edgedof::calcEdgeDoFIndex( microVertexIndexA, microVertexIndexB );

               // maybe not ideal, but a way to check whether we can avoid interpolation
               if ( std::fpclassify( positionFactor - real_c( 1.0 / 3.0 ) ) == FP_ZERO )
               {
                  real_t* blueData = edge->getData( microMesh->p3Mesh()->component( i ).getEdgeDoFFunctionBlue().getEdgeDataID() )
                                         ->getPointer( level );
                  position( Eigen::Index( i ) ) = blueData[edgedof::macroedge::index( level, edgeIndex.x() )];
               }
               else if ( std::fpclassify( positionFactor - real_c( 2.0 / 3.0 ) ) == FP_ZERO )
               {
                  real_t* redData = edge->getData( microMesh->p3Mesh()->component( i ).getEdgeDoFFunctionRed().getEdgeDataID() )
                                        ->getPointer( level );
                  position( Eigen::Index( i ) ) = redData[edgedof::macroedge::index( level, edgeIndex.x() )];
               }
               else
               {
                  WALBERLA_CHECK( microMesh->polynomialDegree() == 3, "Full P3 interpolation not implemented, yet!" );
               }
            }
            else
            {
               WALBERLA_ABORT( "MicroMesh: this should not happen." );
            }
         }
      }
   }
   else if ( storage->faceExistsLocally( primitiveId ) )
   {
      auto face = storage->getFace( primitiveId );
      if ( !microMesh )
      {
         using vertexdof::macroface::coordinateFromIndex;
         Point3D vertexA = coordinateFromIndex( level, *face, microVertexIndexA );
         Point3D vertexB = coordinateFromIndex( level, *face, microVertexIndexB );
         position        = vertexA + positionFactor * ( vertexB - vertexA );
         return applyBlending( position, *face );
      }
      else
      {
         for ( uint_t i = 0; i < microMesh->dimension(); i++ )
         {
            if ( microMesh->polynomialDegree() == 1 )
            {
               real_t* vdata  = face->getData( microMesh->p1Mesh()->component( i ).getFaceDataID() )->getPointer( level );
               real_t  coordA = vdata[vertexdof::macroface::index( level, microVertexIndexA.x(), microVertexIndexA.y() )];
               real_t  coordB = vdata[vertexdof::macroface::index( level, microVertexIndexB.x(), microVertexIndexB.y() )];
               position( Eigen::Index( i ) ) = coordA + positionFactor * ( coordB - coordA );
            }
            else if ( microMesh->polynomialDegree() == 2 )
            {
               real_t* vdata = face->getData( microMesh->p2Mesh()->component( i ).getVertexDoFFunction().getFaceDataID() )
                                   ->getPointer( level );
               real_t* edata =
                   face->getData( microMesh->p2Mesh()->component( i ).getEdgeDoFFunction().getFaceDataID() )->getPointer( level );
               const auto edgeIndex       = edgedof::calcEdgeDoFIndex( microVertexIndexA, microVertexIndexB );
               const auto edgeOrientation = edgedof::calcEdgeDoFOrientation( microVertexIndexA, microVertexIndexB );

               real_t weightA = vdata[vertexdof::macroface::index( level, microVertexIndexA.x(), microVertexIndexA.y() )];
               real_t weightB = edata[edgedof::macroface::index( level, edgeIndex.x(), edgeIndex.y(), edgeOrientation )];
               real_t weightC = vdata[vertexdof::macroface::index( level, microVertexIndexB.x(), microVertexIndexB.y() )];

               const real_t xi = positionFactor;
               position( Eigen::Index( i ) ) =
                   ( real_c( 2 ) * weightA - real_c( 4 ) * weightB + real_c( 2 ) * weightC ) * xi * xi -
                   ( real_c( 3 ) * weightA - real_c( 4 ) * weightB + weightC ) * xi + weightA;
            }
            else if ( microMesh->polynomialDegree() == 3 )
            {
               const auto edgeIndex       = edgedof::calcEdgeDoFIndex( microVertexIndexA, microVertexIndexB );
               const auto edgeOrientation = edgedof::calcEdgeDoFOrientation( microVertexIndexA, microVertexIndexB );

               // maybe not ideal, but a way to check whether we can avoid interpolation
               if ( std::fpclassify( positionFactor - real_c( 1.0 / 3.0 ) ) == FP_ZERO )
               {
                  real_t* blueData = face->getData( microMesh->p3Mesh()->component( i ).getEdgeDoFFunctionBlue().getFaceDataID() )
                                         ->getPointer( level );
                  position( Eigen::Index( i ) ) =
                      blueData[edgedof::macroface::index( level, edgeIndex.x(), edgeIndex.y(), edgeOrientation )];
               }
               else if ( std::fpclassify( positionFactor - real_c( 2.0 / 3.0 ) ) == FP_ZERO )
               {
                  real_t* redData = face->getData( microMesh->p3Mesh()->component( i ).getEdgeDoFFunctionRed().getFaceDataID() )
                                        ->getPointer( level );
                  position( Eigen::Index( i ) ) =
                      redData[edgedof::macroface::index( level, edgeIndex.x(), edgeIndex.y(), edgeOrientation )];
               }
               else
               {
                  WALBERLA_CHECK( microMesh->polynomialDegree() == 3, "Full P3 interpolation not implemented, yet!" );
               }
            }
            else
            {
               WALBERLA_ABORT( "MicroMesh: this should not happen." );
            }
         }
      }
   }
   else if ( storage->cellExistsLocally( primitiveId ) )
   {
      auto cell = storage->getCell( primitiveId );
      if ( !microMesh )
      {
         using vertexdof::macrocell::coordinateFromIndex;
         Point3D vertexA = coordinateFromIndex( level, *cell, microVertexIndexA );
         Point3D vertexB = coordinateFromIndex( level, *cell, microVertexIndexB );
         position        = vertexA + positionFactor * ( vertexB - vertexA );
         return applyBlending( position, *cell );
      }
      else
      {
         for ( uint_t i = 0; i < microMesh->dimension(); i++ )
         {
            if ( microMesh->polynomialDegree() == 1 )
            {
               real_t* vdata  = cell->getData( microMesh->p1Mesh()->component( i ).getCellDataID() )->getPointer( level );
               real_t  coordA = vdata[vertexdof::macrocell::index(
                   level, microVertexIndexA.x(), microVertexIndexA.y(), microVertexIndexA.z() )];
               real_t  coordB = vdata[vertexdof::macrocell::index(
                   level, microVertexIndexB.x(), microVertexIndexB.y(), microVertexIndexB.z() )];
               position( Eigen::Index( i ) ) = coordA + positionFactor * ( coordB - coordA );
            }
            else if ( microMesh->polynomialDegree() == 2 )
            {
               real_t* vdata = cell->getData( microMesh->p2Mesh()->component( i ).getVertexDoFFunction().getCellDataID() )
                                   ->getPointer( level );
               real_t* edata =
                   cell->getData( microMesh->p2Mesh()->component( i ).getEdgeDoFFunction().getCellDataID() )->getPointer( level );
               const auto edgeIndex       = edgedof::calcEdgeDoFIndex( microVertexIndexA, microVertexIndexB );
               const auto edgeOrientation = edgedof::calcEdgeDoFOrientation( microVertexIndexA, microVertexIndexB );

               real_t weightA = vdata[vertexdof::macrocell::index(
                   level, microVertexIndexA.x(), microVertexIndexA.y(), microVertexIndexA.z() )];
               real_t weightB =
                   edata[edgedof::macrocell::index( level, edgeIndex.x(), edgeIndex.y(), edgeIndex.z(), edgeOrientation )];
               real_t weightC = vdata[vertexdof::macrocell::index(
                   level, microVertexIndexB.x(), microVertexIndexB.y(), microVertexIndexB.z() )];

               const real_t xi = positionFactor;
               position( Eigen::Index( i ) ) =
                   ( real_c( 2 ) * weightA - real_c( 4 ) * weightB + real_c( 2 ) * weightC ) * xi * xi -
                   ( real_c( 3 ) * weightA - real_c( 4 ) * weightB + weightC ) * xi + weightA;
            }
            else if ( microMesh->polynomialDegree() == 3 )
            {
               const auto edgeIndex       = edgedof::calcEdgeDoFIndex( microVertexIndexA, microVertexIndexB );
               const auto edgeOrientation = edgedof::calcEdgeDoFOrientation( microVertexIndexA, microVertexIndexB );

               // maybe not ideal, but a way to check whether we can avoid interpolation
               if ( std::fpclassify( positionFactor - real_c( 1.0 / 3.0 ) ) == FP_ZERO )
               {
                  real_t* blueData = cell->getData( microMesh->p3Mesh()->component( i ).getEdgeDoFFunctionBlue().getCellDataID() )
                                         ->getPointer( level );
                  position( Eigen::Index( i ) ) =
                      blueData[edgedof::macrocell::index( level, edgeIndex.x(), edgeIndex.y(), edgeIndex.z(), edgeOrientation )];
               }
               else if ( std::fpclassify( positionFactor - real_c( 2.0 / 3.0 ) ) == FP_ZERO )
               {
                  real_t* redData = cell->getData( microMesh->p3Mesh()->component( i ).getEdgeDoFFunctionRed().getCellDataID() )
                                        ->getPointer( level );
                  position( Eigen::Index( i ) ) =
                      redData[edgedof::macrocell::index( level, edgeIndex.x(), edgeIndex.y(), edgeIndex.z(), edgeOrientation )];
               }
               else
               {
                  WALBERLA_CHECK( microMesh->polynomialDegree() == 3, "Full P3 interpolation not implemented, yet!" );
               }
            }
            else
            {
               WALBERLA_ABORT( "MicroMesh: this should not happen." );
            }
         }
      }
   }
   else
   {
      WALBERLA_ABORT( "MicroMesh: PrimitiveID not existing locally. " )
   }

   return position;
}

std::vector< std::function< real_t( const Point3D& ) > > microMeshMapFromGeometryMap( const GeometryMap& geometryMap )
{
   return { [&]( const Point3D& x ) {
              Point3D xout;
              geometryMap.evalF( x, xout );
              return xout[0];
           },
            [&]( const Point3D& x ) {
               Point3D xout;
               geometryMap.evalF( x, xout );
               return xout[1];
            },
            [&]( const Point3D& x ) {
               Point3D xout;
               geometryMap.evalF( x, xout );
               return xout[2];
            } };
}

/////////////////////////////////////////////
/// Computation of micro-center positions ///
/////////////////////////////////////////////

static Point3D microFaceCenterPositionNoMesh( const std::shared_ptr< PrimitiveStorage >& storage,
                                              PrimitiveID                                faceId,
                                              uint_t                                     level,
                                              const indexing::Index&                     microFaceIndex,
                                              facedof::FaceType                          faceType,
                                              bool                                       withBlending )
{
   WALBERLA_ASSERT( storage->faceExistsLocally( faceId ), "Cannot compute center of micro-face for non-local macro-face" );

   const Face& face = *( storage->getPrimitiveGenerically< Face >( faceId ) );

   std::array< indexing::Index, 3 > vertexIndices = facedof::macroface::getMicroVerticesFromMicroFace( microFaceIndex, faceType );

   Point3D v0 = vertexdof::macroface::coordinateFromIndex( level, face, vertexIndices[0] );
   Point3D v1 = vertexdof::macroface::coordinateFromIndex( level, face, vertexIndices[1] );
   Point3D v2 = vertexdof::macroface::coordinateFromIndex( level, face, vertexIndices[2] );

   Point3D center = ( v0 + v1 + v2 ) / real_c( 3 );

   if ( !withBlending )
   {
      return center;
   }

   return applyBlending( center, face );
}

Point3D microFaceCenterPosition( const std::shared_ptr< PrimitiveStorage >& storage,
                                 PrimitiveID                                faceId,
                                 uint_t                                     level,
                                 const indexing::Index&                     microFaceIndex,
                                 facedof::FaceType                          faceType )

{
   WALBERLA_CHECK( storage->faceExistsLocally( faceId ), "Cannot compute center of micro-face for non-local macro-face" );

   auto microMesh = storage->getMicroMesh();

   if ( !microMesh )
   {
      return microFaceCenterPositionNoMesh( storage, faceId, level, microFaceIndex, faceType, true );
   }

   const uint_t                  dim  = microMesh->dimension();
   std::shared_ptr< const Face > face = storage->getPrimitiveGenerically< Face >( faceId );
   Point3D                       center;

   if ( microMesh->polynomialDegree() == 1 )
   {
      const auto& p1Mesh = *( microMesh->p1Mesh() );

      std::array< uint_t, 3 > vertexDoFIndices;
      vertexdof::getVertexDoFDataIndicesFromMicroFace( microFaceIndex, faceType, level, vertexDoFIndices );

      std::array< Point3D, 3 > coords;

      for ( uint_t k = 0; k < dim; k++ )
      {
         real_t* vdata = face->getData( p1Mesh.component( k ).getFaceDataID() )->getPointer( level );
         for ( uint_t vIdx = 0; vIdx < 3; ++vIdx )
         {
            center( Eigen::Index( k ) ) += vdata[vertexDoFIndices[vIdx]];
         }
         center( Eigen::Index( k ) ) /= real_t( 3 );
      }
   }
   else if ( microMesh->polynomialDegree() == 2 )
   {
      const auto& p2Mesh = *( microMesh->p2Mesh() );

      std::array< uint_t, 3 > vertexDoFIndices;
      vertexdof::getVertexDoFDataIndicesFromMicroFace( microFaceIndex, faceType, level, vertexDoFIndices );

      std::array< uint_t, 3 > edgeDoFIndices;
      edgedof::getEdgeDoFDataIndicesFromMicroFaceFEniCSOrdering( microFaceIndex, faceType, level, edgeDoFIndices );

      std::array< Point3D, 6 > node;

      for ( uint_t k = 0; k < dim; k++ )
      {
         real_t* vData = face->getData( p2Mesh.component( k ).getVertexDoFFunction().getFaceDataID() )->getPointer( level );
         real_t* eData = face->getData( p2Mesh.component( k ).getEdgeDoFFunction().getFaceDataID() )->getPointer( level );

         real_t node0 = vData[vertexDoFIndices[0]];
         real_t node1 = vData[vertexDoFIndices[1]];
         real_t node2 = vData[vertexDoFIndices[2]];

         real_t node3 = eData[edgeDoFIndices[0]];
         real_t node4 = eData[edgeDoFIndices[1]];
         real_t node5 = eData[edgeDoFIndices[2]];

         center( Eigen::Index( k ) ) = ( real_c( 4 ) * ( node3 + node4 + node5 ) - ( node0 + node1 + node2 ) ) / real_c( 9 );
      }
   }
   else if ( microMesh->polynomialDegree() == 3 )
   {
      for ( uint_t k = 0; k < dim; k++ )
      {
         const auto memLayout = microMesh->p3Mesh()->component( k ).getFaceDoFFunction().memoryLayout();
         auto       dofs      = microMesh->p3Mesh()->component( k ).getFaceDoFFunction().dofMemory( faceId, level );

         center( Eigen::Index( k ) ) =
             dofs[volumedofspace::indexing::index( microFaceIndex.x(), microFaceIndex.y(), faceType, 0, 1, level, memLayout )];
      }
   }
   else
   {
      WALBERLA_ABORT( "Invalid polynomial degree " << microMesh->polynomialDegree() << " of MicroMesh." )
   }

   return center;
}

static Point3D microCellCenterPositionNoMesh( const std::shared_ptr< PrimitiveStorage >& storage,
                                              PrimitiveID                                cellId,
                                              uint_t                                     level,
                                              const indexing::Index&                     microCellIndex,
                                              celldof::CellType                          cellType,
                                              bool                                       withBlending )
{
   WALBERLA_ASSERT( storage->cellExistsLocally( cellId ), "Cannot compute center of micro-cell for non-local macro-cell" );

   const Cell& cell = *( storage->getPrimitiveGenerically< Cell >( cellId ) );

   std::array< indexing::Index, 4 > vertexIndices = celldof::macrocell::getMicroVerticesFromMicroCell( microCellIndex, cellType );

   Point3D v0 = vertexdof::macrocell::coordinateFromIndex( level, cell, vertexIndices[0] );
   Point3D v1 = vertexdof::macrocell::coordinateFromIndex( level, cell, vertexIndices[1] );
   Point3D v2 = vertexdof::macrocell::coordinateFromIndex( level, cell, vertexIndices[2] );
   Point3D v3 = vertexdof::macrocell::coordinateFromIndex( level, cell, vertexIndices[3] );

   Point3D center = ( v0 + v1 + v2 + v3 ) / real_c( 4 );

   if ( !withBlending )
   {
      return center;
   }

   return applyBlending( center, cell );
}

Point3D microCellCenterPosition( const std::shared_ptr< PrimitiveStorage >& storage,
                                 PrimitiveID                                cellId,
                                 uint_t                                     level,
                                 const indexing::Index&                     microCellIndex,
                                 celldof::CellType                          cellType )

{
   WALBERLA_CHECK( storage->cellExistsLocally( cellId ), "Cannot compute center of micro-cell for non-local macro-cell" );

   auto microMesh = storage->getMicroMesh();

   if ( !microMesh )
   {
      return microCellCenterPositionNoMesh( storage, cellId, level, microCellIndex, cellType, true );
   }
   else
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "Centroid computation is not accurate in 3D when parametric map is active" );

      const Cell& cell = *( storage->getPrimitiveGenerically< Cell >( cellId ) );

      std::array< indexing::Index, 4 > vertexIndices =
          celldof::macrocell::getMicroVerticesFromMicroCell( microCellIndex, cellType );

      Point3D v0 = micromesh::microVertexPosition( storage, cellId, level, vertexIndices[0] );
      Point3D v1 = micromesh::microVertexPosition( storage, cellId, level, vertexIndices[1] );
      Point3D v2 = micromesh::microVertexPosition( storage, cellId, level, vertexIndices[2] );
      Point3D v3 = micromesh::microVertexPosition( storage, cellId, level, vertexIndices[3] );

      Point3D center = ( v0 + v1 + v2 + v3 ) / real_c( 4 );

      return microCellCenterPositionNoMesh( storage, cellId, level, microCellIndex, cellType, false );
   }
}

} // namespace hyteg::micromesh
