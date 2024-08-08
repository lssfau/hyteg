/*
* Copyright (c) 2017-2024 Nils Kohl.
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
#include "hyteg/geometry/surfaces/Surface.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

namespace hyteg::micromesh {

static Point3D microVertexPositionNoMesh( const std::shared_ptr< PrimitiveStorage >& storage,
                                          PrimitiveID                                primitiveId,
                                          uint_t                                     level,
                                          const indexing::Index&                     microVertexIndex )
{
   if ( storage->cellExistsLocally( primitiveId ) )
   {
      return vertexdof::macrocell::coordinateFromIndex( level, *storage->getCell( primitiveId ), microVertexIndex );
   }
   else if ( storage->faceExistsLocally( primitiveId ) )
   {
      return vertexdof::macroface::coordinateFromIndex( level, *storage->getFace( primitiveId ), microVertexIndex );
   }
   else if ( storage->edgeExistsLocally( primitiveId ) )
   {
      return vertexdof::macroedge::coordinateFromIndex( level, *storage->getEdge( primitiveId ), microVertexIndex );
   }
   else if ( storage->vertexExistsLocally( primitiveId ) )
   {
      return storage->getVertex( primitiveId )->getCoordinates();
   }
   else
   {
      WALBERLA_ABORT( "MicroMesh: Primitive does not exist locally!" );
   }
}

static Point3D microEdgePositionNoMesh( const std::shared_ptr< PrimitiveStorage >& storage,
                                        PrimitiveID                                primitiveId,
                                        uint_t                                     level,
                                        const indexing::Index&                     microEdgeIndex,
                                        edgedof::EdgeDoFOrientation                microEdgeOrientation )
{
   if ( storage->cellExistsLocally( primitiveId ) )
   {
      return edgedof::macrocell::coordinateFromIndex(
          level, *storage->getCell( primitiveId ), microEdgeIndex, microEdgeOrientation );
   }
   else if ( storage->faceExistsLocally( primitiveId ) )
   {
      return edgedof::macroface::coordinateFromIndex(
          level, *storage->getFace( primitiveId ), microEdgeIndex, microEdgeOrientation );
   }
   else if ( storage->edgeExistsLocally( primitiveId ) )
   {
      return edgedof::macroedge::coordinateFromIndex( level, *storage->getEdge( primitiveId ), microEdgeIndex );
   }
   else
   {
      WALBERLA_ABORT( "MicroMesh: Primitive does not exist locally, or you passed a macro-vertex!" );
   }
}

static void initMicroMeshFromMacroMesh( P1VectorFunction< real_t >& p1Mesh, uint_t level )
{
   auto         storage   = p1Mesh.getStorage();
   const uint_t dimension = p1Mesh.getDimension();

   for ( const auto& [pid, vertex] : storage->getVertices() )
   {
      for ( uint_t i = 0; i < dimension; i++ )
      {
         vertex->getData( p1Mesh.component( i ).getVertexDataID() )->getPointer( level )[0] =
             microVertexPositionNoMesh( storage, pid, level, indexing::Index( 0, 0, 0 ) )( Eigen::Index( i ) );
      }
   }

   for ( const auto& [pid, edge] : storage->getEdges() )
   {
      for ( uint_t i = 0; i < dimension; i++ )
      {
         auto data = edge->getData( p1Mesh.component( i ).getEdgeDataID() )->getPointer( level );
         for ( auto idx : vertexdof::macroedge::Iterator( level ) )
         {
            data[vertexdof::macroedge::index( level, idx.x() )] =
                microVertexPositionNoMesh( storage, pid, level, idx )( Eigen::Index( i ) );
         }
      }
   }

   for ( const auto& [pid, face] : storage->getFaces() )
   {
      for ( uint_t i = 0; i < dimension; i++ )
      {
         auto data = face->getData( p1Mesh.component( i ).getFaceDataID() )->getPointer( level );
         for ( auto idx : vertexdof::macroface::Iterator( level ) )
         {
            data[vertexdof::macroface::index( level, idx.x(), idx.y() )] =
                microVertexPositionNoMesh( storage, pid, level, idx )( Eigen::Index( i ) );
         }
      }
   }

   for ( const auto& [pid, cell] : storage->getCells() )
   {
      for ( uint_t i = 0; i < dimension; i++ )
      {
         auto data = cell->getData( p1Mesh.component( i ).getCellDataID() )->getPointer( level );
         for ( auto idx : vertexdof::macrocell::Iterator( level ) )
         {
            data[vertexdof::macrocell::index( level, idx.x(), idx.y(), idx.z() )] =
                microVertexPositionNoMesh( storage, pid, level, idx )( Eigen::Index( i ) );
         }
      }
   }
}

static void initMicroMeshFromMacroMesh( P2VectorFunction< real_t >& p2Mesh, uint_t level )
{
   auto         storage   = p2Mesh.getStorage();
   const uint_t dimension = p2Mesh.getDimension();

   for ( const auto& [pid, vertex] : storage->getVertices() )
   {
      for ( uint_t i = 0; i < dimension; i++ )
      {
         vertex->getData( p2Mesh.component( i ).getVertexDoFFunction().getVertexDataID() )->getPointer( level )[0] =
             microVertexPositionNoMesh( storage, pid, level, indexing::Index( 0, 0, 0 ) )( Eigen::Index( i ) );
      }
   }

   for ( const auto& [pid, edge] : storage->getEdges() )
   {
      for ( uint_t i = 0; i < dimension; i++ )
      {
         auto vdata = edge->getData( p2Mesh.component( i ).getVertexDoFFunction().getEdgeDataID() )->getPointer( level );
         for ( auto idx : vertexdof::macroedge::Iterator( level ) )
         {
            vdata[vertexdof::macroedge::index( level, idx.x() )] =
                microVertexPositionNoMesh( storage, pid, level, idx )( Eigen::Index( i ) );
         }

         auto edata = edge->getData( p2Mesh.component( i ).getEdgeDoFFunction().getEdgeDataID() )->getPointer( level );
         for ( auto idx : edgedof::macroedge::Iterator( level ) )
         {
            edata[edgedof::macroedge::index( level, idx.x() )] =
                microEdgePositionNoMesh( storage, pid, level, idx, edgedof::EdgeDoFOrientation::X )( Eigen::Index( i ) );
         }
      }
   }

   for ( const auto& [pid, face] : storage->getFaces() )
   {
      for ( uint_t i = 0; i < dimension; i++ )
      {
         auto vdata = face->getData( p2Mesh.component( i ).getVertexDoFFunction().getFaceDataID() )->getPointer( level );
         for ( auto idx : vertexdof::macroface::Iterator( level ) )
         {
            vdata[vertexdof::macroface::index( level, idx.x(), idx.y() )] =
                microVertexPositionNoMesh( storage, pid, level, idx )( Eigen::Index( i ) );
         }

         auto edata = face->getData( p2Mesh.component( i ).getEdgeDoFFunction().getFaceDataID() )->getPointer( level );
         for ( auto idx : edgedof::macroface::Iterator( level ) )
         {
            for ( const auto& orientation : edgedof::faceLocalEdgeDoFOrientations )
            {
               edata[edgedof::macroface::index( level, idx.x(), idx.y(), orientation )] =
                   microEdgePositionNoMesh( storage, pid, level, idx, orientation )( Eigen::Index( i ) );
            }
         }
      }
   }

   for ( const auto& [pid, cell] : storage->getCells() )
   {
      for ( uint_t i = 0; i < dimension; i++ )
      {
         auto vdata = cell->getData( p2Mesh.component( i ).getVertexDoFFunction().getCellDataID() )->getPointer( level );
         for ( auto idx : vertexdof::macrocell::Iterator( level ) )
         {
            vdata[vertexdof::macrocell::index( level, idx.x(), idx.y(), idx.z() )] =
                microVertexPositionNoMesh( storage, pid, level, idx )( Eigen::Index( i ) );
         }

         auto edata = cell->getData( p2Mesh.component( i ).getEdgeDoFFunction().getCellDataID() )->getPointer( level );
         for ( auto idx : edgedof::macrocell::Iterator( level ) )
         {
            for ( const auto& orientation : edgedof::allEdgeDoFOrientationsWithoutXYZ )
            {
               edata[edgedof::macrocell::index( level, idx.x(), idx.y(), idx.z(), orientation )] =
                   microEdgePositionNoMesh( storage, pid, level, idx, orientation )( Eigen::Index( i ) );
            }
         }

         for ( auto idx : edgedof::macrocell::IteratorXYZ( level ) )
         {
            edata[edgedof::macrocell::index( level, idx.x(), idx.y(), idx.z(), edgedof::EdgeDoFOrientation::XYZ )] =
                microEdgePositionNoMesh( storage, pid, level, idx, edgedof::EdgeDoFOrientation::XYZ )( Eigen::Index( i ) );
         }
      }
   }
}

MicroMesh::MicroMesh( const std::shared_ptr< PrimitiveStorage >& storage,
                      uint_t                                     minLevel,
                      uint_t                                     maxLevel,
                      uint_t                                     polynomialDegree,
                      uint_t                                     dimension )
{
   WALBERLA_CHECK_GREATER_EQUAL( dimension, 2, "The only supported mesh dimensions are 2 and 3." )
   WALBERLA_CHECK_LESS_EQUAL( dimension, 3, "The only supported mesh dimensions are 2 and 3." )

   WALBERLA_CHECK_LESS_EQUAL( minLevel, maxLevel, "The MicroMesh's min level should be less or equal to its max level." )

   if ( polynomialDegree == 1 )
   {
      p1_ = std::make_shared< P1VectorFunction< real_t > >( "microMeshP1", storage, minLevel, maxLevel, dimension );
      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         initMicroMeshFromMacroMesh( *p1_, level );
      }
   }
   else if ( polynomialDegree == 2 )
   {
      p2_ = std::make_shared< P2VectorFunction< real_t > >( "microMeshP2", storage, minLevel, maxLevel, dimension );
      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         initMicroMeshFromMacroMesh( *p2_, level );
      }
   }
   else
   {
      WALBERLA_ABORT( "MicroMesh with polynomial degree " << polynomialDegree << " not supported." );
   }
}

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
                  const std::vector< std::function< real_t( const Point3D& ) > >& blendingFunction,
                  walberla::uint_t                                                level )
{
   WALBERLA_CHECK_GREATER_EQUAL(
       blendingFunction.size(),
       microMesh.dimension(),
       "MicroMesh: blending function is dimension is smaller than the dimension of the space the mesh is embedded in." )

   if ( microMesh.p1Mesh() )
   {
      microMesh.p1Mesh()->interpolate( blendingFunction, level );
   }
   else if ( microMesh.p2Mesh() )
   {
      microMesh.p2Mesh()->interpolate( blendingFunction, level );
   }
   else
   {
      WALBERLA_ABORT( "MicroMesh has no function allocated :/" )
   }
}

void interpolate( const std::shared_ptr< PrimitiveStorage >&                      storage,
                  const std::vector< std::function< real_t( const Point3D& ) > >& blendingFunction,
                  walberla::uint_t                                                level )
{
   auto microMesh = storage->getMicroMesh();

   WALBERLA_CHECK_NOT_NULLPTR( microMesh, "MicroMesh: Cannot interpolate if no mesh has been added to the PrimitiveStorage!" )

   interpolate( *microMesh, blendingFunction, level );
}

void interpolateAndCommunicate( MicroMesh&                                                      microMesh,
                                const std::vector< std::function< real_t( const Point3D& ) > >& blendingFunction,
                                walberla::uint_t                                                level )
{
   interpolate( microMesh, blendingFunction, level );
   communicate( microMesh, level );
}

void interpolateAndCommunicate( const std::shared_ptr< PrimitiveStorage >&                      storage,
                                const std::vector< std::function< real_t( const Point3D& ) > >& blendingFunction,
                                walberla::uint_t                                                level )
{
   interpolate( storage, blendingFunction, level );
   communicate( storage, level );
}

void interpolateRefinedCoarseMesh( MicroMesh& microMesh, uint_t level )
{
   if ( microMesh.p1Mesh() )
   {
      initMicroMeshFromMacroMesh( *microMesh.p1Mesh(), level );
   }
   else if ( microMesh.p2Mesh() )
   {
      initMicroMeshFromMacroMesh( *microMesh.p2Mesh(), level );
   }
}

void interpolateRefinedCoarseMesh( const std::shared_ptr< PrimitiveStorage >& storage, uint_t level )
{
   auto microMesh = storage->getMicroMesh();

   WALBERLA_CHECK_NOT_NULLPTR( microMesh, "MicroMesh: Cannot interpolate if no mesh has been added to the PrimitiveStorage!" )

   interpolateRefinedCoarseMesh( *microMesh, level );
}

Point3D microVertexPosition( const std::shared_ptr< PrimitiveStorage >& storage,
                             PrimitiveID                                primitiveId,
                             uint_t                                     level,
                             const indexing::Index&                     microVertexIndex )
{
   WALBERLA_CHECK( storage->primitiveExistsLocally( primitiveId ), "Cannot compute micro-mesh index of non-local primitive." );

   auto microMesh = storage->getMicroMesh();

   if ( !microMesh )
   {
      return microVertexPositionNoMesh( storage, primitiveId, level, microVertexIndex );
   }

   WALBERLA_CHECK( microMesh->polynomialDegree() == 1 || microMesh->polynomialDegree() == 2,
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
      return microEdgePositionNoMesh( storage, primitiveId, level, microEdgeIndex, microEdgeOrientation );
   }

   WALBERLA_CHECK( microMesh->polynomialDegree() == 1 || microMesh->polynomialDegree() == 2,
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

} // namespace hyteg::micromesh
