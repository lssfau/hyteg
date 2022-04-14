/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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

#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"

#include "hyteg/HytegDefinitions.hpp"
#include "hyteg/Levelinfo.hpp"
#include "hyteg/gridtransferoperators/generatedKernels/restrict_2D_macroface_P1_pull_additive.hpp"
#include "hyteg/gridtransferoperators/generatedKernels/restrict_3D_macrocell_P1_pull_additive.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"

namespace hyteg {

void P1toP1LinearRestriction::restrict2D( const P1Function< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const
{
  const uint_t destinationLevel  = sourceLevel - 1;
  const auto   storage           = function.getStorage();
  const auto   boundaryCondition = function.getBoundaryCondition();

  function.communicate< Vertex, Edge >( sourceLevel );
  function.communicate< Edge, Face >( sourceLevel );
  function.communicate< Face, Edge >( sourceLevel );
  function.communicate< Edge, Vertex >( sourceLevel );

  for( const auto& it : storage->getVertices() )
  {
    const Vertex& vertex  = *it.second;
    const auto    srcData = vertex.getData( function.getVertexDataID() )->getPointer( sourceLevel );
    auto          dstData = vertex.getData( function.getVertexDataID() )->getPointer( destinationLevel );

    if( testFlag( boundaryCondition.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
    {
      restrictMacroVertex( srcData, dstData, sourceLevel, vertex.getNumNeighborEdges() );
    }
  }

  for( const auto& it : storage->getEdges() )
  {
    const Edge& edge    = *it.second;
    const auto  srcData = edge.getData( function.getEdgeDataID() )->getPointer( sourceLevel );
    auto        dstData = edge.getData( function.getEdgeDataID() )->getPointer( destinationLevel );

    if( testFlag( boundaryCondition.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
    {
      restrictMacroEdge( srcData, dstData, sourceLevel, edge.getNumNeighborFaces() );
    }
  }

  for( const auto& it : storage->getFaces() )
  {
    const Face& face    = *it.second;
    const auto  srcData = face.getData( function.getFaceDataID() )->getPointer( sourceLevel );
    auto        dstData = face.getData( function.getFaceDataID() )->getPointer( destinationLevel );

    if( testFlag( boundaryCondition.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
    {
      restrictMacroFace( srcData, dstData, sourceLevel, 0 );
    }
  }
}

static real_t calculateInverseFactorToScaleNeighborhoodContribution(
  const std::array< real_t, 4 > & invNumNeighborsOfVertex,
  const std::array< real_t, 6 > & invNumNeighborsOfEdge,
  const std::array< real_t, 4 > invNumNeighborsOfFace,
  const indexing::Index & index,
  const uint_t & level )
{
  const auto stencilLeafOnCellVertices = vertexdof::macrocell::isOnCellVertex( index, level );
  const auto stencilLeafOnCellEdges    = vertexdof::macrocell::isOnCellEdge  ( index, level );
  const auto stencilLeafOnCellFaces    = vertexdof::macrocell::isOnCellFace  ( index, level );

  real_t invFactorDueToNeighborhood = real_c( 1 );

  if ( stencilLeafOnCellVertices.size() > 0 )
  {
    WALBERLA_ASSERT_EQUAL( stencilLeafOnCellVertices.size(), 1 );
    const auto localVertexID = *stencilLeafOnCellVertices.begin();
    invFactorDueToNeighborhood = invNumNeighborsOfVertex[ localVertexID ];
  }
  else if ( stencilLeafOnCellEdges.size() > 0 )
  {
    WALBERLA_ASSERT_EQUAL( stencilLeafOnCellEdges.size(), 1 );
    const auto localEdgeID = *stencilLeafOnCellEdges.begin();
    invFactorDueToNeighborhood = invNumNeighborsOfEdge[ localEdgeID ];
  }
  else if ( stencilLeafOnCellFaces.size() > 0 )
  {
    WALBERLA_ASSERT_EQUAL( stencilLeafOnCellFaces.size(), 1 );
    const auto localFaceID = *stencilLeafOnCellFaces.begin();
    invFactorDueToNeighborhood = invNumNeighborsOfFace[ localFaceID ];
  }

  return invFactorDueToNeighborhood;
}

void P1toP1LinearRestriction::restrict2DAdditively( const P1Function< real_t >& function,
                                                    const uint_t&               sourceLevel,
                                                    const DoFType&              flag ) const
{
   /// XOR flag with all to get the DoFTypes that should be excluded
  const DoFType excludeFlag = (flag ^ All);

  const uint_t destinationLevel = sourceLevel - 1;

  function.communicate< Vertex, Edge >  ( sourceLevel );
  function.communicate< Edge,   Face >  ( sourceLevel );

  auto storage = function.getStorage();

  for ( const auto & faceIt : function.getStorage()->getFaces() )
  {
    const auto face = faceIt.second;
    const auto srcData = face->getData( function.getFaceDataID())->getPointer( sourceLevel );
    auto dstData = face->getData( function.getFaceDataID())->getPointer( destinationLevel );

    const auto numNeighborFacesEdge0 =
    static_cast< double >( storage->getEdge( face->neighborEdges().at( 0 ))->getNumNeighborFaces());
    const auto numNeighborFacesEdge1 =
    static_cast< double >( storage->getEdge( face->neighborEdges().at( 1 ))->getNumNeighborFaces());
    const auto numNeighborFacesEdge2 =
    static_cast< double >( storage->getEdge( face->neighborEdges().at( 2 ))->getNumNeighborFaces());
    const auto numNeighborFacesVertex0 =
    static_cast< double >( storage->getVertex( face->neighborVertices().at( 0 ))->getNumNeighborFaces());
    const auto numNeighborFacesVertex1 =
    static_cast< double >( storage->getVertex( face->neighborVertices().at( 1 ))->getNumNeighborFaces());
    const auto numNeighborFacesVertex2 =
    static_cast< double >( storage->getVertex( face->neighborVertices().at( 2 ))->getNumNeighborFaces());

    vertexdof::macroface::generated::restrict_2D_macroface_P1_pull_additive( dstData,
                                                                             srcData,
                                                                             static_cast< int32_t >( destinationLevel ),
                                                                             numNeighborFacesEdge0,
                                                                             numNeighborFacesEdge1,
                                                                             numNeighborFacesEdge2,
                                                                             numNeighborFacesVertex0,
                                                                             numNeighborFacesVertex1,
                                                                             numNeighborFacesVertex2 );
  }

  function.communicateAdditively< Face, Edge >( destinationLevel, excludeFlag, *storage );
  function.communicateAdditively< Face, Vertex >( destinationLevel, excludeFlag, *storage );
}

void P1toP1LinearRestriction::restrict3D( const P1Function< real_t >& function,
                                          const uint_t&               sourceLevel,
                                          const DoFType&              flag ) const
{
   /// XOR flag with all to get the DoFTypes that should be excluded
  const DoFType excludeFlag = (flag ^ All);

  const uint_t destinationLevel = sourceLevel - 1;

  function.communicate< Vertex, Edge >  ( sourceLevel );
  function.communicate< Edge,   Face >  ( sourceLevel );
  function.communicate< Face,   Cell >  ( sourceLevel );

  for ( const auto & cellIt : function.getStorage()->getCells() )
  {
    const auto cell = cellIt.second;
    const auto srcData = cell->getData( function.getCellDataID())->getPointer( sourceLevel );
    auto dstData = cell->getData( function.getCellDataID())->getPointer( destinationLevel );

    if ( globalDefines::useGeneratedKernels )
    {
       auto storage = function.getStorage();

       const double numNeighborCellsFace0 =
           static_cast< double >( storage->getFace( cell->neighborFaces().at( 0 ) )->getNumNeighborCells() );
       const double numNeighborCellsFace1 =
           static_cast< double >( storage->getFace( cell->neighborFaces().at( 1 ) )->getNumNeighborCells() );
       const double numNeighborCellsFace2 =
           static_cast< double >( storage->getFace( cell->neighborFaces().at( 2 ) )->getNumNeighborCells() );
       const double numNeighborCellsFace3 =
           static_cast< double >( storage->getFace( cell->neighborFaces().at( 3 ) )->getNumNeighborCells() );

       const double numNeighborCellsEdge0 =
           static_cast< double >( storage->getEdge( cell->neighborEdges().at( 0 ) )->getNumNeighborCells() );
       const double numNeighborCellsEdge1 =
           static_cast< double >( storage->getEdge( cell->neighborEdges().at( 1 ) )->getNumNeighborCells() );
       const double numNeighborCellsEdge2 =
           static_cast< double >( storage->getEdge( cell->neighborEdges().at( 2 ) )->getNumNeighborCells() );
       const double numNeighborCellsEdge3 =
           static_cast< double >( storage->getEdge( cell->neighborEdges().at( 3 ) )->getNumNeighborCells() );
       const double numNeighborCellsEdge4 =
           static_cast< double >( storage->getEdge( cell->neighborEdges().at( 4 ) )->getNumNeighborCells() );
       const double numNeighborCellsEdge5 =
           static_cast< double >( storage->getEdge( cell->neighborEdges().at( 5 ) )->getNumNeighborCells() );

       const double numNeighborCellsVertex0 =
           static_cast< double >( storage->getVertex( cell->neighborVertices().at( 0 ) )->getNumNeighborCells() );
       const double numNeighborCellsVertex1 =
           static_cast< double >( storage->getVertex( cell->neighborVertices().at( 1 ) )->getNumNeighborCells() );
       const double numNeighborCellsVertex2 =
           static_cast< double >( storage->getVertex( cell->neighborVertices().at( 2 ) )->getNumNeighborCells() );
       const double numNeighborCellsVertex3 =
           static_cast< double >( storage->getVertex( cell->neighborVertices().at( 3 ) )->getNumNeighborCells() );

       vertexdof::macrocell::generated::restrict_3D_macrocell_P1_pull_additive( dstData,
                                                                                srcData,
                                                                                static_cast< int32_t >( destinationLevel ),
                                                                                numNeighborCellsEdge0,
                                                                                numNeighborCellsEdge1,
                                                                                numNeighborCellsEdge2,
                                                                                numNeighborCellsEdge3,
                                                                                numNeighborCellsEdge4,
                                                                                numNeighborCellsEdge5,
                                                                                numNeighborCellsFace0,
                                                                                numNeighborCellsFace1,
                                                                                numNeighborCellsFace2,
                                                                                numNeighborCellsFace3,
                                                                                numNeighborCellsVertex0,
                                                                                numNeighborCellsVertex1,
                                                                                numNeighborCellsVertex2,
                                                                                numNeighborCellsVertex3 );
    }

    else
    {
      // Calculate inverse number of neighboring cells for each neighboring macro-primitive.
      std::array< real_t, 4 > invNumNeighborsOfVertex;
      std::array< real_t, 6 > invNumNeighborsOfEdge;
      std::array< real_t, 4 > invNumNeighborsOfFace;

      invNumNeighborsOfVertex.fill( real_c( 0 ));
      invNumNeighborsOfEdge.fill( real_c( 0 ));
      invNumNeighborsOfFace.fill( real_c( 0 ));

      for ( const auto & neighborVertexID : cell->neighborVertices())
      {
        invNumNeighborsOfVertex[cell->getLocalVertexID( neighborVertexID )] = real_c( 1 ) / real_c( function.getStorage()->getVertex( neighborVertexID )->getNumNeighborCells());
      }
      for ( const auto & neighborEdgeID : cell->neighborEdges())
      {
        invNumNeighborsOfEdge[cell->getLocalEdgeID( neighborEdgeID )] = real_c( 1 ) / real_c( function.getStorage()->getEdge( neighborEdgeID )->getNumNeighborCells());
      }
      for ( const auto & neighborFaceID : cell->neighborFaces())
      {
        invNumNeighborsOfFace[cell->getLocalFaceID( neighborFaceID )] = real_c( 1 ) / real_c( function.getStorage()->getFace( neighborFaceID )->getNumNeighborCells());
      }

      for ( const auto & dstIdx : vertexdof::macrocell::Iterator( destinationLevel ))
      {
        const auto srcIdx = dstIdx * 2;
        const auto arrayIdxDst = vertexdof::macrocell::index( destinationLevel, dstIdx.x(), dstIdx.y(), dstIdx.z());
        const auto arrayIdxSrcCenter = vertexdof::macrocell::index( sourceLevel, srcIdx.x(), srcIdx.y(), srcIdx.z());

        const auto onCellVertices = vertexdof::macrocell::isOnCellVertex( dstIdx, destinationLevel );
        const auto onCellEdges = vertexdof::macrocell::isOnCellEdge( dstIdx, destinationLevel );
        const auto onCellFaces = vertexdof::macrocell::isOnCellFace( dstIdx, destinationLevel );

        // add center with weight one and scale depending on location of dst unknown
        const auto invFactorToScaleContributionCenter = calculateInverseFactorToScaleNeighborhoodContribution( invNumNeighborsOfVertex, invNumNeighborsOfEdge, invNumNeighborsOfFace,
                                                                                                               dstIdx, destinationLevel );

        dstData[arrayIdxDst] = invFactorToScaleContributionCenter * srcData[arrayIdxSrcCenter];

        // add leaves with weight .5 and scale depending on location of dst unknown
        if ( onCellVertices.size() > 0 )
        {
          WALBERLA_ASSERT_EQUAL( onCellVertices.size(), 1 );
          const auto localVertexID = *onCellVertices.begin();

          for ( const auto & dir : vertexdof::macrocell::neighborsOnVertexWithoutCenter[localVertexID] )
          {
            const auto arrayIdxSrcDir = vertexdof::macrocell::indexFromVertex( sourceLevel, srcIdx.x(), srcIdx.y(), srcIdx.z(), dir );
            dstData[arrayIdxDst] += 0.5 * invFactorToScaleContributionCenter * srcData[arrayIdxSrcDir];
          }
        } else if ( onCellEdges.size() > 0 )
        {
          WALBERLA_ASSERT_EQUAL( onCellEdges.size(), 1 );
          const auto localEdgeID = *onCellEdges.begin();

          for ( const auto & dir : vertexdof::macrocell::neighborsOnEdgeWithoutCenter[localEdgeID] )
          {
            const auto arrayIdxSrcDir = vertexdof::macrocell::indexFromVertex( sourceLevel, srcIdx.x(), srcIdx.y(), srcIdx.z(), dir );
            dstData[arrayIdxDst] += 0.5 * invFactorToScaleContributionCenter * srcData[arrayIdxSrcDir];
          }
        } else if ( onCellFaces.size() > 0 )
        {
          WALBERLA_ASSERT_EQUAL( onCellFaces.size(), 1 );
          const auto localFaceID = *onCellFaces.begin();

          for ( const auto & dir : vertexdof::macrocell::neighborsOnFaceWithoutCenter[localFaceID] )
          {
            const auto arrayIdxSrcDir = vertexdof::macrocell::indexFromVertex( sourceLevel, srcIdx.x(), srcIdx.y(), srcIdx.z(), dir );
            dstData[arrayIdxDst] += 0.5 * invFactorToScaleContributionCenter * srcData[arrayIdxSrcDir];
          }
        } else
        {
          for ( const auto & dir : vertexdof::macrocell::neighborsWithoutCenter )
          {
            const auto arrayIdxSrcDir = vertexdof::macrocell::indexFromVertex( sourceLevel, srcIdx.x(), srcIdx.y(), srcIdx.z(), dir );
            dstData[arrayIdxDst] += 0.5 * invFactorToScaleContributionCenter * srcData[arrayIdxSrcDir];
          }
        }
      }
    }
  }

  function.communicateAdditively< Cell, Vertex >( destinationLevel, excludeFlag, *function.getStorage() );
  function.communicateAdditively< Cell, Edge >( destinationLevel, excludeFlag, *function.getStorage() );
  function.communicateAdditively< Cell, Face >( destinationLevel, excludeFlag, *function.getStorage() );
}


void P1toP1LinearRestriction::restrictMacroVertex( const real_t *src, real_t *dst, const uint_t & sourceLevel,
                                                   const uint_t & numNeighborEdges ) const
{
  WALBERLA_UNUSED( sourceLevel );
  dst[0] = src[0];

  for (uint_t i = 0; i < numNeighborEdges; ++i) {
    dst[0] += 0.5*src[i+1];
    i += 1;
  }
}

void P1toP1LinearRestriction::restrictMacroEdge( const real_t *src, real_t *dst, const uint_t & sourceLevel,
                                                 const uint_t & numNeighborFaces ) const
{
  size_t rowsize_c = levelinfo::num_microvertices_per_edge( sourceLevel - 1 );

  idx_t i_c;
  for ( i_c = 1; i_c < idx_t( rowsize_c ) - 1; ++i_c )
  {
    dst[vertexdof::macroedge::indexFromVertex( sourceLevel - 1, i_c, stencilDirection::VERTEX_C )] = 1.0 * src[vertexdof::macroedge::indexFromVertex( sourceLevel,
                                                                                                                                                      2 * i_c,
                                                                                                                                                      stencilDirection::VERTEX_C )];

    for ( auto & neighbor : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF )
    {
      dst[vertexdof::macroedge::indexFromVertex( sourceLevel - 1, i_c, stencilDirection::VERTEX_C )] += 0.5 * src[vertexdof::macroedge::indexFromVertex( sourceLevel,
                                                                                                                                                         2 * i_c,
                                                                                                                                                         neighbor )];
    }

    for ( auto & neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF )
    {
      dst[vertexdof::macroedge::indexFromVertex( sourceLevel - 1, i_c, stencilDirection::VERTEX_C )] += 0.5 * src[vertexdof::macroedge::indexFromVertex( sourceLevel,
                                                                                                                                                         2 * i_c,
                                                                                                                                                         neighbor )];
    }

    if ( numNeighborFaces == 2 )
    {
      for ( auto & neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF )
      {
        dst[vertexdof::macroedge::indexFromVertex( sourceLevel - 1, i_c, stencilDirection::VERTEX_C )] += 0.5 * src[vertexdof::macroedge::indexFromVertex( sourceLevel,
                                                                                                                                                           2 *
                                                                                                                                                           i_c,
                                                                                                                                                           neighbor )];
      }
    }
  }
}

void P1toP1LinearRestriction::restrictMacroFace( const real_t *src, real_t *dst, const uint_t & sourceLevel,
                                                 const uint_t & numNeighborCells ) const
{
  WALBERLA_UNUSED( numNeighborCells );
  uint_t N_c = levelinfo::num_microvertices_per_edge( sourceLevel - 1 );
  uint_t N_c_i = N_c;

  real_t tmp;

  for ( idx_t j = 1; j < idx_t( N_c ) - 2; ++j )
  {
    for ( idx_t i = 1; i < idx_t( N_c_i ) - 2; ++i )
    {

      tmp = src[vertexdof::macroface::indexFromVertex( sourceLevel, 2 * i, 2 * j,
                                                       stencilDirection::VERTEX_C )];

      for ( const auto & neighbor : vertexdof::macroface::neighborsWithoutCenter )
      {
        tmp += 0.5 * src[vertexdof::macroface::indexFromVertex( sourceLevel, 2 * i, 2 * j, neighbor )];
      }

      dst[vertexdof::macroface::indexFromVertex( sourceLevel - 1, i, j, stencilDirection::VERTEX_C )] = tmp;
    }

    --N_c_i;
  }
}

}
