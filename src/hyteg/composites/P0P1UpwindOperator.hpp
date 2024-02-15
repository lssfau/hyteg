/*
 * Copyright (c) 2017-2023 Daniel Drzisga, Dominik Thoennes, Nils Kohl,
 * Marcus Mohr.
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
#pragma once

#include <array>

#include "hyteg/operators/Operator.hpp"
#include "hyteg/p0functionspace/P0Function.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/types/PointND.hpp"

namespace hyteg {

using volumedofspace::indexing::ElementNeighborInfo;
using walberla::real_c;
using Point = Point3D;

// \note: The P0P1UpwindOperator cannot handle boundary conditions, yet, and only works for 2D
class P0P1UpwindOperator : public Operator< P0Function< real_t >, P0Function< real_t > >
{
 public:
   P0P1UpwindOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                       const P1VectorFunction< real_t >&          velocity,
                       uint_t                                     minLevel,
                       uint_t                                     maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , velocity_( velocity )
   {
      WALBERLA_ASSERT( velocity.getDimension() == 2 );
   }

   ~P0P1UpwindOperator() override = default;

   void apply( const P0Function< real_t >& src,
               const P0Function< real_t >& dst,
               uint_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const
   {
      // need to update faces on the ghost-layers
      src.communicate( level );

      velocity_[0].startCommunication< Edge, Vertex >( level );
      velocity_[1].startCommunication< Edge, Vertex >( level );

      velocity_[0].startCommunication< Face, Edge >( level );
      velocity_[1].startCommunication< Face, Edge >( level );

      velocity_[0].endCommunication< Edge, Vertex >( level );
      velocity_[1].endCommunication< Edge, Vertex >( level );

      velocity_[0].endCommunication< Face, Edge >( level );
      velocity_[1].endCommunication< Face, Edge >( level );

      // it is sufficient to only work with the macro-faces
      for ( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if ( testFlag( faceBC, flag ) )
         {
            treatMacroFace( level, face, storage_, src, dst, updateType );
         }
      }

      // I do not think that we need this, we always assume that functions are dirty before using them
      // dst.communicate( level );
   }

 private:
   P1VectorFunction< real_t > velocity_;

   void treatMacroFace( const uint_t&                              level,
                        Face&                                      face,
                        const std::shared_ptr< PrimitiveStorage >& storage,
                        const P0Function< real_t >&                src,
                        const P0Function< real_t >&                dst,
                        UpdateType                                 updateType ) const
   {
      using namespace vertexdof::macroface;

      // get velocity memories
      auto u = face.getData( velocity_[0].getFaceDataID() )->getPointer( level );
      auto v = face.getData( velocity_[1].getFaceDataID() )->getPointer( level );

      // get p0 memories
      const auto faceID       = face.getID();
      const auto srcMemLayout = src.getDGFunction()->volumeDoFFunction()->memoryLayout();
      const auto dstMemLayout = dst.getDGFunction()->volumeDoFFunction()->memoryLayout();
      auto       srcDoFs      = src.getDGFunction()->volumeDoFFunction()->dofMemory( faceID, level );
      auto       dstDoFs      = dst.getDGFunction()->volumeDoFFunction()->dofMemory( faceID, level );

      real_t faceArea    = std::pow( real_c( 4 ), -real_c( level ) ) * face.getArea();
      real_t faceAreaInv = real_c( 1 ) / faceArea;

      Point uEdge;

      uint_t offsetToCenter = 0u;
      static bool didWarnAlready = false;  // warn only once about Dirichlet BC per macro face

      // got blue and gray triangles
      for ( uint_t microVolType = 0; microVolType < 2; microVolType++ )
      {
         auto faceType = facedof::allFaceTypes[microVolType];

         for ( auto itFace = facedof::macroface::Iterator( level, faceType, offsetToCenter ).begin(); itFace != itFace.end();
               ++itFace )
         {
            Index               elementIdx = *itFace;
            ElementNeighborInfo neighborInfo =
                ElementNeighborInfo( elementIdx, faceType, level, src.getBoundaryCondition(), faceID, storage_ );

            // for accumulating the new value to assign or add to dst
            real_t aux{ real_c( 0 ) };

            // loop over the three micro-edges
            for ( uint_t k = 0; k < 3; ++k )
            {
               const std::vector< Point >& vertex     = neighborInfo.interfaceVertexCoords( k );
               real_t                      edgeLength = ( vertex[1] - vertex[0] ).norm();
               auto                        normal     = neighborInfo.outwardNormal( k );

               // evaluate velocity on edge
               const std::vector< Index > vertexID = neighborInfo.interfaceVertexIndices( k );
               uEdge[0] =
                   real_c( 0.5 ) * ( u[indexFromVertex( level, vertexID[0].x(), vertexID[0].y(), stencilDirection::VERTEX_C )] +
                                     u[indexFromVertex( level, vertexID[1].x(), vertexID[1].y(), stencilDirection::VERTEX_C )] );
               uEdge[1] =
                   real_c( 0.5 ) * ( v[indexFromVertex( level, vertexID[0].x(), vertexID[0].y(), stencilDirection::VERTEX_C )] +
                                     v[indexFromVertex( level, vertexID[1].x(), vertexID[1].y(), stencilDirection::VERTEX_C )] );
               uEdge[2] = real_c( 0.0 );

               real_t normalFlow = edgeLength * uEdge.dot( normal );
               if ( normalFlow >= real_c( 0 ) )
               {
                  aux += normalFlow * srcDoFs[volumedofspace::indexing::index(
                                          elementIdx.x(), elementIdx.y(), faceType, 0, 1, level, srcMemLayout )];
               }
               else
               {
                  // case 1: neighbour face is not on ghost layer
                  if ( !neighborInfo.atMacroBoundary( k ) )
                  {
                     auto nbrID   = neighborInfo.neighborElementIndices( k );
                     auto nbrType = neighborInfo.neighborFaceType( k );
                     aux += normalFlow *
                            srcDoFs[volumedofspace::indexing::index( nbrID.x(), nbrID.y(), nbrType, 0, 1, level, srcMemLayout )];
                  }

                  // case 2: neighbour face is on ghost layer (no domain boundary)
                  else if ( neighborInfo.neighborBoundaryType( k ) == Inner )
                  {
                     uint_t idx = volumedofspace::indexing::indexNeighborInGhostLayer(
                         neighborInfo.macroBoundaryID( k ), elementIdx.x(), elementIdx.y(), faceType, 0, 1, level, srcMemLayout );
                     real_t* glMemory =
                         src.getDGFunction()->volumeDoFFunction()->glMemory( faceID, level, neighborInfo.macroBoundaryID( k ) );

                     aux += normalFlow * glMemory[idx];
                  }

                  // case 3: no neighbour face because we are at a domain boundary
                  else
                  {
                     // need to implement something reasonable here
                     if ( neighborInfo.neighborBoundaryType( k ) == DirichletBoundary && !didWarnAlready )
                     {
                        WALBERLA_LOG_WARNING_ON_ROOT( "P0P1UpwindOperator ignores Dirichlet boundary condition!" );
                        didWarnAlready = true;
                     }
                  }
               }
            }

            if ( updateType == Replace )
            {
               dstDoFs[volumedofspace::indexing::index( elementIdx.x(), elementIdx.y(), faceType, 0, 1, level, dstMemLayout )] =
                   aux * faceAreaInv;
            }
            else if ( updateType == Add )
            {
               dstDoFs[volumedofspace::indexing::index( elementIdx.x(), elementIdx.y(), faceType, 0, 1, level, dstMemLayout )] +=
                   aux * faceAreaInv;
            }
         }
      }
   }
};

} // namespace hyteg
