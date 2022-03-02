/*
 * Copyright (c) 2017-2019 Christoph Schwarzmeier, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include <array>
#include <vector>

#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"

#include "hyteg/mesh/MeshInfo.hpp"

using walberla::real_c;
using walberla::math::pi;

namespace hyteg {

MeshInfo MeshInfo::meshAnnulus( const real_t      rhoMin,
                                const real_t      rhoMax,
                                const real_t      phiLeft,
                                const real_t      phiRight,
                                const meshFlavour flavour,
                                uint_t            nTan,
                                uint_t            nRad )
{
   WALBERLA_ASSERT_LESS( rhoMin, rhoMax );
   WALBERLA_ASSERT_LESS( phiLeft, phiRight );
   WALBERLA_ASSERT_GREATER( rhoMin, 1e-8 );
   WALBERLA_ASSERT_GREATER( nTan, 0 );
   WALBERLA_ASSERT_GREATER( nRad, 0 );

   // mesh partial annulus in polar coordinates
   MeshInfo meshInfo =
       MeshInfo::meshRectangle( Point2D( {rhoMin, phiLeft} ), Point2D( {rhoMax, phiRight} ), flavour, nRad, nTan );

   // map vertex coordinates to cartesian domain
   Point3D node;
   node[2] = 0.0;
   uint_t boundaryFlag;
   real_t rho, phi;
   for ( size_t id = 0; id < meshInfo.vertices_.size(); ++id )
   {
      node                   = meshInfo.vertices_[id].getCoordinates();
      boundaryFlag           = meshInfo.vertices_[id].getBoundaryFlag();
      rho                    = node[0];
      phi                    = node[1];
      node[0]                = rho * std::cos( phi );
      node[1]                = rho * std::sin( phi );
      meshInfo.vertices_[id] = MeshInfo::Vertex( id, Point3D( node ), boundaryFlag );
   }

   return meshInfo;
}

MeshInfo MeshInfo::meshAnnulus( const real_t rhoMin, const real_t rhoMax, const meshFlavour flavour, uint_t nTan, uint_t nRad )
{
   WALBERLA_ASSERT_LESS( rhoMin, rhoMax );
   WALBERLA_ASSERT_GREATER( rhoMin, 1e-8 );
   WALBERLA_ASSERT_GREATER( nTan, 0 );
   WALBERLA_ASSERT_GREATER( nRad, 0 );

   // mesh a rectangle representing the annulus in polar coordinates
   MeshInfo meshInfo = MeshInfo::meshRectangle(
       Point2D( {rhoMin, real_c( 0.0 )} ), Point2D( {rhoMax, real_c( 2.0 ) * pi} ), flavour, nRad, nTan );

   // determine some tolerances for further operations
   const real_t tolFactor  = 0.1;
   const real_t xTolerance = ( ( rhoMax - rhoMin ) / real_c( nRad ) ) * tolFactor;
   const real_t yTolerance = ( 2.0 * pi / real_c( nTan ) ) * tolFactor;

   // set boundary flags for vertices;
   // note that the interior of the left and right edge of the rectangle will be glued together
   for ( auto& it : meshInfo.vertices_ )
   {
      real_t rho = it.second.getCoordinates()[0];
      if ( std::abs( rhoMax - rho ) < xTolerance )
      {
         it.second.setBoundaryFlag( flagOuterBoundary );
      }
      else if ( std::abs( rhoMin - rho ) < xTolerance )
      {
         it.second.setBoundaryFlag( flagInnerBoundary );
      }
      else
      {
         it.second.setBoundaryFlag( flagInterior );
      }
   }

   // now do the same for edges
   meshInfo.deduceEdgeFlagsFromVertices( flagInterior );

   // and (re)set flags for all faces
   for ( auto& it : meshInfo.faces_ )
   {
      it.second.setBoundaryFlag( flagInterior );
   }

   // close domain at phi = 0 = 2*pi:

   // Gather all vertices on the top boundary (phi = 2*pi) and matching
   // vertices on the bottom boundary (phi=0).
   std::map< IDType, Vertex > bottomBoundaryVertices; // map from top ID to bottom vertex instance
   for ( const auto& it : meshInfo.vertices_ )
   {
      auto vertex = it.second;

      if ( std::abs( vertex.getCoordinates()[1] - 2.0 * pi ) < yTolerance )
      {
         for ( const auto& itInner : meshInfo.vertices_ )
         {
            auto vertexBottom = itInner.second;

            if ( std::abs( vertexBottom.getCoordinates()[1] ) < yTolerance &&
                 std::abs( vertex.getCoordinates()[0] - vertexBottom.getCoordinates()[0] ) < xTolerance )
            {
               bottomBoundaryVertices[vertex.getID()] = vertexBottom;
            }
         }
      }
   }

   // Gather all top boundary edges from mesh.
   std::vector< std::array< IDType, 2 > > topBoundaryEdges;
   for ( const auto& it : meshInfo.edges_ )
   {
      auto edge = it.second;
      auto v0   = meshInfo.vertices_[edge.getVertices()[0]];
      auto v1   = meshInfo.vertices_[edge.getVertices()[1]];

      if ( std::abs( v0.getCoordinates()[1] - 2.0 * pi ) < yTolerance &&
           std::abs( v1.getCoordinates()[1] - 2.0 * pi ) < yTolerance )
      {
         topBoundaryEdges.push_back( {edge.getVertices()} );
      }
   }

   // Remove all top boundary vertices from mesh.
   for ( const auto& it : bottomBoundaryVertices )
   {
      meshInfo.vertices_.erase( it.first );
   }

   // Remove all top boundary edges from mesh.
   for ( const auto& it : topBoundaryEdges )
   {
      meshInfo.edges_.erase( it );
   }

   // Replace vertex IDs in all primitives.
   std::vector< std::array< IDType, 2 > >    edgesToRemove;
   std::map< std::array< IDType, 2 >, Edge > newEdges;

   for ( const auto& it : meshInfo.edges_ )
   {
      auto oldEdgeID = it.first;

      std::array< IDType, 2 > newKey      = oldEdgeID;
      bool                    replaceEdge = false;

      if ( bottomBoundaryVertices.count( oldEdgeID[0] ) > 0 )
      {
         newKey[0]   = bottomBoundaryVertices[oldEdgeID[0]].getID();
         replaceEdge = true;
      }
      if ( bottomBoundaryVertices.count( oldEdgeID[1] ) > 0 )
      {
         newKey[1]   = bottomBoundaryVertices[oldEdgeID[1]].getID();
         replaceEdge = true;
      }

      if ( replaceEdge )
      {
         edgesToRemove.push_back( it.first );
         Edge newEdge( newKey, it.second.getBoundaryFlag() );
         newEdges[newKey] = newEdge;
      }
   }

   for ( const auto& it : edgesToRemove )
   {
      meshInfo.edges_.erase( it );
   }

   for ( const auto& it : newEdges )
   {
      meshInfo.edges_[it.first] = it.second;
   }

   std::vector< std::vector< IDType > >    facesToRemove;
   std::map< std::vector< IDType >, Face > newFaces;

   for ( const auto& it : meshInfo.faces_ )
   {
      auto oldFaceID = it.first;

      std::vector< IDType > newKey      = oldFaceID;
      bool                  replaceFace = false;

      for ( uint_t i = 0; i < 3; i++ )
      {
         if ( bottomBoundaryVertices.count( oldFaceID[i] ) > 0 )
         {
            newKey[i]   = bottomBoundaryVertices[oldFaceID[i]].getID();
            replaceFace = true;
         }
      }

      if ( replaceFace )
      {
         facesToRemove.push_back( it.first );
         Face newFace( newKey, flagInterior );
         newFaces[newKey] = newFace;
      }
   }

   for ( const auto& it : facesToRemove )
   {
      meshInfo.faces_.erase( it );
   }

   for ( const auto& it : newFaces )
   {
      meshInfo.faces_[it.first] = it.second;
   }

   // map vertex coordinates to cartesian domain
   Point3D node;
   node[2] = 0.0;
   uint_t boundaryFlag;
   real_t rho, phi;
   for ( size_t id = 0; id < meshInfo.vertices_.size(); ++id )
   {
      node                   = meshInfo.vertices_[id].getCoordinates();
      boundaryFlag           = meshInfo.vertices_[id].getBoundaryFlag();
      rho                    = node[0];
      phi                    = node[1];
      node[0]                = rho * std::cos( phi );
      node[1]                = rho * std::sin( phi );
      meshInfo.vertices_[id] = MeshInfo::Vertex( id, Point3D( node ), boundaryFlag );
   }

   return meshInfo;
}

} // namespace hyteg
