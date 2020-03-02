/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Marcus Mohr.
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

#include <cmath>

#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

#include "GeometryMap.hpp"

#define SHELL_MAP_DEBUG
#define SHELL_MAP_LOG( STR ) WALBERLA_LOG_INFO_ON_ROOT( STR );

namespace hyteg {

/// Class providing geometry mapping for a facetted isosahedral shell
///
///
///
///
///
///
class IcosahedralShellMap : public GeometryMap
{
 public:
   IcosahedralShellMap( const Cell& cell, const SetupPrimitiveStorage& storage )
   {
      SHELL_MAP_LOG( "---------------------------------------------" );
      SHELL_MAP_LOG( "Initialising Shell map for cellID: " << cell.getID() );
      SHELL_MAP_LOG( "---------------------------------------------" );
      classifyVertices( cell, storage );
   }

   IcosahedralShellMap( const Face& face, const SetupPrimitiveStorage& storage )
   {
      SHELL_MAP_LOG( "---------------------------------------------" );
      SHELL_MAP_LOG( "Initialising Shell map for faceID: " << face.getID() );
      SHELL_MAP_LOG( "---------------------------------------------" );

      std::vector< PrimitiveID > neighborCells;
      face.getNeighborCells( neighborCells );
      WALBERLA_ASSERT_GREATER( neighborCells.size(), 0 );
      const Cell& cell = *storage.getCell( neighborCells[0] );
      classifyVertices( cell, storage );
   }

   IcosahedralShellMap( const Edge& edge, const SetupPrimitiveStorage& storage )
   {
      SHELL_MAP_LOG( "---------------------------------------------" );
      SHELL_MAP_LOG( "Initialising Shell map for edgeID: " << edge.getID() );
      SHELL_MAP_LOG( "---------------------------------------------" );

      std::vector< PrimitiveID > neighborCells;
      edge.getNeighborCells( neighborCells );
      WALBERLA_ASSERT_GREATER( neighborCells.size(), 0 );
      const Cell& cell = *storage.getCell( neighborCells[0] );
      classifyVertices( cell, storage );
   }

   void evalF( const Point3D& xold, Point3D& xnew ) const
   {
      // determine barycentric coordinate w.r.t. vertex refVertex_
      real_t tmp0  = -rayVertex_[2];
      real_t tmp1  = refVertex_[2] + tmp0;
      real_t tmp2  = -rayVertex_[0];
      real_t tmp3  = thrVertex_[0] + tmp2;
      real_t tmp4  = -rayVertex_[1];
      real_t tmp5  = forVertex_[1] + tmp4;
      real_t tmp6  = tmp3 * tmp5;
      real_t tmp7  = refVertex_[1] + tmp4;
      real_t tmp8  = thrVertex_[2] + tmp0;
      real_t tmp9  = forVertex_[0] + tmp2;
      real_t tmp10 = tmp8 * tmp9;
      real_t tmp11 = refVertex_[0] + tmp2;
      real_t tmp12 = thrVertex_[1] + tmp4;
      real_t tmp13 = forVertex_[2] + tmp0;
      real_t tmp14 = tmp12 * tmp13;
      real_t tmp15 = tmp13 * tmp3;
      real_t tmp16 = tmp12 * tmp9;
      real_t tmp17 = tmp5 * tmp8;
      real_t tmp18 = tmp0 + xold[2];
      real_t tmp19 = tmp4 + xold[1];
      real_t tmp20 = tmp2 + xold[0];

      real_t volT = -tmp1 * tmp16 + tmp1 * tmp6 + tmp10 * tmp7 + tmp11 * tmp14 - tmp11 * tmp17 - tmp15 * tmp7;
      real_t volX = tmp10 * tmp19 + tmp14 * tmp20 - tmp15 * tmp19 - tmp16 * tmp18 - tmp17 * tmp20 + tmp18 * tmp6;
      real_t bary = std::abs( volX / volT );

      // compute new coordinates
      real_t oldRad = std::sqrt( xold[0] * xold[0] + xold[1] * xold[1] + xold[2] * xold[2] );
      real_t newRad = radRayVertex_ + bary * ( radRefVertex_ - radRayVertex_ );
      xnew[0]       = xold[0] * newRad / oldRad;
      xnew[1]       = xold[1] * newRad / oldRad;
      xnew[2]       = xold[2] * newRad / oldRad;

      // SHELL_MAP_LOG( "Mapped: " << xold << " --> " << xnew );
   }

   void evalDF( const Point3D& x, Matrix2r& DFx ) const { WALBERLA_ABORT( "evalDF unimplemented!" ); }

   void evalDFinv( const Point3D& x, Matrix2r& DFinvx ) const { WALBERLA_ABORT( "evalDFinv unimplemented!" ); }

   void serializeSubClass( walberla::mpi::SendBuffer& sendBuffer ) const
   {
      sendBuffer << Type::ICOSAHEDRAL_SHELL_MAP << rayVertex_ << refVertex_ << thrVertex_ << forVertex_ << radRefVertex_
                 << radRayVertex_;
   }

   static void setMap( SetupPrimitiveStorage& setupStorage )
   {
      for ( auto it : setupStorage.getCells() )
      {
         Cell& cell = *it.second;
         setupStorage.setGeometryMap( cell.getID(), std::make_shared< IcosahedralShellMap >( cell, setupStorage ) );
      }

      for ( auto it : setupStorage.getFaces() )
      {
         Face& face = *it.second;
         setupStorage.setGeometryMap( face.getID(), std::make_shared< IcosahedralShellMap >( face, setupStorage ) );
      }

      for ( auto it : setupStorage.getEdges() )
      {
         Edge& edge = *it.second;
         setupStorage.setGeometryMap( edge.getID(), std::make_shared< IcosahedralShellMap >( edge, setupStorage ) );
      }
   }

 private:
   /// \name Classified vertices of macro triangle
   ///
   /// Each macro triangle of the annulus has two vertices which lie on a ray coming from the origin and
   /// two with the same distance from the origin. The vertex opposite to the edge formed by the latter
   /// two is stored as rayVertex_, the one on the ray with the rayVertex_ is stored as refVertex_.
   /// The third vertex then is stored as thrVertex_.
   ///@{
   Point3D refVertex_;
   Point3D rayVertex_;
   Point3D thrVertex_;
   Point3D forVertex_;
   ///@}

   /// tolerance for comparing numerical values for equality
   const real_t tol = 1e-14;

   /// distance from origin of vertex rayVertex_
   real_t radRefVertex_;

   /// distance from origin of vertex refVertex_
   real_t radRayVertex_;

   /// internal enumeration class for classifying tetrahedra
   enum class tetType
   {
      TET_INWARDS,
      TET_OUTWARDS,
      TET_SKEW
   };

   /// method for classifying the vertices of the macro tetrahedron
   tetType classifyTet( const Cell& cell, std::array< real_t, 4 >& radius )
   {
      const std::array< Point3D, 4 >& coords = cell.getCoordinates();
      real_t                          innerRad, outerRad;
      innerRad = std::numeric_limits< real_t >::max();
      outerRad = 0.0;

      for ( uint_t k = 0; k < 4; k++ )
      {
         radius[k] = std::sqrt( coords[k].normSq() );
         innerRad  = radius[k] < innerRad ? radius[k] : innerRad;
         outerRad  = radius[k] > outerRad ? radius[k] : outerRad;
      }

      SHELL_MAP_LOG( "outer radius = " << outerRad );
      SHELL_MAP_LOG( "inner radius = " << innerRad );

      uint_t nOuterNodes = 0;
      for ( uint_t k = 0; k < 4; k++ )
      {
         SHELL_MAP_LOG( "radius[" << k << "] = " << radius[k] );
         nOuterNodes += ( outerRad - radius[k] ) < tol ? 1 : 0;
      }

      SHELL_MAP_LOG( "classifyTet: nOuterNodes = " << nOuterNodes );
      tetType thisTetType;
      switch ( nOuterNodes )
      {
      case 1:
         thisTetType = tetType::TET_OUTWARDS;
         SHELL_MAP_LOG( " -> TET_OUTWARDS" );
         break;
      case 2:
         thisTetType = tetType::TET_SKEW;
         SHELL_MAP_LOG( " -> TET_SKEW" );
         break;
      case 3:
         thisTetType = tetType::TET_INWARDS;
         SHELL_MAP_LOG( " -> TET_INWARDS" );
         break;
      default:
         WALBERLA_ABORT( "Houston we have a problem! Cannot classify macro tetrahedron!" );
      }

      return thisTetType;
   }

   void classifyVertices( const Cell& cell, const SetupPrimitiveStorage& storage )
   {
      const std::array< Point3D, 4 >& coords = cell.getCoordinates();

      real_t aux;
      real_t pairFound = false;

      uint_t idxRefVertex;
      uint_t idxRayVertex;
      uint_t idxThrVertex;
      uint_t idxForVertex;

      SHELL_MAP_LOG( "micro-vertex 0 = " << coords[0] );
      SHELL_MAP_LOG( "micro-vertex 1 = " << coords[1] );
      SHELL_MAP_LOG( "micro-vertex 2 = " << coords[2] );
      SHELL_MAP_LOG( "micro-vertex 3 = " << coords[3] );

      // determine type of macro-tet
      std::array< real_t, 4 > radius;
      tetType                 thisTetType = classifyTet( cell, radius );

      // determine the two vertices lying on a radial ray
      for ( uint_t k = 0; k < 4 && !pairFound; k++ )
      {
         for ( uint_t j = k + 1; j < 4 && !pairFound; j++ )
         {
            SHELL_MAP_LOG( "Testing cross product v[" << k << "] x v[" << j << "]" );

            // x-component of cross-product
            aux = coords[k][1] * coords[j][2] - coords[k][2] * coords[j][1];
            if ( std::abs( aux ) > tol )
               continue;

            // y-component of cross-product
            aux = coords[k][2] * coords[j][0] - coords[k][0] * coords[j][2];
            if ( std::abs( aux ) > tol )
               continue;

            // z-component of cross-product
            aux = coords[k][0] * coords[j][1] - coords[k][1] * coords[j][0];
            if ( std::abs( aux ) > tol )
               continue;

            // still here, so we found the vertex pair (ordering does not matter)
            idxRefVertex = k;
            idxRayVertex = j;
            pairFound    = true;
         }
      }

      if ( !pairFound )
      {
         WALBERLA_ABORT( "Error in finding vertex pair on radial ray!!!" );
      }

      // ------------------------------------------------------------------
      //  for skew tets we have a problem, so we need to find another tet
      //  from the same prims that is non-skew
      // ------------------------------------------------------------------
      if ( thisTetType == tetType::TET_SKEW )
      {
         PrimitiveID altCellID = findNonSkewTetInPrism( cell, storage, idxRefVertex, idxRayVertex );
         classifyVertices( *storage.getCell( altCellID ), storage );
      }

      else
      {
         // remember the indices of the two remaining nodes
         for ( uint_t k = 0; k < 4; k++ )
         {
            if ( k != idxRefVertex && k != idxRayVertex )
            {
               idxThrVertex = k;
            }
         }
         for ( uint_t k = 0; k < 4; k++ )
         {
            if ( k != idxRefVertex && k != idxRayVertex && k != idxThrVertex )
            {
               idxForVertex = k;
            }
         }

         // now sort ref and ray vertices depending on tet type
         switch ( thisTetType )
         {
         case tetType::TET_OUTWARDS:
            if ( radius[idxRefVertex] < radius[idxRayVertex] )
            {
               uint_t swp   = idxRayVertex;
               idxRayVertex = idxRefVertex;
               idxRefVertex = swp;
            }
            break;

         case tetType::TET_INWARDS:
            if ( radius[idxRefVertex] > radius[idxRayVertex] )
            {
               uint_t swp   = idxRayVertex;
               idxRayVertex = idxRefVertex;
               idxRefVertex = swp;
            }
            break;

         default:
            WALBERLA_ABORT( "We should not have ended up in this case!" );
         }

         // store vertex coordinates and radii
         refVertex_ = coords[idxRefVertex];
         rayVertex_ = coords[idxRayVertex];
         thrVertex_ = coords[idxThrVertex];
         forVertex_ = coords[idxForVertex];

         radRefVertex_ = std::sqrt( refVertex_.normSq() );
         radRayVertex_ = std::sqrt( rayVertex_.normSq() );

         SHELL_MAP_LOG( "refVertex = " << refVertex_ );
         SHELL_MAP_LOG( "rayVertex = " << rayVertex_ );
         SHELL_MAP_LOG( "thrVertex = " << thrVertex_ );
         SHELL_MAP_LOG( "forVertex = " << forVertex_ );
      }
   }

   /// For a "skew" tetrahedron find a non-skew neighbour from the same prism
   PrimitiveID
       findNonSkewTetInPrism( const Cell& cell, const SetupPrimitiveStorage& storage, uint_t& idxRefVertex, uint_t& idxRayVertex )
   {
      std::vector< PrimitiveID > verts;
      cell.getNeighborVertices( verts );

      PrimitiveID rayNode = verts[idxRayVertex];
      PrimitiveID refNode = verts[idxRefVertex];

      WALBERLA_ASSERT_EQUAL( idxRayVertex, cell.getLocalVertexID( rayNode ) );
      WALBERLA_ASSERT_EQUAL( idxRefVertex, cell.getLocalVertexID( refNode ) );

      std::vector< PrimitiveID > nbrFaces;
      cell.getNeighborFaces( nbrFaces );
      SHELL_MAP_LOG( "Skew cell has " << nbrFaces.size() << " face neighbours" );

      for ( uint_t k = 0; k < nbrFaces.size(); k++ )
      {
         const Face* candidate = storage.getFace( nbrFaces[k] );
         candidate->getNeighborVertices( verts );

         // check that not both vertices are part of the face
         auto itRef = std::find( verts.begin(), verts.end(), refNode );
         auto itRay = std::find( verts.begin(), verts.end(), rayNode );
         if ( !( itRef != verts.end() && itRay != verts.end() ) )
         {
            SHELL_MAP_LOG( "Found fitting face: id = " << nbrFaces[k] );

            // now select the correct cell
            std::vector< PrimitiveID > nbrCells;
            candidate->getNeighborCells( nbrCells );
            SHELL_MAP_LOG( "  Candidate face has " << nbrCells.size() << " cell neighbours" );
            for ( uint_t j = 0; j < nbrCells.size(); j++ )
            {
               if ( nbrCells[j] != cell.getID() )
               {
                  SHELL_MAP_LOG( "    Cell we need has ID = " << nbrCells[j] );
                  return nbrCells[j];
               }
            }
         }
      }
   }
};

} // namespace hyteg
