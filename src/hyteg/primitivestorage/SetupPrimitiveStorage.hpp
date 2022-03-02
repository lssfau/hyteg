/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include <map>
#include <set>
#include <tuple>
#include <vector>

#include "core/debug/Debug.h"

#include "hyteg/PrimitiveID.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitives/Cell.hpp"
#include "hyteg/primitives/Edge.hpp"
#include "hyteg/primitives/Face.hpp"
#include "hyteg/primitives/Primitive.hpp"
#include "hyteg/primitives/Vertex.hpp"

namespace hyteg {

using walberla::memory_t;
using walberla::real_t;

class SetupPrimitiveStorage
{
 public:
   typedef std::map< PrimitiveID::IDType, std::shared_ptr< Primitive > > PrimitiveMap;
   typedef std::map< PrimitiveID::IDType, std::shared_ptr< Vertex > >    VertexMap;
   typedef std::map< PrimitiveID::IDType, std::shared_ptr< Edge > >      EdgeMap;
   typedef std::map< PrimitiveID::IDType, std::shared_ptr< Face > >      FaceMap;
   typedef std::map< PrimitiveID::IDType, std::shared_ptr< Cell > >      CellMap;

   SetupPrimitiveStorage( const MeshInfo& meshInfo, const uint_t& numberOfProcesses );

   SetupPrimitiveStorage( const VertexMap& vertices,
                          const EdgeMap&   edges,
                          const FaceMap&   faces,
                          const CellMap&   cells,
                          const uint_t&    numberOfProcesses );

   void toStream( std::ostream& os, bool verbose = false ) const;

   uint_t getNumberOfProcesses() const { return numberOfProcesses_; }
   uint_t getNumberOfEmptyProcesses() const;

   bool primitiveExists( const PrimitiveID& id ) const
   {
      return vertexExists( id ) || edgeExists( id ) || faceExists( id ) || cellExists( id );
   }
   bool vertexExists( const PrimitiveID& id ) const { return vertices_.count( id.getID() ) > 0; }
   bool edgeExists( const PrimitiveID& id ) const { return edges_.count( id.getID() ) > 0; }
   bool faceExists( const PrimitiveID& id ) const { return faces_.count( id.getID() ) > 0; }
   bool cellExists( const PrimitiveID& id ) const { return cells_.count( id.getID() ) > 0; }

   Primitive*       getPrimitive( const PrimitiveID& id );
   const Primitive* getPrimitive( const PrimitiveID& id ) const;
   Vertex*          getVertex( const PrimitiveID& id ) { return vertexExists( id ) ? vertices_.at( id.getID() ).get() : nullptr; }
   const Vertex*    getVertex( const PrimitiveID& id ) const
   {
      return vertexExists( id ) ? vertices_.at( id.getID() ).get() : nullptr;
   }
   Edge*       getEdge( const PrimitiveID& id ) { return edgeExists( id ) ? edges_.at( id.getID() ).get() : nullptr; }
   const Edge* getEdge( const PrimitiveID& id ) const { return edgeExists( id ) ? edges_.at( id.getID() ).get() : nullptr; }
   Face*       getFace( const PrimitiveID& id ) { return faceExists( id ) ? faces_.at( id.getID() ).get() : nullptr; }
   const Face* getFace( const PrimitiveID& id ) const { return faceExists( id ) ? faces_.at( id.getID() ).get() : nullptr; }
   Cell*       getCell( const PrimitiveID& id ) { return cellExists( id ) ? cells_.at( id.getID() ).get() : nullptr; }
   const Cell* getCell( const PrimitiveID& id ) const { return cellExists( id ) ? cells_.at( id.getID() ).get() : nullptr; }

   void getSetupPrimitives( PrimitiveMap& setupPrimitiveMap ) const;

   uint_t getNumberOfPrimitives() const
   {
      return getNumberOfVertices() + getNumberOfEdges() + getNumberOfFaces() + getNumberOfCells();
   }
   uint_t getNumberOfVertices() const { return vertices_.size(); }
   uint_t getNumberOfEdges() const { return edges_.size(); }
   uint_t getNumberOfFaces() const { return faces_.size(); }
   uint_t getNumberOfCells() const { return cells_.size(); }

   /// Returns a reference to a map of \ref Vertex instances
   const VertexMap& getVertices() const { return vertices_; }

   /// Returns a reference to a map of \ref Edge instances
   const EdgeMap& getEdges() const { return edges_; }

   /// Returns a reference to a map of \ref Face instances
   const FaceMap& getFaces() const { return faces_; }

   /// Returns a reference to a map of \ref Cell instances
   const CellMap& getCells() const { return cells_; }

   void setTargetRank( const PrimitiveID& primitiveID, const uint_t& targetRank )
   {
      primitiveIDToTargetRankMap_[primitiveID.getID()] = targetRank;
   }
   uint_t getTargetRank( const PrimitiveID& primitiveID ) const { return primitiveIDToTargetRankMap_.at( primitiveID.getID() ); }

   uint_t getNumCellsOnRank( uint_t rank ) const;
   uint_t getNumFacesOnRank( uint_t rank ) const;
   uint_t getNumEdgesOnRank( uint_t rank ) const;
   uint_t getNumVerticesOnRank( uint_t rank ) const;

   void setGeometryMap( const PrimitiveID& primitiveID, const std::shared_ptr< GeometryMap >& map )
   {
      getPrimitive( primitiveID )->geometryMap_ = map;
   }

   void setMeshBoundaryFlag( const PrimitiveID& primitiveID, const uint_t& meshBoundaryFlag )
   {
      getPrimitive( primitiveID )->meshBoundaryFlag_ = meshBoundaryFlag;
   }

   /// Sets the mesh boundary flag of the primitives to a specified value if they are located on the boundary of the domain
   /// \param meshBoundaryFlagOnBoundary the flag the primitives are set to if they are located on the boundary
   /// \param meshBoundaryFlagInner the flag the primitives are set to if they are not located on the boundary
   /// \param highestDimensionAlwaysInner if true, cells in 3D meshes and faces in 2D meshes are treated as inner primitives
   void setMeshBoundaryFlagsOnBoundary( const uint_t& meshBoundaryFlagOnBoundary,
                                        const uint_t& meshBoundaryFlagInner,
                                        const bool&   highestDimensionAlwaysInner );

   /// Sets the mesh boundary flag of the primitives to a specified value if they are NOT located on the boundary of the domain.
   /// This does not change mesh boundary flags on primitives that are not located on the boundary
   /// \param meshBoundaryFlagInner the flag the primitives are set to if they are not located on the boundary
   /// \param highestDimensionAlwaysInner if true, cells in 3D meshes and faces in 2D meshes are treated as inner primitives
   void setMeshBoundaryFlagsInner( const uint_t& meshBoundaryFlagInner, const bool& highestDimensionAlwaysInner );

   /// Every primitive for which onBoundary() returns true for all / any (if allVertices == true / false) of the primitve's vertices is assigned the passed mesh boundary flag.
   void setMeshBoundaryFlagsByVertexLocation( const uint_t&                                    meshBoundaryFlag,
                                              const std::function< bool( const Point3D& x ) >& onBoundary,
                                              const bool&                                      allVertices = true );

   /// Every primitive for which onBoundary() returns true for the primitive's centroid is assigned the passed mesh boundary flag.
   void setMeshBoundaryFlagsByCentroidLocation( const uint_t&                                    meshBoundaryFlag,
                                                const std::function< bool( const Point3D& x ) >& onBoundary,
                                                bool                                             useGeometryMap = true );

   /// Returns true, if the primitive lies on the boundary.
   /// \param primitiveID the ID of the primitive to be tested
   /// \param highestDimensionAlwaysInner if true, this method always returns false if the targeted primitive is of highest dimension (cells for 3D, faces for 2D meshes)
   bool onBoundary( const PrimitiveID& primitiveID, const bool& highestDimensionAlwaysInner = false ) const;
   /// Returns the number of primitives that lie on the boundary
   uint_t getNumVerticesOnBoundary() const;
   uint_t getNumEdgesOnBoundary() const;
   uint_t getNumFacesOnBoundary() const;
   uint_t getNumCellsOnBoundary() const;

 private:
   typedef std::map< uint_t, std::vector< PrimitiveID::IDType > > RankToSetupPrimitivesMap;

   PrimitiveID generatePrimitiveID() const;

   void assembleRankToSetupPrimitivesMap( RankToSetupPrimitivesMap& rankToSetupPrimitivesMap ) const;

   /// Returns the number of primitives on the target rank with the least number of primitives
   uint_t getMinPrimitivesPerRank() const;
   /// Returns the number of primitives on the target rank with the largest number of primitives
   uint_t getMaxPrimitivesPerRank() const;
   /// Returns the average number of primitives per rank
   real_t getAvgPrimitivesPerRank() const;

   uint_t numberOfProcesses_;

   VertexMap vertices_;
   EdgeMap   edges_;
   FaceMap   faces_;
   CellMap   cells_;

   std::map< PrimitiveID::IDType, uint_t > primitiveIDToTargetRankMap_;
};

inline std::ostream& operator<<( std::ostream& os, const SetupPrimitiveStorage& storage )
{
   storage.toStream( os );
   return os;
}

} // namespace hyteg
