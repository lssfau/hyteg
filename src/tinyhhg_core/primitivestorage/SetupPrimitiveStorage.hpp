
#pragma once

#include "core/debug/Debug.h"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitiveid.hpp"
#include "tinyhhg_core/primitives/SetupEdge.hpp"
#include "tinyhhg_core/primitives/SetupFace.hpp"
#include "tinyhhg_core/primitives/SetupVertex.hpp"

#include <map>
#include <set>
#include <tuple>
#include <vector>

namespace hhg {

class SetupPrimitiveStorage
{
public:

  /// Function definition for loadbalancing callbacks
  /// \param storage the \ref SetupPrimitiveStorage to be balanced
  /// \param numberOfProcesses the overall number of processes available
  /// \param perProcessMemoryLimit the main memory limit per process
  /// \return number of processes that have been assigned at least one \ref Primitive, \n
  ///         this number can be less than numberOfProcesses
  typedef std::function< uint_t ( SetupPrimitiveStorage & storage,
				  const memory_t & perProcessMemoryLimit ) >
          TargetProcessAssignmentFunction;

  typedef std::map< PrimitiveID::IDType, SetupPrimitive* > SetupPrimitiveMap;
  typedef std::map< PrimitiveID::IDType, SetupVertex* >    SetupVertexMap;
  typedef std::map< PrimitiveID::IDType, SetupEdge* >      SetupEdgeMap;
  typedef std::map< PrimitiveID::IDType, SetupFace* >      SetupFaceMap;

  SetupPrimitiveStorage( const MeshInfo & meshInfo, const uint_t & numberOfProcesses );

  void toStream( std::ostream & os ) const;

  uint_t getNumberOfProcesses() const { return numberOfProcesses_; }
  uint_t getNumberOfEmptyProcesses() const;

  bool primitiveExists( const PrimitiveID & id ) const { return vertexExists( id ) || edgeExists( id ) || faceExists( id ); }
  bool vertexExists   ( const PrimitiveID & id ) const { return vertices_.count( id.getID() ) > 0; }
  bool edgeExists     ( const PrimitiveID & id ) const { return edges_.count( id.getID() )    > 0; }
  bool faceExists     ( const PrimitiveID & id ) const { return faces_.count( id.getID() )    > 0; }

  const SetupVertex * getVertex( const PrimitiveID & id ) const { return vertexExists( id ) ? vertices_.at( id.getID() ) : nullptr; }
  const SetupEdge   * getEdge  ( const PrimitiveID & id ) const { return edgeExists( id )   ? edges_.at( id.getID() )    : nullptr; }
  const SetupFace   * getFace  ( const PrimitiveID & id ) const { return faceExists( id )   ? faces_.at( id.getID() )    : nullptr; }

  void getSetupPrimitives( SetupPrimitiveMap & setupPrimitiveMap ) const;
  uint_t getNumberOfPrimitives() const;

  SetupVertexMap::iterator beginVertices() { return vertices_.begin(); }
  SetupVertexMap::iterator endVertices()   { return vertices_.end(); }

  SetupEdgeMap::iterator beginEdges()      { return edges_.begin(); }
  SetupEdgeMap::iterator endEdges()        { return edges_.end(); }

  SetupFaceMap::iterator beginFaces()      { return faces_.begin(); }
  SetupFaceMap::iterator endFaces()        { return faces_.end(); }

  SetupVertexMap::const_iterator beginVertices() const { return vertices_.begin(); }
  SetupVertexMap::const_iterator endVertices()   const { return vertices_.end(); }

  SetupEdgeMap::const_iterator beginEdges()      const { return edges_.begin(); }
  SetupEdgeMap::const_iterator endEdges()        const { return edges_.end(); }

  SetupFaceMap::const_iterator beginFaces()      const { return faces_.begin(); }
  SetupFaceMap::const_iterator endFaces()        const { return faces_.end(); }

  /// Searches an edge with the respective vertices by ID\n
  /// \param edge is set to the ID of the edge if one was found
  /// \return true, if an edge was found, false otherwise
  bool findEdgeByVertexIDs( const PrimitiveID & vertexID0, const PrimitiveID & vertexID1, PrimitiveID & edge ) const;

  void balanceLoad( const TargetProcessAssignmentFunction & loadbalanceCallback,
		    const memory_t & perProcessMemoryLimit );

  void   setTargetRank( const PrimitiveID & primitiveID, const uint_t & targetRank )       { primitiveIDToTargetRankMap_[ primitiveID.getID() ] = targetRank; }
  uint_t getTargetRank( const PrimitiveID & primitiveID )                            const { return primitiveIDToTargetRankMap_.at( primitiveID.getID() ); }

private:

  typedef std::map< uint_t, std::vector< PrimitiveID::IDType > > RankToSetupPrimitivesMap;

  PrimitiveID generatePrimitiveID() const;

  void assembleRankToSetupPrimitivesMap( RankToSetupPrimitivesMap & rankToSetupPrimitivesMap ) const;

  /// Returns the number of primitives on the target rank with the least number of primitives
  uint_t getMinPrimitivesPerRank() const;
  /// Returns the number of primitives on the target rank with the largest number of primitives
  uint_t getMaxPrimitivesPerRank() const;
  /// Returns the average number of primitives per rank
  real_t getAvgPrimitivesPerRank() const;

  uint_t numberOfProcesses_;

  SetupVertexMap vertices_;
  SetupEdgeMap   edges_;
  SetupFaceMap   faces_;

  std::map< PrimitiveID::IDType, uint_t > primitiveIDToTargetRankMap_;

};

inline std::ostream & operator<<( std::ostream & os, const SetupPrimitiveStorage & storage )
{
  storage.toStream( os );
  return os;
}

} // namespace hhg
