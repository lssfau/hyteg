
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
  ///         this number can be less than numberOfProcesses if some processes do not obtain a primitive
  typedef std::function< uint_t ( SetupPrimitiveStorage & storage, const uint_t & numberOfProcesses,
				  const memory_t & perProcessMemoryLimit ) >
          TargetProcessAssignmentFunction;

  typedef std::map< PrimitiveID::IDType, SetupPrimitive* > SetupPrimitiveMap;
  typedef std::map< PrimitiveID::IDType, SetupVertex* >    SetupVertexMap;
  typedef std::map< PrimitiveID::IDType, SetupEdge* >      SetupEdgeMap;
  typedef std::map< PrimitiveID::IDType, SetupFace* >      SetupFaceMap;

  SetupPrimitiveStorage( const MeshInfo & meshInfo );

  void toStream( std::ostream & os ) const;

  void getSetupPrimitives( SetupPrimitiveMap & setupPrimitiveMap ) const;

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
		    const uint_t & numberOfProcesses,
		    const memory_t & perProcessMemoryLimit );

private:

  PrimitiveID generatePrimitiveID();

  SetupVertexMap vertices_;
  SetupEdgeMap   edges_;
  SetupFaceMap   faces_;

};

inline std::ostream & operator<<( std::ostream & os, const SetupPrimitiveStorage & storage )
{
  storage.toStream( os );
  return os;
}

} // namespace hhg
