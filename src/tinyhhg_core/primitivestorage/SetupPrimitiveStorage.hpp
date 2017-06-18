
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

  SetupPrimitiveStorage( const MeshInfo & meshInfo );

  void toStream( std::ostream & os ) const;

  void getVertices( std::map< PrimitiveID::IDType, SetupVertex* > & vertices ) const { vertices = vertices_; }
  void getEdges( std::map< PrimitiveID::IDType, SetupEdge* > & edges ) const { edges = edges_; }
  void getFaces( std::map< PrimitiveID::IDType, SetupFace* > & faces ) const { faces = faces_; }

  /// Searches an edge with the respective vertices by ID\n
  /// \param edge is set to the ID of the edge if one was found
  /// \return true, if an edge was found, false otherwise
  bool findEdgeByVertexIDs( const PrimitiveID & vertexID0, const PrimitiveID & vertexID1, PrimitiveID & edge ) const;

private:

  PrimitiveID generatePrimitiveID();

  void checkIDConsistency() {};

  std::map< PrimitiveID::IDType, SetupVertex* > vertices_;
  std::map< PrimitiveID::IDType, SetupEdge*  >  edges_;
  std::map< PrimitiveID::IDType, SetupFace* >   faces_;

};

inline std::ostream & operator<<( std::ostream & os, const SetupPrimitiveStorage & storage )
{
  storage.toStream( os );
  return os;
}

} // namespace hhg
