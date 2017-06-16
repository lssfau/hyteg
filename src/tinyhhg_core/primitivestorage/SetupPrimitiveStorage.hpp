
#pragma once

#include "core/debug/Debug.h"
#include "tinyhhg_core/types/pointnd.hpp"
#include "tinyhhg_core/types/flags.hpp"
#include "tinyhhg_core/primitiveid.hpp"
#include "tinyhhg_core/primitives/SetupEdge.hpp"
#include "tinyhhg_core/primitives/SetupFace.hpp"
#include "tinyhhg_core/primitives/SetupVertex.hpp"

#include <map>
#include <set>
#include <tuple>
#include <vector>

namespace hhg {

using walberla::uint_t;
using walberla::real_t;

/// \brief Contains information about a mesh
/// \author Daniel Drzisga (drzisga@ma.tum.de)
///         Nils Kohl (nils.kohl@fau.de)
///
/// The \ref MeshInfo class is used to store information about meshes that are
/// typically constructed by reading mesh files. It provides static methods to construct
/// meshed from files.
///
/// \ref MeshInfo instance can then be used to create domains via \ref SetupPrimitiveStorage.
class MeshInfo
{
public:

  /// Construct a MeshInfo from a file in Gmsh format
  static MeshInfo fromGmshFile( const std::string & meshFileName );

  /// Returns vertices of the mesh
  /// \return map that maps ( vertex ID -> vertex coordinate )
  const std::map< uint_t, Point3D > & getVertices() const { return vertices_; };

  /// Returns edges of the mesh
  /// \return set of a pair of \n
  ///         - a pair of vertex indices (matching indices from the vertices from \ref getVertices \n
  ///         - and the corresponding \ref DoFType
  const std::set< std::pair< std::pair< uint_t, uint_t >, DoFType > > & getEdges() const { return edges_; };

  /// Returns faces of the mesh
  /// \return set of 3-tuples of vertex indices from the vertices from \ref getVertices
  const std::set< std::tuple< uint_t, uint_t, uint_t > > & getFaces() const { return faces_; }

private:

  MeshInfo() {};

  /// Adds edge in ascending index order and performs checks
  void addEdge( uint_t idx0, uint_t idx1, DoFType dofType );

  /// Vertices: ( ID, coordinate )
  std::map< uint_t, Point3D > vertices_;

  /// Edges: ( ( VertexID_0, VertexID_1 ), DoFType )
  std::set< std::pair< std::pair< uint_t, uint_t >, DoFType > > edges_;

  /// Faces: ( VertexID_0 )
  std::set< std::tuple< uint_t, uint_t, uint_t > > faces_;

};


class SetupPrimitiveStorage
{
public:

  SetupPrimitiveStorage( const MeshInfo & meshInfo )
  {

  }

private:

  void generatePrimitiveID();

  std::map< PrimitiveID::IDType, SetupVertex* > vertices_;
  std::map< PrimitiveID::IDType, SetupEdge*  >  edges_;
  std::map< PrimitiveID::IDType, SetupFace* >   faces_;

};

} // namespace hhg
