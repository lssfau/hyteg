
#pragma once

#include "core/debug/Debug.h"
#include "tinyhhg_core/types/pointnd.hpp"
#include "tinyhhg_core/types/flags.hpp"

#include <array>
#include <map>
#include <set>
#include <vector>

namespace hhg {

using walberla::uint_t;
using walberla::real_t;

/// \brief Contains information about a mesh
/// \author Daniel Drzisga (drzisga@ma.tum.de) \n
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

  /// Vertices: ( ID, coordinate )
  typedef std::map< uint_t, Point3D >                       VertexContainer;
  /// Edges: ( ( VertexID_0, VertexID_1 ), DoFType )
  typedef std::map< std::pair< uint_t, uint_t >, DoFType >  EdgeContainer;
  /// Faces: ( ( VertexID_0, VertexID_1, VertexID_2 )
  typedef std::set< std::array< uint_t, 3 > >               FaceContainer;

  /// Construct empty MeshInfo (for testing)
  static MeshInfo emptyMeshInfo() { return MeshInfo(); }
  /// Construct a MeshInfo from a file in Gmsh format
  static MeshInfo fromGmshFile( const std::string & meshFileName );

  /// Returns vertices of the mesh
  /// \return Map that maps ( vertex ID -> vertex coordinate )
  const VertexContainer & getVertices() const { return vertices_; };

  /// Returns edges of the mesh
  /// \return Map that maps \n
  ///         - a pair of vertex indices (matching indices from the vertices from \ref getVertices \n
  ///         - to the corresponding \ref DoFType
  const EdgeContainer & getEdges() const { return edges_; };

  /// Returns faces of the mesh
  /// \return Set of 3-tuples of vertex indices from the vertices from \ref getVertices
  const FaceContainer & getFaces() const { return faces_; }

private:

  MeshInfo() {};

  /// Adds edge in ascending index order and performs checks
  void addEdge( uint_t idx0, uint_t idx1, DoFType dofType );

  VertexContainer vertices_;
  EdgeContainer edges_;
  FaceContainer faces_;

};

}
