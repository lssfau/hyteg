
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

  typedef uint_t IDType;

  class Vertex
  {
  public:
    Vertex() : id_( 0 ), coordinates_( Point3D() ), dofType_( Inner ) {};
    Vertex( const IDType & id, const Point3D & coordinates, const DoFType & dofType ) :
      id_( id ), coordinates_( coordinates ), dofType_( dofType )
    {}

    const IDType  getID()          const { return id_; }
    const Point3D getCoordinates() const { return coordinates_; }
    const DoFType getDoFType()     const { return dofType_; }

  private:
    IDType  id_;
    Point3D coordinates_;
    DoFType dofType_;
  };

  class Edge
  {
  public:
    Edge() : dofType_( Inner ) {};
    Edge( const std::array< IDType, 2 > & vertices, const DoFType & dofType ) :
      vertices_( vertices ), dofType_( dofType )
    {}

    const std::array< IDType, 2 > getVertices() const { return vertices_; }
    const DoFType                 getDoFType()  const { return dofType_; }

  private:
    std::array< IDType, 2 > vertices_;
    DoFType                 dofType_;
  };

  class Face
  {
  public:
    Face() : dofType_( Inner ) {};
    Face( const std::vector< IDType > & vertices, const DoFType & dofType ) :
      vertices_( vertices ), dofType_( dofType )
    {}

    const std::vector< IDType > getVertices() const { return vertices_; }
    const DoFType               getDoFType()  const { return dofType_; }

  private:
    std::vector< IDType > vertices_;
    DoFType               dofType_;
  };

  class Cell
  {
  public:
    Cell() : dofType_( Inner ) {};
    Cell( const std::vector< IDType > & vertices, const DoFType & dofType ) :
      vertices_( vertices ), dofType_( dofType )
    {}

    const std::vector< IDType > getVertices() const { return vertices_; }
    const DoFType               getDoFType()  const { return dofType_; }

  private:
    std::vector< IDType > vertices_;
    DoFType               dofType_;
  };


  typedef std::map< IDType,                  Vertex > VertexContainer;
  typedef std::map< std::array< IDType, 2 >, Edge   > EdgeContainer;
  typedef std::map< std::vector< IDType >,   Face   > FaceContainer;
  typedef std::map< std::vector< IDType >,   Cell   > CellContainer;

  /// Construct empty MeshInfo (for testing)
  static MeshInfo emptyMeshInfo() { return MeshInfo(); }
  /// Construct a MeshInfo from a file in Gmsh format
  static MeshInfo fromGmshFile( const std::string & meshFileName );

  /// Returns vertices of the mesh
  const VertexContainer & getVertices() const { return vertices_; };

  /// Returns edges of the mesh
  const EdgeContainer & getEdges() const { return edges_; };

  /// Returns faces of the mesh
  const FaceContainer & getFaces() const { return faces_; }

  /// Returns cells of the mesh
  const CellContainer & getCells() const { return cells_; }

private:

  MeshInfo() {};

  /// Adds edge in ascending index order and performs checks
  void addEdge( const Edge & edge );

  /// Adds face in ascending index order and performs checks
  void addFace( const Face & face );

  VertexContainer vertices_;
  EdgeContainer   edges_;
  FaceContainer   faces_;
  CellContainer   cells_;

};

}
