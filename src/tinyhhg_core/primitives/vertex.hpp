#ifndef VERTEX_HPP
#define VERTEX_HPP

#include <fmt/ostream.h>
#include "tinyhhg_core/types/pointnd.hpp"
#include <core/DataTypes.h>
#include <core/Deprecated.h>
#include <tinyhhg_core/types/flags.hpp>
#include <tinyhhg_core/primitives/Primitive.hpp>
#include <tinyhhg_core/mesh/MeshInfo.hpp>

#include <vector>


namespace hhg
{

class Edge;
class Face;

class VertexMemory;

/// \brief  Macro-Vertex primitive
/// \author Daniel Drzisga (drzisga@ma.tum.de)
/// \date   March, 2017
///
/// The Vertex class represents a Macro-Vertex primitve. It saves geometrical and topological information
/// as well as pointers to memory reserved for \ref Function and \ref Operator.
class Vertex : public Primitive
{
public:

  friend class SetupPrimitiveStorage;
  friend class PrimitiveStorage;

  /// Constructs a vertex with given id and coordinates
  /// \param id Id of vertex
  /// \param coords Spatial coordinates of vertex
  Vertex(size_t id, const Point3D& coords);
  Vertex( const PrimitiveID & primitiveID, const Point3D & coordinates );

  /// Adds given edge to \ref edges
  /// \param edge Pointer to edge which will be added
  void addEdge(Edge* edge);

  /// Adds given face to \ref faces
  /// \param face Pointer to face which will be added
  void addFace(Face* face);

  /// Returns the index of \p edge within \ref edges
  /// \param edge Edge
  /// \returns Index of \p edge within \ref edges
  size_t edge_index(const Edge& edge) const;

  /// Returns the index of \p face within \ref faces
  /// \param face Face
  /// \returns Index of \p face within \ref faces
  size_t face_index(const Face& face) const;

  /// Id of vertex
  size_t id;

  /// Processor rank this vertex belongs to
  walberla::uint_t rank;

  /// DoF type of vertex
  DoFType type;

  /// Spatial coordinates of vertex
  Point3D coords;

  /// Pointers to edges adjacent to vertex
  std::vector<Edge*> edges;

  /// Pointers to faces adjacent to vertex
  std::vector<Face*> faces;


  /// Vector containing pointers to memory used by \ref Function and \ref Operator
  /// The std::vector corresponds to the memory id of a \ref Function or an \Operator
  /// This replaces the old data and opr_data vectors
  std::vector<VertexMemory*> memory;


  /// Method overload for string formatting
  friend std::ostream &operator<<(std::ostream &os, const Vertex &vertex);

  const Point3D getCoordinates() const
  {
    return coords;
  }

  /// Returns a pointer to the data that belongs to the passed \ref PrimitiveDataID.
  /// \param index the \ref PrimitiveDataID of the data that should be returned
  template< typename DataType >
  DataType* getData( const PrimitiveDataID< DataType, Vertex > & index ) const
  {
    return genericGetData< DataType >( index );
  }

  virtual void getLowerDimNeighbors ( std::vector< PrimitiveID > & lowerDimNeighbors )  const { lowerDimNeighbors.clear(); }
  virtual void getHigherDimNeighbors( std::vector< PrimitiveID > & higherDimNeighbors ) const { getNeighborEdges( higherDimNeighbors ); }

  virtual const std::vector< PrimitiveID > & getLowerDimNeighbors()  const { WALBERLA_ASSERT_EQUAL( getNumNeighborVertices(), 0 ); return neighborVertices(); }
  virtual const std::vector< PrimitiveID > & getHigherDimNeighbors() const { return neighborEdges(); }

  virtual uint_t getNumLowerDimNeighbors()  const { return 0; }
  virtual uint_t getNumHigherDimNeighbors() const { return getNumNeighborEdges(); }

protected:

  /// Not public in order to guarantee that data is only added through the governing structure.
  /// This ensures valid DataIDs.
  template< typename DataType >
  inline void addData( const PrimitiveDataID< DataType, Vertex > & index,
	               const PrimitiveDataHandling< DataType, Vertex > & dataHandling )
  {
    genericAddData( index, dataHandling, this );
  }

private:

  void addEdge( const PrimitiveID & edgeID ) { neighborEdges_.push_back( edgeID ); }

};


class VertexMemory
{
public:
  const MemoryType type;
  virtual ~VertexMemory() { }

protected:
  VertexMemory(MemoryType t) : type(t) { ; }
};

}

#endif /* VERTEX_HPP */
