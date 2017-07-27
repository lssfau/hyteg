#ifndef FACE_HPP
#define FACE_HPP

#include <tinyhhg_core/types/pointnd.hpp>
#include <tinyhhg_core/types/flags.hpp>
#include <tinyhhg_core/primitives/Primitive.hpp>
#include <tinyhhg_core/primitivestorage/PrimitiveStorage.hpp>
#include <tinyhhg_core/math.hpp>
#include <core/DataTypes.h>
#include <core/Deprecated.h>

#include <array>
#include <vector>

namespace hhg
{

class Vertex;
class Edge;

class Face : public Primitive
{
public:

  friend class SetupPrimitiveStorage;
  friend class PrimitiveStorage;

  Face(size_t id, Edge* edges[3]);

  Face( const PrimitiveID & primitiveID,
        const PrimitiveID & edgeID0,
        const PrimitiveID & edgeID1,
        const PrimitiveID & edgeID2,
        const std::array< int, 3 >     & edgeOrientation,
        const std::array< Point3D, 3 > & coordinates ) :
    Primitive( primitiveID ), id( primitiveID.getID() ), rank( 0 ), type( Inner ), edge_orientation( edgeOrientation ), coords( coordinates )
  {
    neighborEdges_.push_back( edgeID0 );
    neighborEdges_.push_back( edgeID1 );
    neighborEdges_.push_back( edgeID2 );
    std::array<Point3D, 2> B({{coords[1]-coords[0], coords[2] - coords[0]}});
    area = std::abs(0.5 * math::det2(B));
  }

  size_t vertex_index(const Vertex& vertex) const;
  size_t edge_index(const Edge& edge) const;

  std::vector<Edge*> adjacent_edges(const Vertex& vertex) const;
  Vertex* get_vertex_opposite_to_edge(const Edge& edge) const;

  size_t id;
  walberla::uint_t rank;
  DoFType type;
  real_t area;
  Edge* edges[3];
  std::vector<Vertex*> vertices;
  std::array<int, 3> edge_orientation;
  std::array<Point3D, 3> coords;

  std::vector<void*> memory;

  friend std::ostream &operator<<(std::ostream &os, const Face &face);

  std::array< int, 3 >     getEdgeOrientation() const { return edge_orientation; }
  std::array< Point3D, 3 > getCoordinates()     const { return coords; }

  const PrimitiveID & getEdgeID0() const { WALBERLA_ASSERT_EQUAL( getNumNeighborEdges(), 3 ); return neighborEdges_[0]; }
  const PrimitiveID & getEdgeID1() const { WALBERLA_ASSERT_EQUAL( getNumNeighborEdges(), 3 ); return neighborEdges_[1]; }
  const PrimitiveID & getEdgeID2() const { WALBERLA_ASSERT_EQUAL( getNumNeighborEdges(), 3 ); return neighborEdges_[2]; }

  /// Returns a pointer to the data that belongs to the passed \ref PrimitiveDataID.
  /// \param index the \ref PrimitiveDataID of the data that should be returned
  template< typename DataType >
  DataType* getData( const PrimitiveDataID< DataType, Face > & index ) const
  {
    return genericGetData< DataType >( index );
  }

  virtual void getLowerDimNeighbors ( std::vector< PrimitiveID > & lowerDimNeighbors )  const { getNeighborEdges( lowerDimNeighbors ); }
  virtual void getHigherDimNeighbors( std::vector< PrimitiveID > & higherDimNeighbors ) const { higherDimNeighbors.clear(); }

  virtual const std::vector< PrimitiveID > & getLowerDimNeighbors()  const { return neighborEdges(); }
  virtual const std::vector< PrimitiveID > & getHigherDimNeighbors() const { WALBERLA_ASSERT_EQUAL( getNumNeighborFaces(), 0 ); return neighborFaces(); }

  virtual uint_t getNumLowerDimNeighbors()  const { return getNumNeighborEdges(); }
  virtual uint_t getNumHigherDimNeighbors() const { return 0; }

protected:

  /// Not public in order to guarantee that data is only added through the governing structure.
  /// This ensures valid DataIDs.
  template< typename DataType >
  inline void addData( const PrimitiveDataID< DataType, Face > & index,
		       const PrimitiveDataHandling< DataType, Face > & dataHandling )
  {
    genericAddData( index, dataHandling, this );
  }

};

}

#endif /* FACE_HPP */
