
#pragma once

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

  Face( const PrimitiveID & primitiveID,
        const std::array<PrimitiveID, 3>& vertexIDs,
        const std::array<PrimitiveID, 3>& edgeIDs,
        const std::array< int, 3 >     & edgeOrientation,
        const std::array< Point3D, 3 > & coordinates ) :
    Primitive( primitiveID ), type( Inner ), edge_orientation( edgeOrientation ), coords( coordinates )
  {
    neighborVertices_.push_back( vertexIDs[0] );
    neighborVertices_.push_back( vertexIDs[1] );
    neighborVertices_.push_back( vertexIDs[2] );

    neighborEdges_.push_back( edgeIDs[0] );
    neighborEdges_.push_back( edgeIDs[1] );
    neighborEdges_.push_back( edgeIDs[2] );
    std::array<Point3D, 2> B({{coords[1]-coords[0], coords[2] - coords[0]}});
    area = std::abs(0.5 * math::det2(B));
  }

  Face( walberla::mpi::RecvBuffer & recvBuffer ) : Primitive( recvBuffer ) { deserializeSubclass( recvBuffer ); }

  uint_t vertex_index(const PrimitiveID& vertex) const;
  uint_t edge_index(const PrimitiveID& edge) const;

  std::vector<PrimitiveID> adjacent_edges(const PrimitiveID& vertex) const;
  PrimitiveID get_vertex_opposite_to_edge(const PrimitiveID& edge) const;

  DoFType type;
  real_t area;
  std::array<int, 3> edge_orientation;
  std::array<Point3D, 3> coords;

  friend std::ostream &operator<<(std::ostream &os, const Face &face);

  std::array< int, 3 >     getEdgeOrientation() const { return edge_orientation; }
  std::array< Point3D, 3 > getCoordinates()     const { return coords; }

  const PrimitiveID & getVertexID0() const { WALBERLA_ASSERT_EQUAL( getNumNeighborVertices(), 3 ); return neighborVertices_[0]; }
  const PrimitiveID & getVertexID1() const { WALBERLA_ASSERT_EQUAL( getNumNeighborVertices(), 3 ); return neighborVertices_[1]; }
  const PrimitiveID & getVertexID2() const { WALBERLA_ASSERT_EQUAL( getNumNeighborVertices(), 3 ); return neighborVertices_[2]; }

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

  using Primitive::getData;

  virtual void getLowerDimNeighbors ( std::vector< PrimitiveID > & lowerDimNeighbors )  const { getNeighborEdges( lowerDimNeighbors ); }
  virtual void getHigherDimNeighbors( std::vector< PrimitiveID > & higherDimNeighbors ) const { higherDimNeighbors.clear(); }

  virtual const std::vector< PrimitiveID > & getLowerDimNeighbors()  const { return neighborEdges(); }
  virtual const std::vector< PrimitiveID > & getHigherDimNeighbors() const { WALBERLA_ASSERT_EQUAL( getNumNeighborFaces(), 0 ); return neighborFaces(); }

  virtual uint_t getNumLowerDimNeighbors()  const { return getNumNeighborEdges(); }
  virtual uint_t getNumHigherDimNeighbors() const { return 0; }

protected:

  /// Not public in order to guarantee that data is only added through the governing structure.
  /// This ensures valid DataIDs.
  template< typename DataType, typename DataHandlingType >
  inline void addData( const PrimitiveDataID< DataType, Face > & index,
                       const std::shared_ptr< DataHandlingType > & dataHandling )
  {
    genericAddData( index, dataHandling, this );
  }

  virtual void   serializeSubclass ( walberla::mpi::SendBuffer & sendBuffer ) const;
  virtual void deserializeSubclass ( walberla::mpi::RecvBuffer & recvBuffer );

};

}

