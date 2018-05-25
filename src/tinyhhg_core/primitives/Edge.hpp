
#pragma once

#include "tinyhhg_core/types/pointnd.hpp"

#include <vector>

#include <tinyhhg_core/types/flags.hpp>
#include <tinyhhg_core/primitives/Primitive.hpp>

#include <core/DataTypes.h>
#include <core/Deprecated.h>

namespace hhg
{

class Vertex;
class Face;

class Edge : public Primitive
{
public:

  friend class SetupPrimitiveStorage;
  friend class PrimitiveStorage;

  Edge( const PrimitiveID & primitiveID,
        const PrimitiveID & vertexID0,
        const PrimitiveID & vertexID1,
        const std::array<Point3D, 2>& coords) :
    Primitive( primitiveID ), coordinates_(coords)
  {
    neighborVertices_.push_back( vertexID0 );
    neighborVertices_.push_back( vertexID1 );

    direction_ = coordinates_[1] - coordinates_[0];
    length_ = direction_.norm();
    tangent_ = direction_ / length_;

    const std::array<walberla::real_t,3> init{{tangent_[1], -tangent_[0], 0.0}};
    normal2D_ = Point3D(init);
  }

  Edge( walberla::mpi::RecvBuffer & recvBuffer ) : Primitive( recvBuffer ) { deserializeSubclass( recvBuffer ); }

  uint_t vertex_index(const PrimitiveID& vertex) const;
  uint_t face_index(const PrimitiveID& face) const;
  uint_t cell_index(const PrimitiveID& cell) const;

  PrimitiveID get_opposite_vertex(const PrimitiveID& vertex) const;
  bool opposite_face_exists(const PrimitiveID& face) const;
  PrimitiveID get_opposite_face(const PrimitiveID& face) const;

  friend std::ostream &operator<<(std::ostream &os, const Edge &edge);

  const PrimitiveID & getVertexID0() const { WALBERLA_ASSERT_EQUAL( getNumNeighborVertices(), 2 ); return neighborVertices_[0]; }
  const PrimitiveID & getVertexID1() const { WALBERLA_ASSERT_EQUAL( getNumNeighborVertices(), 2 ); return neighborVertices_[1]; }

  const std::array< Point3D, 2 > & getCoordinates() const { return coordinates_; }
  const Point3D                  & getDirection()   const { return direction_; }
  const real_t                   & getLength()      const { return length_; }
  const Point3D                  & getTangent()     const { return tangent_; }
  const Point3D                  & get2DNormal()    const { return normal2D_; }

  /// Returns a pointer to the data that belongs to the passed \ref PrimitiveDataID.
  /// \param index the \ref PrimitiveDataID of the data that should be returned
  template< typename DataType >
  DataType* getData( const PrimitiveDataID< DataType, Edge > & index ) const
  {
    return genericGetData< DataType >( index );
  }

  using Primitive::getData;

  virtual void getLowerDimNeighbors ( std::vector< PrimitiveID > & lowerDimNeighbors )  const { getNeighborVertices( lowerDimNeighbors ); }
  virtual void getHigherDimNeighbors( std::vector< PrimitiveID > & higherDimNeighbors ) const { getNeighborFaces( higherDimNeighbors ); }

  virtual const std::vector< PrimitiveID > & getLowerDimNeighbors()  const { return neighborVertices(); }
  virtual const std::vector< PrimitiveID > & getHigherDimNeighbors() const { return neighborFaces(); }

  virtual uint_t getNumLowerDimNeighbors()  const { return getNumNeighborVertices(); }
  virtual uint_t getNumHigherDimNeighbors() const { return getNumNeighborFaces(); }

protected:

  /// Not public in order to guarantee that data is only added through the governing structure.
  /// This ensures valid DataIDs.
  template< typename DataType, typename DataHandlingType >
  inline void addData( const PrimitiveDataID< DataType, Edge > & index,
                       const std::shared_ptr< DataHandlingType > & dataHandling )
  {
    genericAddData( index, dataHandling, this );
  }

  virtual void   serializeSubclass ( walberla::mpi::SendBuffer & sendBuffer ) const;
  virtual void deserializeSubclass ( walberla::mpi::RecvBuffer & recvBuffer );

private:

  void addFace( const PrimitiveID & faceID ) { neighborFaces_.push_back( faceID ); }
  void addCell( const PrimitiveID & cellID ) { neighborCells_.push_back( cellID ); }

  std::array<Point3D, 2> coordinates_;
  Point3D direction_;
  real_t  length_;
  Point3D tangent_;
  Point3D normal2D_;

};

}

