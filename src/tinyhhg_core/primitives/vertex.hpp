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
  Vertex( const PrimitiveID & primitiveID, const Point3D & coordinates );

  Vertex( walberla::mpi::RecvBuffer & recvBuffer ) : Primitive( recvBuffer ) { deserializeSubclass( recvBuffer ); }

  /// Returns the index of \p edge within \ref edges
  /// \param edge Edge
  /// \returns Index of \p edge within \ref edges
  uint_t edge_index(const PrimitiveID& edge) const;

  /// Returns the index of \p face within \ref faces
  /// \param face Face
  /// \returns Index of \p face within \ref faces
  uint_t face_index(const PrimitiveID& face) const;

  /// DoF type of vertex
  DoFType getDoFType() const { return dofType_; }

  /// Method overload for string formatting
  friend std::ostream &operator<<(std::ostream &os, const Vertex &vertex);

  /// Spatial coordinates of vertex
  const Point3D getCoordinates() const { return coordinates_; }

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
  template< typename DataType, typename DataHandlingType >
  inline void addData( const PrimitiveDataID< DataType, Vertex > & index,
                       const std::shared_ptr< DataHandlingType > & dataHandling )
  {
    genericAddData( index, dataHandling, this );
  }

  virtual void   serializeSubclass ( walberla::mpi::SendBuffer & sendBuffer ) const;
  virtual void deserializeSubclass ( walberla::mpi::RecvBuffer & recvBuffer );

private:

  void addEdge( const PrimitiveID & edgeID ) { neighborEdges_.push_back( edgeID ); }
  void addFace( const PrimitiveID & faceID ) { neighborFaces_.push_back( faceID ); }

  DoFType dofType_;
  Point3D coordinates_;
};

}

#endif /* VERTEX_HPP */
