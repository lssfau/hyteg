
#pragma once

#include "tinyhhg_core/primitives/Primitive.hpp"

namespace hhg {

class Cell : public Primitive
{
public:

  Cell( const PrimitiveID                & primitiveID,
        const std::vector< PrimitiveID > & vertexIDs,
        const std::vector< PrimitiveID > & edgeIDs,
        const std::vector< PrimitiveID > & faceIDs );

  virtual void getLowerDimNeighbors ( std::vector< PrimitiveID > & lowerDimNeighbors )  const { getNeighborFaces( lowerDimNeighbors ); }
  virtual void getHigherDimNeighbors( std::vector< PrimitiveID > & higherDimNeighbors ) const { higherDimNeighbors.clear(); }

  virtual const std::vector< PrimitiveID > & getLowerDimNeighbors()  const { return neighborFaces(); }
  virtual const std::vector< PrimitiveID > & getHigherDimNeighbors() const { WALBERLA_ASSERT_EQUAL( getNumNeighborCells(), 0 ); return neighborCells(); }

  virtual uint_t getNumLowerDimNeighbors()  const { return getNumNeighborFaces(); }
  virtual uint_t getNumHigherDimNeighbors() const { return 0; }

protected:

  virtual void   serializeSubclass ( walberla::mpi::SendBuffer & sendBuffer ) const { WALBERLA_UNUSED( sendBuffer ); }
  virtual void deserializeSubclass ( walberla::mpi::RecvBuffer & recvBuffer )       { WALBERLA_UNUSED( recvBuffer ); };

};

}
