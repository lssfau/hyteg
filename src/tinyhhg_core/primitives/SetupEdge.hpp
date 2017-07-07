
#pragma once

#include "tinyhhg_core/primitives/SetupPrimitive.hpp"

namespace hhg {

class SetupEdge : public SetupPrimitive
{
public:

  friend class SetupPrimitiveStorage;

  SetupEdge( const SetupPrimitiveStorage & storage,
		     const PrimitiveID & id,
	         const PrimitiveID & vertexID0,
	         const PrimitiveID & vertexID1,
	         const DoFType & dofType,
	         const Point3D & direction ) :
    SetupPrimitive( storage, id ), dofType_( dofType ), direction_( direction )
  {
	lowerDimNeighbors_.push_back( vertexID0 );
	lowerDimNeighbors_.push_back( vertexID1 );
  }

  const PrimitiveID & getVertexID0() const { WALBERLA_ASSERT_EQUAL( lowerDimNeighbors_.size(), 2 ); return lowerDimNeighbors_[0]; }
  const PrimitiveID & getVertexID1() const { WALBERLA_ASSERT_EQUAL( lowerDimNeighbors_.size(), 2 ); return lowerDimNeighbors_[1]; }

  const DoFType & getDoFType() const   { return dofType_; }
  const Point3D & getDirection() const { return direction_; }

private:

  void addFace( const PrimitiveID & faceID ) { higherDimNeighbors_.push_back( faceID ); }

  DoFType dofType_;
  Point3D direction_;
};

}
