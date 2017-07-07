
#pragma once

#include "tinyhhg_core/primitives/SetupPrimitive.hpp"

namespace hhg {

class SetupFace : public SetupPrimitive
{
public:

  SetupFace( const SetupPrimitiveStorage & storage,
		     const PrimitiveID & id,
             const PrimitiveID & edgeID0,
             const PrimitiveID & edgeID1,
             const PrimitiveID & edgeID2,
             const std::array< int, 3 > edgeOrientation,
             const std::array< Point3D, 3 > coordinates ) :
    SetupPrimitive( storage, id ), edgeOrientation_( edgeOrientation ), coordinates_( coordinates )
  {
	lowerDimNeighbors_.push_back( edgeID0 );
	lowerDimNeighbors_.push_back( edgeID1 );
	lowerDimNeighbors_.push_back( edgeID2 );
  }

  std::array< int, 3 >     getEdgeOrientation() const { return edgeOrientation_; }
  std::array< Point3D, 3 > getCoordinates()     const { return coordinates_; }

  const PrimitiveID & getEdgeID0() const { WALBERLA_ASSERT_EQUAL( lowerDimNeighbors_.size(), 3 ); return lowerDimNeighbors_[0]; }
  const PrimitiveID & getEdgeID1() const { WALBERLA_ASSERT_EQUAL( lowerDimNeighbors_.size(), 3 ); return lowerDimNeighbors_[1]; }
  const PrimitiveID & getEdgeID2() const { WALBERLA_ASSERT_EQUAL( lowerDimNeighbors_.size(), 3 ); return lowerDimNeighbors_[2]; }

private:

  std::array< int, 3 >     edgeOrientation_;
  std::array< Point3D, 3 > coordinates_;

};

}
