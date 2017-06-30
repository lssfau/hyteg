
#pragma once

#include "tinyhhg_core/primitives/SetupPrimitive.hpp"

namespace hhg {

class SetupFace : public SetupPrimitive
{
public:

  SetupFace( const PrimitiveID & id,
	     const PrimitiveID & edgeID0,
	     const PrimitiveID & edgeID1,
	     const PrimitiveID & edgeID2 ) :
    SetupPrimitive( id ), edgeIDs_( { edgeID0, edgeID1, edgeID2 } )
  {}

  const PrimitiveID & getEdgeID0() const { return edgeIDs_[0]; }
  const PrimitiveID & getEdgeID1() const { return edgeIDs_[1]; }
  const PrimitiveID & getEdgeID2() const { return edgeIDs_[2]; }

  virtual PrimitiveID::const_iterator beginLowerDimNeighbors() const { return edgeIDs_.begin(); }
  virtual PrimitiveID::const_iterator endLowerDimNeighbors()   const { return edgeIDs_.end(); }

  virtual PrimitiveID::const_iterator beginHigherDimNeighbors() const { return volumeIDs_.begin(); }
  virtual PrimitiveID::const_iterator endHigherDimNeighbors()   const { return volumeIDs_.end(); }

private:

  std::vector< PrimitiveID > edgeIDs_;
  std::vector< PrimitiveID > volumeIDs_;

};

}
