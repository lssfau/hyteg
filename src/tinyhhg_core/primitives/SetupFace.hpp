
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
    SetupPrimitive( id ), edgeID0_( edgeID0 ),
    edgeID1_( edgeID1 ), edgeID2_( edgeID2 )
  {}

private:

  PrimitiveID edgeID0_;
  PrimitiveID edgeID1_;
  PrimitiveID edgeID2_;

};

}
