
#pragma once

#include "tinyhhg_core/primitives/SetupPrimitive.hpp"

namespace hhg {

class SetupEdge : public SetupPrimitive
{
public:

  SetupEdge( const PrimitiveID & id,
	     const PrimitiveID & vertexID0,
	     const PrimitiveID & vertexID1,
	     const DoFType & dofType ) :
    SetupPrimitive( id ), vertexID0_( vertexID0 ), vertexID1_( vertexID1 ),
    dofType_( dofType )
  {}

private:

  PrimitiveID vertexID0_;
  PrimitiveID vertexID1_;

  DoFType dofType_;
};

}
