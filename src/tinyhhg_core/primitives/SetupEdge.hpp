
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

  const PrimitiveID & getVertexID0() const { return vertexID0_; }
  const PrimitiveID & getVertexID1() const { return vertexID1_; }

  const DoFType & getDoFType() const { return dofType_; }

private:

  PrimitiveID vertexID0_;
  PrimitiveID vertexID1_;

  DoFType dofType_;
};

}
