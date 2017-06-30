
#pragma once

#include "tinyhhg_core/primitives/SetupPrimitive.hpp"

namespace hhg {

class SetupEdge : public SetupPrimitive
{
public:

  friend class SetupPrimitiveStorage;

  SetupEdge( const PrimitiveID & id,
	     const PrimitiveID & vertexID0,
	     const PrimitiveID & vertexID1,
	     const DoFType & dofType ) :
    SetupPrimitive( id ), vertexIDs_( { vertexID0, vertexID1 } ),
    dofType_( dofType )
  {}

  const PrimitiveID & getVertexID0() const { return vertexIDs_[0]; }
  const PrimitiveID & getVertexID1() const { return vertexIDs_[1]; }

  const DoFType & getDoFType() const { return dofType_; }

  virtual PrimitiveID::const_iterator beginLowerDimNeighbors() const { return vertexIDs_.begin(); }
  virtual PrimitiveID::const_iterator endLowerDimNeighbors()   const { return vertexIDs_.end(); }

  virtual PrimitiveID::const_iterator beginHigherDimNeighbors() const { return faceIDs_.begin(); }
  virtual PrimitiveID::const_iterator endHigherDimNeighbors()   const { return faceIDs_.end(); }

private:

  void addFace( const PrimitiveID & faceID ) { faceIDs_.push_back( faceID ); }

  std::vector< PrimitiveID > vertexIDs_;
  std::vector< PrimitiveID > faceIDs_;

  DoFType dofType_;
};

}
