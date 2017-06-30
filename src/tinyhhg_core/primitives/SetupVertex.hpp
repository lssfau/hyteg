
#pragma once

#include "tinyhhg_core/primitives/SetupPrimitive.hpp"

namespace hhg {

class SetupVertex : public SetupPrimitive
{
public:

  SetupVertex( const PrimitiveID & id, const Point3D& coordinates ) :
    SetupPrimitive( id ), coordinates_( coordinates )
  {}

  const Point3D & getCoordinates() const
  {
    return coordinates_;
  }

  virtual PrimitiveID::const_iterator beginLowerDimNeighbors() const { return empty_.begin(); }
  virtual PrimitiveID::const_iterator endLowerDimNeighbors()   const { return empty_.end(); }

  virtual PrimitiveID::const_iterator beginHigherDimNeighbors() const { return edgeIDs_.begin(); }
  virtual PrimitiveID::const_iterator endHigherDimNeighbors()   const { return edgeIDs_.end(); }

private:

  Point3D coordinates_;

  std::vector< PrimitiveID > edgeIDs_;
  const std::vector< PrimitiveID > empty_;

};

}
