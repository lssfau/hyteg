
#pragma once

#include "tinyhhg_core/primitives/SetupPrimitive.hpp"

namespace hhg {

class SetupVertex : public SetupPrimitive
{
public:

  friend class SetupPrimitiveStorage;

  SetupVertex( const PrimitiveID & id, const Point3D& coordinates ) :
    SetupPrimitive( id ), coordinates_( coordinates )
  {}

  const Point3D & getCoordinates() const
  {
    return coordinates_;
  }

private:

  void addEdge( const PrimitiveID & edgeID ) { higherDimNeighbors_.push_back( edgeID ); }

  Point3D coordinates_;

};

}
