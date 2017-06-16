
#pragma once

namespace hhg {

class SetupVertex : public SetupPrimitive
{
public:

  SetupVertex( const PrimitiveID & id, const Point3D& coordinates ) :
    SetupPrimitive( id ), coordinates_( coordinates )
  {}

private:

  Point3D coordinates_;

};

}
