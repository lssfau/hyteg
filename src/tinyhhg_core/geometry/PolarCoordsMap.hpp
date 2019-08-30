#pragma once

#include <cmath>
#include "GeometryMap.hpp"

namespace hyteg {

  /// Class providing geometry mapping based on polar coordinates; convention is x[0] = $\rho$, x[1] = $\varphi$
  class PolarCoordsMap : public GeometryMap
  {
  public:

    PolarCoordsMap(){}

    void evalF( const Point3D& x, Point3D& Fx ) const {
      Fx[0] = x[0] * std::cos( x[1] );
      Fx[1] = x[0] * std::sin( x[1] );
    }

    void evalDF( const Point3D& x, Matrix2r& DFx ) const
    {
      DFx( 0, 0 ) =          std::cos( x[1] );
      DFx( 0, 1 ) = - x[0] * std::sin( x[1] );
      DFx( 1, 0 ) =          std::sin( x[1] );
      DFx( 1, 1 ) =   x[0] * std::cos( x[1] );
    }

    void evalDFinv( const Point3D& x, Matrix2r& DFinvx ) const
    {
      DFinvx( 0, 0 ) =   std::cos( x[1] );
      DFinvx( 0, 1 ) =   std::sin( x[1] );
      DFinvx( 1, 0 ) = - std::sin( x[1] ) / x[0];
      DFinvx( 1, 1 ) =   std::cos( x[1] ) / x[0];
    }

    void serializeSubClass( walberla::mpi::SendBuffer& sendBuffer ) const {
      sendBuffer << Type::POLAR_COORDS;
    }

  };

} // namespace hyteg
