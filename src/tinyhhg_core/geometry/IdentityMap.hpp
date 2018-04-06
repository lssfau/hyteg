#pragma once

#include "GeometryMap.hpp"

namespace hhg {

class IdentityMap : public GeometryMap
{
 public:
   IdentityMap() {}

   void evalF( const Point3D& x, Point3D& Fx ) const { Fx = x; }

   void evalDF( const Point3D&, Matrix2r& DFx ) const
   {
      DFx( 0, 0 ) = 1.0;
      DFx( 0, 1 ) = 0.0;
      DFx( 1, 0 ) = 0.0;
      DFx( 1, 1 ) = 1.0;
   }

   void evalDFinv( const Point3D&, Matrix2r& DFinvx ) const
   {
      DFinvx( 0, 0 ) = 1.0;
      DFinvx( 0, 1 ) = 0.0;
      DFinvx( 1, 0 ) = 0.0;
      DFinvx( 1, 1 ) = 1.0;
   }

   void serialize( walberla::mpi::SendBuffer& sendBuffer ) { sendBuffer << Type::IDENTITY; }
};

} // namespace hhg
