#pragma once

#include "tinyhhg_core/types/matrix.hpp"
#include "tinyhhg_core/types/pointnd.hpp"

namespace hhg {

class GeometryMap
{
 public:
   enum class Type : uint_t
   {
      IDENTITY = 0,
      CIRCULAR = 1
   };

   virtual ~GeometryMap(){};

   virtual void evalF( const Point3D&, Point3D& ) const      = 0;
   virtual void evalDF( const Point3D&, Matrix2r& ) const    = 0;
   virtual void evalDFinv( const Point3D&, Matrix2r& ) const = 0;

   virtual void                          serialize( walberla::mpi::SendBuffer& sendBuffer ) = 0;
   static std::shared_ptr< GeometryMap > deserialize( walberla::mpi::RecvBuffer& recvBuffer );
};

} // namespace hhg
