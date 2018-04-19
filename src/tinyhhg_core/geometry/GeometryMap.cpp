#include "AffineMap.hpp"
#include "CircularMap.hpp"
#include "IdentityMap.hpp"

namespace hhg {

real_t GeometryMap::evalDetDF( const Point3D& x )
{
   Matrix2r DF;
   evalDF( x, DF );
   return DF.det();
}

void GeometryMap::serialize( const std::shared_ptr< GeometryMap >& map, walberla::mpi::SendBuffer& sendBuffer )
{
   map->serializeSubClass( sendBuffer );
}

std::shared_ptr< GeometryMap > GeometryMap::deserialize( walberla::mpi::RecvBuffer& recvBuffer )
{
   Type type;
   recvBuffer >> type;

   switch( type )
   {
   case Type::IDENTITY:
      return std::make_shared< IdentityMap >();
   case Type::AFFINE:
      return std::make_shared< AffineMap >( recvBuffer );
   case Type::CIRCULAR:
      return std::make_shared< CircularMap >( recvBuffer );
   default:
      WALBERLA_ABORT( "Error in deserializing GeometryMap: Unknown Type" )
   }
}

} // namespace hhg
