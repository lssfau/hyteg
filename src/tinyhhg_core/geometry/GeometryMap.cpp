#include "CircularMap.hpp"
#include "IdentityMap.hpp"

namespace hhg {

void GeometryMap::serialize( const std::shared_ptr<GeometryMap>& map, walberla::mpi::SendBuffer& sendBuffer ) {
   map->serializeSubClass(sendBuffer);
}

std::shared_ptr< GeometryMap > GeometryMap::deserialize( walberla::mpi::RecvBuffer& recvBuffer )
{
   Type type;
   recvBuffer >> type;

   switch( type )
   {
   case Type::IDENTITY:
      return std::make_shared< IdentityMap >();
   case Type::CIRCULAR:
      return std::make_shared< CircularMap >( recvBuffer );
   default:
      WALBERLA_ABORT( "Error in deserializing GeometryMap: Unknown Type" )
   }
}

} // namespace hhg
