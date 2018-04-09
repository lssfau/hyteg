#include "CircularMap.hpp"
#include "IdentityMap.hpp"

namespace hhg {

// TODO: remove me before merging into master
void GeometryMap::evalTensorCoeff(const Point3D& x, Matrix2r& coeff) {
  Point3D xPhy;
  evalF(x, xPhy);
  Matrix2r Khat;
  // TODO: remove hardcoded TENSOR
  Khat(0,0) = 1.0;
  Khat(0,1) = 0.0;
  Khat(1,0) = Khat(0,1);
  Khat(1,1) = 1.0;
  evalDFinv(x, coeff);
  real_t invDet = coeff.det();
  Matrix2r coeffT = coeff.transpose();
  coeff = coeff.mul(Khat.mul(coeffT));
  coeff *= 1.0 / std::abs(invDet);
}


real_t GeometryMap::evalDetDF( const Point3D& x)  {
   Matrix2r DF;
   evalDF( x, DF );
   return DF.det();
}


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
