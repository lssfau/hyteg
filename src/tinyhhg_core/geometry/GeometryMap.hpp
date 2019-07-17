#pragma once

#include "tinyhhg_core/types/matrix.hpp"
#include "tinyhhg_core/types/pointnd.hpp"

namespace hhg {

/// Class describing a geometrical mapping from reference to physical space
class GeometryMap
{
 public:
   enum class Type : uint_t
   {
      IDENTITY = 0,
      AFFINE   = 1,
      CIRCULAR = 2,
      POLAR_COORDS = 3
   };

   virtual ~GeometryMap(){};

   /// Mapping of reference coordinates \p x to physical coordinates \p Fx
   /// \param x Reference input coordinates
   /// \param Fx Physical output coordinates
   virtual void evalF( const Point3D& x, Point3D& Fx ) const          = 0;

   /// Evaluation of the Jacobian matrix at reference position \p x
   /// \param x Reference input coordinates
   /// \param DFx Jacobian matrix
   virtual void evalDF( const Point3D& x, Matrix2r& DFx ) const       = 0;

   /// Evaluation of the Jacobian matrix at reference position \p x
   /// \param x Reference input coordinates
   /// \param DFinvx Inverse of the Jacobian matrix
   virtual void evalDFinv( const Point3D& x, Matrix2r& DFinvx ) const = 0;

   /// Evaluation of the determinant of the Jacobian matrix at reference position \p x
   /// \param x Reference input coordinates
   real_t evalDetDF( const Point3D& x );

   /// Serialization of a GeometryMap object
   static void serialize( const std::shared_ptr< GeometryMap >& map, walberla::mpi::SendBuffer& sendBuffer );

   /// Deserialization of a GeometryMap object
   static std::shared_ptr< GeometryMap > deserialize( walberla::mpi::RecvBuffer& recvBuffer );

 protected:
   virtual void serializeSubClass( walberla::mpi::SendBuffer& sendBuffer ) const = 0;
};

} // namespace hhg
