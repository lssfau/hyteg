/*
 * Copyright (c) 2017-2024 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Andreas Burkhart.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include "core/Abort.h"

#include "hyteg/types/Matrix.hpp"
#include "hyteg/types/PointND.hpp"

using walberla::real_c;

namespace hyteg {

/// Class describing a geometrical mapping from reference to physical space
///
/// An object of type GeometryMap describes a blending map, i.e. a local
/// coordinate transformation, given by
/// \f[
/// \Psi(\xi_0,\xi_1,\ldots,\xi_{d-1}) = \left(\begin{array}{c}
/// x_0(\xi_0,\xi_1,\ldots,\xi_{d-1})\\
/// x_1(\xi_0,\xi_1,\ldots,\xi_{d-1})\\
/// \vdots\\
/// x_{d-1}(\xi_0,\xi_1,\ldots,\xi_{d-1})
/// \end{array}\right)
/// \f]
/// Blending maps are HyTeG's approach to work with domains whose boundary
/// is not polyhedral. In this case we distinguish between a computational
/// domain \f$T\f$ and a physical domain \f$\hat{T}\f$. The latter two are
/// connected via a blending map:
///
/// <center>
/// <img src="BlendingMaps.png" width="600" height="600"/>
/// </center>
///
/// In order to allow integral transformations to the reference element
/// \f$T_0\f$ the requirements on a blending map are:
/// - The local mapping \f$\Psi\Bigm|_{M} : M \rightarrow \hat{T}\f$ must be
///   a diffeomorphism for every macro element \f$M\in T\f$.
/// - The global mapping \f$\Psi : T \rightarrow \hat{T}\f$ must be a
///   homeomorphism.
///
/// Note that this allows to construct \f$\Psi\f$ in a piecewise fashion.
///
/// If you want to use the method evalDFinvDF(), then, obviously, higher
/// local smoothness requirements apply.
class GeometryMap
{
 public:
   /// These enum values are use in the serialisation of blending maps to mark
   /// the type of map contained in the stream to allow correct deserialisation
   enum class Type : uint_t
   {
      IDENTITY          = 0,
      AFFINE            = 1, // [keep for backward compatibility or remove?]
      CIRCULAR          = 2,
      POLAR_COORDS      = 3,
      ANNULUS           = 4,
      ICOSAHEDRAL_SHELL = 5,
      AFFINE_2D         = 6,
      AFFINE_3D         = 7,
      THIN_SHELL        = 8,
      TOKAMAK           = 9,
      SPHERICAL_COORDS  = 10,
      TORUS             = 11,
   };

   virtual ~GeometryMap(){};

   /// Mapping of reference coordinates \p x to physical coordinates \p Fx
   /// \param x Reference input coordinates
   /// \param Fx Physical output coordinates
   virtual void evalF( const Point3D& x, Point3D& Fx ) const = 0;

   /// Maps point from physical back to computational domain (inverse blending)
   /// \param xPhys coordinates of point in physical domain
   /// \param xComp coordinates of point in computational domain
   virtual void evalFinv( const Point3D& xPhys, Point3D& xComp ) const
   {
      WALBERLA_UNUSED( xPhys );
      WALBERLA_UNUSED( xComp );
      WALBERLA_ABORT( "GeometryMap::evalFinv() not implemented for current child map" );
   }

   /// Evaluation of the Jacobian matrix at reference position \p x
   /// \param x Reference input coordinates
   /// \param DFx Jacobian matrix
   virtual void evalDF( const Point3D& x, Matrix2r& DFx ) const = 0;

   /// Evaluation of the Jacobian matrix at reference position \p x
   /// \param x Reference input coordinates
   /// \param DFx Jacobian matrix
   /// \return value of Jacobian determinant
   virtual real_t evalDF( const Point3D& x, Matrix3r& DFx ) const
   {
      WALBERLA_UNUSED( x );
      WALBERLA_UNUSED( DFx );
      WALBERLA_ABORT( "GeometryMap::evalDF() not implemented for 3D in child class!" );
   };

   ///@{
   /// \brief Evaluation of the derivatives of the inverse Jacobian matrix at reference position \p x
   ///
   /// Denote by
   /// \f[
   /// J_{\Psi} =  \left(\begin{array}{ccc}
   /// \displaystyle\frac{\partial x_0}{\partial \xi_0} & \cdots &
   /// \displaystyle\frac{\partial x_0}{\partial \xi_{d-1}} \\
   /// \vdots & \ddots & \vdots \\
   /// \displaystyle\frac{\partial x_{d-1}}{\partial \xi_0} & \cdots &
   /// \displaystyle\frac{\partial x_{d-1}}{\partial \xi_{d-1}}
   /// \end{array}\right)
   /// \f]
   /// the Jacobian matrix of the blending map. The method evalDFinvDF() computes the derivative
   /// of the inverse of this Jacobian. The returned array DFinvDFx is of dimension \f$(d \times d^2)\f$,
   /// where \f$d\f$ is the spatial dimension, and has entries given by
   ///
   /// \f[
   /// \left(\text{DFinvDFx}\right)_{i,k} = \frac{\partial \left(J^{-1}_{\Psi}\right)_{i,k\bmod d}}{\partial \xi_j}
   /// \f]
   ///
   /// for \f$i\in\{0,\ldots,d-1\}, k\in\{0,\ldots,d^2-1\}\f$ with
   /// \f$j = \displaystyle\left\lfloor \frac{k}{d} \right\rfloor\f$.
   ///
   /// \param x Reference input coordinates
   /// \param DFinvDFx Result matrix of size dim x dim*dim
   virtual void evalDFinvDF( const Point3D& x, Matrixr< 2, 4 >& DFinvDFx ) const
   {
      WALBERLA_UNUSED( x );
      WALBERLA_UNUSED( DFinvDFx );
      WALBERLA_ABORT( "GeometryMap::evalDFinvDF() not implemented for 2D in child class!" );
   };
   virtual void evalDFinvDF( const Point3D& x, Matrixr< 3, 9 >& DFinvDFx ) const
   {
      WALBERLA_UNUSED( x );
      WALBERLA_UNUSED( DFinvDFx );
      WALBERLA_ABORT( "GeometryMap::evalDFinvDF() not implemented for 3D in child class!" );
   };
   ///@}

   /// Evaluation of the Jacobian matrix at reference position \p x
   /// \param x Reference input coordinates
   /// \param DFinvx Inverse of the Jacobian matrix
   virtual void evalDFinv( const Point3D& x, Matrix2r& DFinvx ) const = 0;

   /// Evaluation of the Jacobian matrix at reference position \p x
   /// \param x Reference input coordinates
   /// \param DFinvx Inverse of the Jacobian matrix
   virtual void evalDFinv( const Point3D& x, Matrix3r& DFinvx ) const
   {
      WALBERLA_UNUSED( x );
      WALBERLA_UNUSED( DFinvx );
      WALBERLA_ABORT( "GeometryMap::evalDFinv() not implemented for 3D in child class!" );
   };

   /// Evaluation of the determinant of the Jacobian matrix at reference position \p x
   /// \param x Reference input coordinates
   real_t evalDetDF( const Point3D& x );

   /// Verify that the given point from the physical domain is connected to the given
   /// point on the computational domain.
   ///
   /// We require a blending map to be a homeomorphism globally. However, we allow that global blending
   /// maps are constructed in a piecewise fashion from maps that differ locally between primitives.
   /// These need to be homeomorphisms locally on the respective primitives, but not globally. This
   /// can cause problems when we try to map a point back from the physical to the computational domain,
   /// but do not know to which primitive the point will belong on the computational domain.
   ///
   /// This method is supposed to test, whether mapping a point from the computational domain to the
   /// physical domain gives the desired coordinates. The base class implements this by mapping the
   /// point from the computation to the physical domain and checking the distance of the two points
   /// via the given threshold. The method is virtual so that children can override it, if need be.
   ///
   /// For details and motivation see issue 184.
#ifdef WALBERLA_DOUBLE_ACCURACY
   virtual bool verifyPointPairing( const Point3D& computationalCoordinates,
                                    const Point3D& physicalCoordinates,
                                    real_t         threshold = real_c( 1e-12 ) ) const;
#else
   virtual bool verifyPointPairing( const Point3D& computationalCoordinates,
                                    const Point3D& physicalCoordinates,
                                    real_t         threshold = real_c( 1e-6 ) ) const;
#endif

   /// Returns whether this Geometry map implements the identity.
   /// Intended to enable optimizations.
   virtual bool isIdentity() const { return false; }

   /// Returns whether this Geometry map implements an affine mapping.
   /// Intended to enable optimizations.
   virtual bool isAffine() const { return false; }

   /// Serialization of a GeometryMap object
   static void serialize( const std::shared_ptr< GeometryMap >& map, walberla::mpi::SendBuffer& sendBuffer );

   /// Deserialization of a GeometryMap object
   static std::shared_ptr< GeometryMap > deserialize( walberla::mpi::RecvBuffer& recvBuffer );

 protected:
   virtual void serializeSubClass( walberla::mpi::SendBuffer& sendBuffer ) const = 0;
};

} // namespace hyteg
