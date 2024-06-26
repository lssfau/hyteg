/*
* Copyright (c) 2017-2022 Nils Kohl.
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

#include "core/DataTypes.h"

#include "hyteg/types/PointND.hpp"

#include "Eigen/Eigen"

namespace hyteg {
namespace dg {

/// \brief Interface for the implementation of the basis functions for DG discretizations.
///
/// The interface provides methods to gather information about the type of basis functions and available polynomial degree.
/// Most crucial is the evaluate() method that evaluates the local polynomial on an element given the corresponding DoFs.
class DGBasisInfo
{
 public:
   /// \brief Returns the minimum polynomial degree for which basis functions are implemented.
   virtual uint_t minPolynomialDegree() const = 0;

   /// \brief Returns the maximum polynomial degree for which basis functions are implemented.
   virtual uint_t maxPolynomialDegree() const = 0;

   /// \brief Returns the number of DoFs per element for the passed polynomial degree and type of the element.
   virtual uint_t numDoFsPerElement( uint_t dim, uint_t degree ) const = 0;

   /// \brief Returns the degree of the quadrature rule that is used for the evaluation of the linear functional \f$\int_T f \phi_i\f$.
   virtual uint_t quadratureDegreeForLinearFunctional() const = 0;

   /// \brief Evaluates the polynomial on the reference triangle.
   ///
   /// \param degree degree of the piecewise polynomials
   /// \param pos    where to evaluate on the micro-element (in reference space)
   /// \param dofs   DoFs that correspond to the basis functions
   /// \param value  value of the polynomial
   virtual void evaluate( uint_t degree, const Point2D& pos, const std::vector< real_t >& dofs, real_t& value ) const = 0;

   /// \brief Evaluates the polynomial on the reference tetrahedron.
   ///
   /// \param degree degree of the piecewise polynomials
   /// \param pos    where to evaluate on the micro-element (in reference space)
   /// \param dofs   DoFs that correspond to the basis functions
   /// \param value  value of the polynomial
   virtual void evaluate( uint_t degree, const Point3D& pos, const std::vector< real_t >& dofs, real_t& value ) const = 0;

   /// \brief Evaluates the linear functional
   ///
   /// \f{align*}{
   ///   l( v ) = \int_T f v
   /// \f}
   ///
   /// by quadrature over the passed face element \f$T\f$ for all basis functions, i.e. it returns an approximation to
   ///
   /// \f{align*}{
   ///   \int_T f \phi_i.
   /// \f}
   ///
   /// \param degree           degree of the piecewise polynomials
   /// \param coords           coordinates of the affine element (computational space)
   /// \param f                function to multiply with the basis functions
   /// \param value            result of the integration (all DoFs)
   virtual void integrateBasisFunction( uint_t                                           degree,
                                        const std::array< Point2D, 3 >&                  coords,
                                        const std::function< real_t( const Point3D& ) >& f,
                                        std::vector< real_t >&                           values ) = 0;

   /// \brief Evaluates the linear functional
   ///
   /// \f{align*}{
   ///   l( v ) = \int_T f v
   /// \f}
   ///
   /// by quadrature over the passed cell element \f$T\f$ for all basis functions, i.e. it returns an approximation to
   ///
   /// \f{align*}{
   ///   \int_T f \phi_i.
   /// \f}
   ///
   /// \param degree           degree of the piecewise polynomials
   /// \param coords           coordinates of the affine element (computational space)
   /// \param f                function to multiply with the basis functions
   /// \param value            result of the integration (all DoFs)
   virtual void integrateBasisFunction( uint_t                                           degree,
                                        const std::array< Point3D, 4 >&                  coords,
                                        const std::function< real_t( const Point3D& ) >& f,
                                        std::vector< real_t >&                           values ) = 0;
};

} // namespace dg
} // namespace hyteg