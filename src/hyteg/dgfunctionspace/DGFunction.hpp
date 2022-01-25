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

#include "hyteg/dgfunctionspace/DGForm.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"
#include "hyteg/volumedofspace/VolumeDoFFunction.hpp"

namespace hyteg {
namespace dg {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

template < typename ValueType >
class DGFunction final : public Function< DGFunction< ValueType > >
{
 public:
   typedef ValueType valueType;

   template < typename VType >
   using FunctionType = DGFunction< VType >;

   /// \brief Instantiates a DGFunction.
   ///
   /// \param name              name of the function for output / debugging purposes
   /// \param storage           PrimitiveStorage that this lives on
   /// \param minLevel          min refinement level to allocate
   /// \param maxLevel          max refinement level to allocate
   /// \param basis             polynomial basis information through corresponding DGBasisInfo object
   /// \param initialPolyDegree polynomial degree to initialize with globally
   DGFunction( const std::string&                         name,
               const std::shared_ptr< PrimitiveStorage >& storage,
               uint_t                                     minLevel,
               uint_t                                     maxLevel,
               const std::shared_ptr< DGBasisInfo >&      basis,
               int                                        initialPolyDegree );

   /// \brief Assigns a linear combination of multiple VolumeDoFFunctions to this.
   void assign( const std::vector< ValueType >&                                               scalars,
                const std::vector< std::reference_wrapper< const DGFunction< ValueType > > >& functions,
                uint_t                                                                        level );

   /// \brief Evaluate finite element function at a specific coordinates.
   ///
   /// In a parallel setting, the specified coordinate might not lie in the local subdomain.
   ///
   /// Evaluation is performed in two steps:
   ///
   ///   1. For all volume primitives of the local subdomain:
   ///      If a point-tet (point-triangle in 2D) inclusion test succeeds, the function returns true
   ///      and the finite-element function is evaluated.
   ///
   ///   2. Skipped, if radius is negative.
   ///      For all volume primitives of the local subdomain:
   ///      A sphere-tet (circle-triangle) intersection
   ///      test is performed, if successful returns true, and the finite-element function is extrapolated
   ///      to the specified coordinate and evaluated, depending on the radius, this might introduce (large) errors.
   ///
   /// If both tests fail, this function returns false, and no evaluation is performed (i.e. the returned, evaluated
   /// value is not set to anything meaningful).
   ///
   /// Note that two parallel processes that return true, may return _different_ values.
   ///
   /// No communication is performed in this function.
   /// -> Does not need to be called collectively.
   /// -> Different values are returned on each process.
   ///
   /// \param coordinates where the function shall be evaluated
   /// \param level refinement level
   /// \param value function value at the coordinate if search was successful
   /// \param searchToleranceRadius radius of the sphere (circle) for the second search phase, skipped if negative
   /// \return true if the function was evaluated successfully, false otherwise
   ///
   bool evaluate( const Point3D& coordinates, uint_t level, ValueType& value, real_t searchToleranceRadius = 1e-05 ) const;

   /// Evaluates the linear functional
   ///
   ///   l( v ) = \int_\Omega f * v
   ///
   /// by integration over the local basis functions and writes the result into the vector, i.e.
   ///
   ///   u_i <- \int_T f * \phi_i
   ///
   /// where \phi_i is the basis function associated with the DoF u_i and f a given analytical function.
   void evaluateLinearFunctional( const std::function< real_t( const Point3D& ) >& f, uint_t level );

   /// \brief Returns the internally stored VolumeDoFFunction.
   [[nodiscard]] std::shared_ptr< volumedofspace::VolumeDoFFunction< ValueType > > volumeDoFFunction() const
   {
      return volumeDoFFunction_;
   }

   /// \brief Returns the associated DGBasisInfo object.
   [[nodiscard]] std::shared_ptr< DGBasisInfo > basis() const { return basis_; }

   /// \brief Returns the polynomial degree of the function on the passed primitive.
   [[nodiscard]] int polynomialDegree( const PrimitiveID& primitiveID ) const
   {
      WALBERLA_ASSERT( storage_->primitiveExistsLocally( primitiveID ),
                       "No information on the polynomial degree on this primitive (primitive does not exist locally)." );
      WALBERLA_ASSERT( polyDegreesPerPrimitive_.count( primitiveID ) > 0,
                       "No information on the polynomial degree on this primitive (no information stored in the map)." );
      return polyDegreesPerPrimitive_.at( primitiveID );
   }

   /// \brief Assigns unique values to all data points. Increments by 1 per DoF.
   ///
   /// Mainly used for sparse matrix assembly to get global integer identifiers for all DoFs.
   ///
   /// \param level refinement level to enumerate
   void enumerate( uint_t level ) const;

   /// \brief Like enumerate() but starting with a value given by the offset parameter.
   ///
   /// \param level  refinement level to enumerate
   /// \param offset value to start from, the next "free" value is returned
   void enumerate( uint_t level, ValueType& offset ) const;

   template < typename OtherValueType >
   void copyBoundaryConditionFromFunction( const DGFunction< OtherValueType >& other )
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "DGFunction: BCs are not copied!" );
   }

   /// \brief Returns the number of DoFs that are allocated on this process.
   uint_t getNumberOfLocalDoFs( uint_t level ) const;

   /// \brief Returns the number of DoFs. Performs global reduction, must be called collectively.
   uint_t getNumberOfGlobalDoFs( uint_t level ) const;

 private:
   using Function< DGFunction< ValueType > >::communicators_;
   using Function< DGFunction< ValueType > >::additiveCommunicators_;

   std::string                         name_;
   std::shared_ptr< PrimitiveStorage > storage_;
   uint_t                              minLevel_;
   uint_t                              maxLevel_;
   std::shared_ptr< DGBasisInfo >      basis_;

   std::shared_ptr< volumedofspace::VolumeDoFFunction< ValueType > > volumeDoFFunction_;

   std::map< PrimitiveID, uint_t > polyDegreesPerPrimitive_;
};

} // namespace dg

void createVectorFromFunction( const dg::DGFunction< real_t >&       function,
                               const dg::DGFunction< idx_t >&        numerator,
                               const std::shared_ptr< VectorProxy >& vec,
                               uint_t                                level,
                               DoFType                               flag );

void createFunctionFromVector( const dg::DGFunction< real_t >&       function,
                               const dg::DGFunction< idx_t >&        numerator,
                               const std::shared_ptr< VectorProxy >& vec,
                               uint_t                                level,
                               DoFType                               flag );

void applyDirichletBC( const dg::DGFunction< idx_t >& numerator, std::vector< idx_t >& mat, uint_t level );

} // namespace hyteg