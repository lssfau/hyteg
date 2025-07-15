/*
* Copyright (c) 2023-2025 Andreas Burkhart
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

// This file has been generated with the AHFC. If buggy try fixing the generator itself.

#pragma once
#include "elementwise_dof_value_operator/P2ElementwiseDoFValueOperator.hpp"

namespace hyteg {

using walberla::real_t;

/// Adiabatic heating operator for the TALA
///
/// Intended for LHS use.
///
/// The gravity vector is defined via std::functions. For the annulus / spherical shell:
/// g := -x / norm(x).
///
/// Weak formulation:
///
///    T: trial function (scalar space: P2)
///    w: test function (scalar space: P2)
///    u: coefficient (vector space: P2)
///    g: arbitrary std::function vector coefficient
///    k: arbitrary std::function coefficient
///
///    - ∫ k T(u · g) w
///
/// Blending: True
/// Quadrature degree: 6
class p2_adiabatic_heat_u_p2_blending_q6_ElementwiseOperator
: public P2ElementwiseDoFValueOperator< P2Function< real_t >, P2Function< real_t >, P2Function< real_t > >
{
   typedef P2ElementwiseDoFValueOperator< P2Function< real_t >, P2Function< real_t >, P2Function< real_t > > baseClass;
   typedef P2Function< real_t >                                                                              srcFunctionType;
   typedef P2Function< real_t >                                                                              dstFunctionType;

 public:
   p2_adiabatic_heat_u_p2_blending_q6_ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                           size_t                                     minLevel,
                                                           size_t                                     maxLevel,
                                                           P2Function< real_t >&                      uX,
                                                           P2Function< real_t >&                      uY,
                                                           P2Function< real_t >&                      uZ,
                                                           std::function< real_t( const Point3D& ) >  gX,
                                                           std::function< real_t( const Point3D& ) >  gY,
                                                           std::function< real_t( const Point3D& ) >  gZ,
                                                           std::function< real_t( const Point3D& ) >  adiabatic_heat_scaling )
   : p2_adiabatic_heat_u_p2_blending_q6_ElementwiseOperator( storage,
                                                             minLevel,
                                                             maxLevel,
                                                             false,
                                                             uX,
                                                             uY,
                                                             uZ,
                                                             gX,
                                                             gY,
                                                             gZ,
                                                             adiabatic_heat_scaling )
   {}
   p2_adiabatic_heat_u_p2_blending_q6_ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                           size_t                                     minLevel,
                                                           size_t                                     maxLevel,
                                                           bool                                       needsInverseDiagEntries,
                                                           P2Function< real_t >&                      uX,
                                                           P2Function< real_t >&                      uY,
                                                           P2Function< real_t >&                      uZ,
                                                           std::function< real_t( const Point3D& ) >  gX,
                                                           std::function< real_t( const Point3D& ) >  gY,
                                                           std::function< real_t( const Point3D& ) >  gZ,
                                                           std::function< real_t( const Point3D& ) >  adiabatic_heat_scaling )
   : P2ElementwiseDoFValueOperator< P2Function< real_t >, P2Function< real_t >, P2Function< real_t > >( storage,
                                                                                                        minLevel,
                                                                                                        maxLevel,
                                                                                                        false,
                                                                                                        uX,
                                                                                                        uY,
                                                                                                        uZ )
   , gX_( gX )
   , gY_( gY )
   , gZ_( gZ )
   , adiabatic_heat_scaling_( adiabatic_heat_scaling )
   {
      if ( needsInverseDiagEntries )
      {
         computeInverseDiagonalOperatorValues();
      }
   }

   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 6, 6 >& elMat ) const override;

   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 6 >& elMat ) const override;

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 10, 10 >& elMat ) const override;

   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 10 >& elMat ) const override;

   // Expose functions to SFINAE.hpp
   // PLEASE DO NOT REMOVE!

   void computeDiagonalOperatorValues() override { baseClass::computeDiagonalOperatorValues(); }

   void computeInverseDiagonalOperatorValues() override { baseClass::computeInverseDiagonalOperatorValues(); }

   void computeLumpedDiagonalOperatorValues() override { baseClass::computeLumpedDiagonalOperatorValues(); }

   void computeLumpedInverseDiagonalOperatorValues() override { baseClass::computeLumpedInverseDiagonalOperatorValues(); }

   void computeInverseDiagonalOperatorValuesScaled( const real_t& alpha ) override
   {
      baseClass::computeInverseDiagonalOperatorValuesScaled( alpha );
   }

   void computeLumpedInverseDiagonalOperatorValuesScaled( const real_t& alpha ) override
   {
      baseClass::computeLumpedInverseDiagonalOperatorValuesScaled( alpha );
   }

   void computeDiagonalOperatorValuesScaled( const real_t& alpha ) override
   {
      baseClass::computeDiagonalOperatorValuesScaled( alpha );
   }

   void computeLumpedDiagonalOperatorValuesScaled( const real_t& alpha ) override
   {
      baseClass::computeLumpedDiagonalOperatorValuesScaled( alpha );
   }

   std::shared_ptr< srcFunctionType > getDiagonalValues() const override { return baseClass::getDiagonalValues(); }

   std::shared_ptr< srcFunctionType > getInverseDiagonalValues() const override { return baseClass::getInverseDiagonalValues(); }

   std::shared_ptr< srcFunctionType > getLumpedDiagonalValues() const override { return baseClass::getLumpedDiagonalValues(); }

   std::shared_ptr< srcFunctionType > getLumpedInverseDiagonalValues() const override
   {
      return baseClass::getLumpedInverseDiagonalValues();
   }

   void smooth_jac_scaled( const real_t&          alpha,
                           const srcFunctionType& dst,
                           const srcFunctionType& rhs,
                           const srcFunctionType& src,
                           real_t                 omega,
                           uint_t                 level,
                           DoFType                flag ) const override
   {
      baseClass::smooth_jac_scaled( alpha, dst, rhs, src, omega, level, flag );
   }

   void smooth_jac( const srcFunctionType& dst,
                    const srcFunctionType& rhs,
                    const srcFunctionType& src,
                    real_t                 omega,
                    uint_t                 level,
                    DoFType                flag ) const override
   {
      baseClass::smooth_jac( dst, rhs, src, omega, level, flag );
   }

   void gemv( const real_t&          alpha,
              const srcFunctionType& src,
              const real_t&          beta,
              const dstFunctionType& dst,
              uint_t                 level,
              DoFType                flag ) const override
   {
      baseClass::gemv( alpha, src, beta, dst, level, flag );
   }

   void applyScaled( const real_t&          alpha,
                     const srcFunctionType& src,
                     const dstFunctionType& dst,
                     uint_t                 level,
                     DoFType                flag,
                     UpdateType             updateType = Replace ) const override
   {
      baseClass::applyScaled( alpha, src, dst, level, flag, updateType );
   }

   void apply( const srcFunctionType& src,
               const dstFunctionType& dst,
               uint_t                 level,
               DoFType                flag,
               UpdateType             updateType = Replace ) const override
   {
      baseClass::apply( src, dst, level, flag, updateType );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >&                            mat,
                  const typename srcFunctionType::template FunctionType< hyteg::idx_t >& src,
                  const typename dstFunctionType::template FunctionType< hyteg::idx_t >& dst,
                  uint_t                                                                 level,
                  DoFType                                                                flag ) const override
   {
      baseClass::toMatrix( mat, src, dst, level, flag );
   }

   void toMatrixScaled( const real_t&                                                          alpha,
                        const std::shared_ptr< SparseMatrixProxy >&                            mat,
                        const typename srcFunctionType::template FunctionType< hyteg::idx_t >& src,
                        const typename dstFunctionType::template FunctionType< hyteg::idx_t >& dst,
                        uint_t                                                                 level,
                        DoFType                                                                flag ) const override
   {
      baseClass::toMatrixScaled( alpha, mat, src, dst, level, flag );
   }

   void computeAndStoreLocalElementMatrices() override { baseClass::computeAndStoreLocalElementMatrices(); }

 private:
   std::function< real_t( const Point3D& ) > gX_;
   std::function< real_t( const Point3D& ) > gY_;
   std::function< real_t( const Point3D& ) > gZ_;
   std::function< real_t( const Point3D& ) > adiabatic_heat_scaling_;
};

} // namespace hyteg