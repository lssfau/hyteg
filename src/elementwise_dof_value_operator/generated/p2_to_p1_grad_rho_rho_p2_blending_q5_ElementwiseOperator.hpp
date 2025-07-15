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
#include "elementwise_dof_value_operator/P2ToP1ElementwiseDoFValueOperator.hpp"

namespace hyteg {

using walberla::real_t;

/// Rho stokes operator for the compressible case.
///
/// Intended for RHS use.
///
/// Weak formulation:
///
///     u: trial function (vectorial space: P2)
///     v: test function (scalar space: P1)
///     rho: coefficient (scalar space: P2)
///
///     ∫ ( ∇rho/rho · u ) v
///
/// Blending: True
/// Quadrature degree: 5
class p2_to_p1_grad_rho_rho_p2_0_blending_q5_ElementwiseOperator
: public P2ToP1ElementwiseDoFValueOperator< P2Function< real_t > >
{
   typedef P2ToP1ElementwiseDoFValueOperator< P2Function< real_t > > baseClass;
   typedef P2Function< real_t >                                      srcFunctionType;
   typedef P1Function< real_t >                                      dstFunctionType;

 public:
   p2_to_p1_grad_rho_rho_p2_0_blending_q5_ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                               size_t                                     minLevel,
                                                               size_t                                     maxLevel,
                                                               P2Function< real_t >&                      rho )
   : P2ToP1ElementwiseDoFValueOperator< P2Function< real_t > >( storage, minLevel, maxLevel, rho )
   {}

   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 6 >& elMat ) const override;

   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 6 >& elMat ) const override;

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 10 >& elMat ) const override;

   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 10 >& elMat ) const override;

   // Expose functions to SFINAE.hpp
   // PLEASE DO NOT REMOVE!

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
};

/// Rho stokes operator for the compressible case.
///
/// Intended for RHS use.
///
/// Weak formulation:
///
///     u: trial function (vectorial space: P2)
///     v: test function (scalar space: P1)
///     rho: coefficient (scalar space: P2)
///
///     ∫ ( ∇rho/rho · u ) v
///
/// Blending: True
/// Quadrature degree: 5
class p2_to_p1_grad_rho_rho_p2_1_blending_q5_ElementwiseOperator
: public P2ToP1ElementwiseDoFValueOperator< P2Function< real_t > >
{
   typedef P2ToP1ElementwiseDoFValueOperator< P2Function< real_t > > baseClass;
   typedef P2Function< real_t >                                      srcFunctionType;
   typedef P1Function< real_t >                                      dstFunctionType;

 public:
   p2_to_p1_grad_rho_rho_p2_1_blending_q5_ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                               size_t                                     minLevel,
                                                               size_t                                     maxLevel,
                                                               P2Function< real_t >&                      rho )
   : P2ToP1ElementwiseDoFValueOperator< P2Function< real_t > >( storage, minLevel, maxLevel, rho )
   {}

   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 6 >& elMat ) const override;

   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 6 >& elMat ) const override;

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 10 >& elMat ) const override;

   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 10 >& elMat ) const override;

   // Expose functions to SFINAE.hpp
   // PLEASE DO NOT REMOVE!

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
};

/// Rho stokes operator for the compressible case.
///
/// Intended for RHS use.
///
/// Weak formulation:
///
///     u: trial function (vectorial space: P2)
///     v: test function (scalar space: P1)
///     rho: coefficient (scalar space: P2)
///
///     ∫ ( ∇rho/rho · u ) v
///
/// Blending: True
/// Quadrature degree: 5
class p2_to_p1_grad_rho_rho_p2_2_blending_q5_ElementwiseOperator
: public P2ToP1ElementwiseDoFValueOperator< P2Function< real_t > >
{
   typedef P2ToP1ElementwiseDoFValueOperator< P2Function< real_t > > baseClass;
   typedef P2Function< real_t >                                      srcFunctionType;
   typedef P1Function< real_t >                                      dstFunctionType;

 public:
   p2_to_p1_grad_rho_rho_p2_2_blending_q5_ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                               size_t                                     minLevel,
                                                               size_t                                     maxLevel,
                                                               P2Function< real_t >&                      rho )
   : P2ToP1ElementwiseDoFValueOperator< P2Function< real_t > >( storage, minLevel, maxLevel, rho )
   {}

   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 6 >& elMat ) const override;

   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 6 >& elMat ) const override;

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 10 >& elMat ) const override;

   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 10 >& elMat ) const override;

   // Expose functions to SFINAE.hpp
   // PLEASE DO NOT REMOVE!

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
};

} // namespace hyteg