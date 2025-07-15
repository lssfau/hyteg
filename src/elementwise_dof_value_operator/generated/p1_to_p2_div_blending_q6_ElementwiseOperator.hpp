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
#include "elementwise_dof_value_operator/P1ToP2ElementwiseDoFValueOperator.hpp"

// This operator is used to test the functionality of the P1ToP2ElementwiseDoFValueOperators. It should not be used elsewhere due to unnecessary dummy function usage. DO NOT REMOVE.

namespace hyteg {

using walberla::real_t;

/// Divergence.
/// 
/// Weak formulation:
/// 
///     u: trial function (vectorial space: P1)
///     v: test function  (scalar space:    P2)
/// 
///     ∫ - ( ∇ · u ) v
/// 
/// Blending: True
/// Quadrature degree: 6
class p1_to_p2_div_0_blending_q6_ElementwiseOperator : public P1ToP2ElementwiseDoFValueOperator< P2Function< real_t > >
{
   typedef P1ToP2ElementwiseDoFValueOperator< P2Function< real_t > > baseClass;
   typedef P1Function< real_t >                                      srcFunctionType;
   typedef P2Function< real_t >                                      dstFunctionType;

 public:
   p1_to_p2_div_0_blending_q6_ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                   size_t                                     minLevel,
                                                   size_t                                     maxLevel,
                                                   P2Function< real_t >&                      dummy )
   : P1ToP2ElementwiseDoFValueOperator< P2Function< real_t > >( storage, minLevel, maxLevel, dummy )
   {}

   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 6, 3 >& elMat ) const override;

   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const override;

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 10, 4 >& elMat ) const override;

   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const override;

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

/// Divergence.
/// 
/// Weak formulation:
/// 
///     u: trial function (vectorial space: P1)
///     v: test function  (scalar space:    P2)
/// 
///     ∫ - ( ∇ · u ) v
/// 
/// Blending: True
/// Quadrature degree: 6
class p1_to_p2_div_1_blending_q6_ElementwiseOperator : public P1ToP2ElementwiseDoFValueOperator< P2Function< real_t > >
{
   typedef P1ToP2ElementwiseDoFValueOperator< P2Function< real_t > > baseClass;
   typedef P1Function< real_t >                                      srcFunctionType;
   typedef P2Function< real_t >                                      dstFunctionType;

 public:
   p1_to_p2_div_1_blending_q6_ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                   size_t                                     minLevel,
                                                   size_t                                     maxLevel,
                                                   P2Function< real_t >&                      dummy )
   : P1ToP2ElementwiseDoFValueOperator< P2Function< real_t > >( storage, minLevel, maxLevel, dummy )
   {}

   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 6, 3 >& elMat ) const override;

   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const override;

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 10, 4 >& elMat ) const override;

   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const override;

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

/// Divergence.
/// 
/// Weak formulation:
/// 
///     u: trial function (vectorial space: P1)
///     v: test function  (scalar space:    P2)
/// 
///     ∫ - ( ∇ · u ) v
/// 
/// Blending: True
/// Quadrature degree: 6
class p1_to_p2_div_2_blending_q6_ElementwiseOperator : public P1ToP2ElementwiseDoFValueOperator< P2Function< real_t > >
{
   typedef P1ToP2ElementwiseDoFValueOperator< P2Function< real_t > > baseClass;
   typedef P1Function< real_t >                                      srcFunctionType;
   typedef P2Function< real_t >                                      dstFunctionType;

 public:
   p1_to_p2_div_2_blending_q6_ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                   size_t                                     minLevel,
                                                   size_t                                     maxLevel,
                                                   P2Function< real_t >&                      dummy )
   : P1ToP2ElementwiseDoFValueOperator< P2Function< real_t > >( storage, minLevel, maxLevel, dummy )
   {}

   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 6, 3 >& elMat ) const override;

   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const override;

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 10, 4 >& elMat ) const override;

   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const override;

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