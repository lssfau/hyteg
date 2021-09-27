/*
 * Copyright (c) 2020 Marcus Mohr.
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

#include "hyteg/operators/Operator.hpp"
#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/forms/P1RowSumForm.hpp"
#include "hyteg/forms/P2RowSumForm.hpp"
#include "hyteg/functions/FunctionIterator.hpp"
#include "hyteg/p1functionspace/P1Petsc.hpp"
#include "hyteg/p1functionspace/VertexDoFFunction.hpp"
#include "hyteg/p2functionspace/P2Petsc.hpp"

namespace hyteg {

using walberla::real_t;

#ifdef HYTEG_BUILD_WITH_PETSC

// As long as we cannot use FunctionIterator< P2Function > we use a specialised external template function
namespace workaround {

template < typename func_T >
void externalDiagonalAssembly( const std::shared_ptr< SparseMatrixProxy >&            mat,
                               const func_T&                                          diagVals,
                               const typename func_T::template FunctionType< idx_t >& numerator,
                               uint_t                                                 level,
                               DoFType                                                flag );

} // namespace workaround

#endif

/// Provides an operator with only "diagonal" values that may change from DoF to DoF
///
/// The DiagonalNonConstantOperator provides an implementation of an operator that is
/// in a matrix sense diagonal. An example would be a lumped mass matrix. The operator
/// wraps itself around an ElementwiseOperator and uses the latter to assemble and obtain the
/// diagonal values.
///
/// The DiagonalNonConstantOperator is intended to work together with RowSumForms, but
/// will work for any other form, too, by just disregarding the non-diagonal matrix entries.
///
/// While currently combined with an underlying ElementwiseOperator the
/// DiagonalNonConstantOperator could also be combined with any other operator,
/// as long as that provides access to the diagonal matrix entries in the form of
/// a HyTeG function.
template < template < class > class opType, class formType, bool InvertDiagonal = false >
class DiagonalNonConstantOperator : public Operator< typename opType< formType >::srcType, typename opType< formType >::dstType >
{
   typedef typename opType< formType >::srcType funcType;

 public:
   DiagonalNonConstantOperator() = delete;

   DiagonalNonConstantOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                size_t                                     minLevel,
                                size_t                                     maxLevel,
                                const std::shared_ptr< formType >&         form )
   : Operator< funcType, funcType >( storage, minLevel, maxLevel )
   , form_( form )
   {
      // generate the operator around which we wrap ourselves
      oper_ = std::make_unique< opType< formType > >( storage, minLevel, maxLevel, *form, false );

      // compute diagonal values and invert if required
      reassemble();
   }

   ~DiagonalNonConstantOperator() override = default;

   /// Reassemble diagonal operator entries
   ///
   /// As the name indicates calling this method triggers
   /// recompuation of the diagonal entries of the matrix
   /// corresponsing to the operator
   void reassemble()
   {
      if ( !InvertDiagonal )
      {
         oper_->computeDiagonalOperatorValues();
      }
      else
      {
         oper_->computeInverseDiagonalOperatorValues();
      }
   }

   void apply( const funcType& src,
               const funcType& dst,
               size_t          level,
               DoFType         flag,
               UpdateType      updateType = Replace ) const override final
   {
      std::shared_ptr< funcType > opVals = InvertDiagonal ? oper_->getInverseDiagonalValues() : oper_->getDiagonalValues();
      if ( updateType == Replace )
      {
         dst.multElementwise( {*opVals, src}, level, flag );
      }
      else
      {
         WALBERLA_ABORT( "DiagonalNonConstantOperator::apply does not support updateType = Add, yet" );
      }
   }

#ifdef HYTEG_BUILD_WITH_PETSC
   void assembleLocalMatrix( const std::shared_ptr< SparseMatrixProxy >&                                 mat,
                             const typename opType< formType >::srcType::template FunctionType< idx_t >& numerator,
                             uint_t                                                                      level,
                             DoFType                                                                     flag ) const
   {
      std::shared_ptr< funcType > opVals = InvertDiagonal ? oper_->getInverseDiagonalValues() : oper_->getDiagonalValues();
      workaround::externalDiagonalAssembly< funcType >( mat, *opVals, numerator, level, flag );
   }
#endif

 private:
   std::shared_ptr< formType >           form_;
   std::unique_ptr< opType< formType > > oper_;
};

/// Diagonal operator for P1 HyTeG forms potentially including blending and/or variable coefficients
typedef DiagonalNonConstantOperator< P1ElementwiseOperator, P1RowSumForm, false > P1BlendingLumpedDiagonalOperator;

/// Diagonal inverse operator for P1 HyTeG forms potentially including blending and/or variable coefficients
typedef DiagonalNonConstantOperator< P1ElementwiseOperator, P1RowSumForm, true > P1BlendingLumpedInverseDiagonalOperator;

/// Diagonal operator for P2 HyTeG forms potentially including blending and/or variable coefficients
typedef DiagonalNonConstantOperator< P2ElementwiseOperator, P2RowSumForm, false > P2BlendingLumpedDiagonalOperator;

/// Diagonal inverse operator for P2 HyTeG forms potentially including blending and/or variable coefficients
typedef DiagonalNonConstantOperator< P2ElementwiseOperator, P2RowSumForm, true > P2BlendingLumpedInverseDiagonalOperator;

// ========================
//  Sparse Matrix Assembly
// ========================
#ifdef HYTEG_BUILD_WITH_PETSC

namespace petsc {

/// Version of createMatrix function for DiagonalNonConstantOperator
template < template < class > class opType, class formType, bool InvertDiagonal = false >
inline void createMatrix(
    const DiagonalNonConstantOperator< opType, formType, InvertDiagonal >&                                          opr,
    const typename DiagonalNonConstantOperator< opType, formType, InvertDiagonal >::template FunctionType< idx_t >& src,
    const typename DiagonalNonConstantOperator< opType, formType, InvertDiagonal >::template FunctionType< idx_t >& dst,
    const std::shared_ptr< SparseMatrixProxy >&                                                                     mat,
    uint_t                                                                                                          level,
    DoFType                                                                                                         flag )
{
   // we don't need to pass same info twice
   WALBERLA_UNUSED( dst );

   // delegate work to operator method
   opr.assembleLocalMatrix( mat, src, level, flag );
}

// We use this specialisations, as the definition above is not sufficient for the
// compiler to recognise that it should call the above and not the templated createMatrix
// function defined in P1Petsc.hpp. Why it works for the createMatrix template functions
// for the P[12]ElementwiseOperator in ElementwiseOperatorPetsc.hpp I do not understand.
// Maybe a missing include somewhere?
template <>
inline void createMatrix( const P1BlendingLumpedDiagonalOperator&     opr,
                          const P1Function< idx_t >&                  src,
                          const P1Function< idx_t >&                  dst,
                          const std::shared_ptr< SparseMatrixProxy >& mat,
                          uint_t                                      level,
                          DoFType                                     flag )
{
   WALBERLA_UNUSED( dst );
   opr.assembleLocalMatrix( mat, src, level, flag );
}

template <>
inline void createMatrix( const P1BlendingLumpedInverseDiagonalOperator& opr,
                          const P1Function< idx_t >&                     src,
                          const P1Function< idx_t >&                     dst,
                          const std::shared_ptr< SparseMatrixProxy >&    mat,
                          uint_t                                         level,
                          DoFType                                        flag )
{
   WALBERLA_UNUSED( dst );
   opr.assembleLocalMatrix( mat, src, level, flag );
}

template <>
inline void createMatrix( const P2BlendingLumpedDiagonalOperator&     opr,
                          const P2Function< idx_t >&                  src,
                          const P2Function< idx_t >&                  dst,
                          const std::shared_ptr< SparseMatrixProxy >& mat,
                          uint_t                                      level,
                          DoFType                                     flag )
{
   WALBERLA_UNUSED( dst );
   opr.assembleLocalMatrix( mat, src, level, flag );
}

template <>
inline void createMatrix( const P2BlendingLumpedInverseDiagonalOperator& opr,
                          const P2Function< idx_t >&                     src,
                          const P2Function< idx_t >&                     dst,
                          const std::shared_ptr< SparseMatrixProxy >&    mat,
                          uint_t                                         level,
                          DoFType                                        flag )
{
   WALBERLA_UNUSED( dst );
   opr.assembleLocalMatrix( mat, src, level, flag );
}

} // namespace petsc
#endif

} // namespace hyteg
