/*
 * Copyright (c) 2024-2025 Andreas Burkhart.
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

#include "hyteg/functions/FunctionWrapper.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/types/types.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

// These macros allow you to create constexpr functions that check if a given operator has a method via an SFINAE approach (see examples below).
// This should be easily adaptable to other contexts

#define UNPACK_MACRO_ARGUMENTS( ... ) __VA_ARGS__
#define DEFINE_CLASS_METHOD_CHECK( functionName, returnType, functionSignature, ... )                                    \
   template < typename OperatorType, typename = void >                                                                   \
   struct check_##functionName##_Struct : std::false_type                                                                \
   {};                                                                                                                   \
                                                                                                                         \
   template < typename OperatorType >                                                                                    \
   struct check_##functionName##_Struct<                                                                                 \
       OperatorType,                                                                                                     \
       typename std::enable_if< std::is_same< decltype( &OperatorType::functionName ),                                   \
                                              returnType ( OperatorType::* )( UNPACK_MACRO_ARGUMENTS functionSignature ) \
                                                  __VA_ARGS__ >::value >::type > : std::true_type                        \
   {};                                                                                                                   \
                                                                                                                         \
   template < class OperatorType >                                                                                       \
   constexpr bool has_##functionName()                                                                                   \
   {                                                                                                                     \
      return check_##functionName##_Struct< OperatorType >::value;                                                       \
   }

namespace hyteg {

// Some SFINAE tricks that can be used to determine if an operator (or class in general) provides a given method
namespace SFINAE {

DEFINE_CLASS_METHOD_CHECK( applyScaled,
                           void,
                           ( const typename OperatorType::srcType::valueType&,
                             const typename OperatorType::srcType&,
                             const typename OperatorType::dstType&,
                             uint_t,
                             hyteg::DoFType,
                             hyteg::UpdateType ),
                           const );

DEFINE_CLASS_METHOD_CHECK( gemv,
                           void,
                           ( const typename OperatorType::srcType::valueType&,
                             const typename OperatorType::srcType&,
                             const typename OperatorType::dstType::valueType&,
                             const typename OperatorType::dstType&,
                             uint_t,
                             hyteg::DoFType ),
                           const );

DEFINE_CLASS_METHOD_CHECK( toMatrixScaled,
                           void,
                           ( const typename OperatorType::srcType::valueType&,
                             const std::shared_ptr< hyteg::SparseMatrixProxy >&,
                             const typename OperatorType::srcType::template FunctionType< hyteg::idx_t >&,
                             const typename OperatorType::dstType::template FunctionType< hyteg::idx_t >&,
                             uint_t,
                             hyteg::DoFType ),
                           const );

DEFINE_CLASS_METHOD_CHECK( computeInverseDiagonalOperatorValuesScaled,
                           void,
                           (const typename OperatorType::srcType::valueType&), );

DEFINE_CLASS_METHOD_CHECK( computeInverseDiagonalOperatorValues, void, (), );

DEFINE_CLASS_METHOD_CHECK( getInverseDiagonalValues, std::shared_ptr< typename OperatorType::srcType >, (), const );

DEFINE_CLASS_METHOD_CHECK( computeLumpedInverseDiagonalOperatorValues, void, (), );

DEFINE_CLASS_METHOD_CHECK( getLumpedInverseDiagonalValues, std::shared_ptr< typename OperatorType::srcType >, (), const );

DEFINE_CLASS_METHOD_CHECK( smooth_jac,
                           void,
                           ( const typename OperatorType::srcType&,
                             const typename OperatorType::srcType&,
                             const typename OperatorType::srcType&,
                             typename OperatorType::srcType::valueType,
                             uint_t,
                             hyteg::DoFType ),
                           const );

DEFINE_CLASS_METHOD_CHECK( smooth_jac_scaled,
                           void,
                           ( const typename OperatorType::srcType::valueType&,
                             const typename OperatorType::srcType&,
                             const typename OperatorType::srcType&,
                             const typename OperatorType::srcType&,
                             typename OperatorType::srcType::valueType,
                             uint_t,
                             hyteg::DoFType ),
                           const );

} // namespace SFINAE

} // namespace hyteg