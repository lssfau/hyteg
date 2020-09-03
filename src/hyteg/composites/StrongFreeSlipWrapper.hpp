/*
 * Copyright (c) 2020 Andreas Wagner
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

#include <type_traits>

#include "hyteg/Operator.hpp"

namespace hyteg {

/// Concatenates a HyTeG operators with the strong free-slip boundary operator.
/// This allows us to enforce a strong free-slip boundary condition, without modifying the operators directly.
///
/// Usage Example:
/// \code
///     using StokesOperator = hyteg::StrongFreeSlipWrapper< hyteg::P1StokesOperator, hyteg::P1ProjectNormalOperator, true >;
///     auto stokes = std::make_shared< hyteg::P1StokesOperator > ( storage, minLevel, maxLevel );
///     auto normals = [](auto, Point3D & n) { n = Point3D({0, -1}); };
///     auto projection = std::make_shared< hyteg::P1ProjectNormalOperator > ( storage, minLevel, maxLevel, normals );
///     StokesOperator L( stokes, projection, FreeslipBoundary );
///
///     auto solver = hyteg::solvertemplates::stokesMinResSolver< StokesOperator >( storage, maxLevel, 1e-6, 1000 );
/// \endcode
///
/// \tparam OpType          Type of operator on which we want to impose free-slip.
/// \tparam ProjOpType      Type of the normal projection operator, which enforces free slip
/// \tparam PreProjection   Should the projection operator also already be applied to the src vector? This is not necessary if it is already consistent!
template < typename OpType, typename ProjOpType, bool PreProjection = false >
class StrongFreeSlipWrapper : public Operator< typename OpType::srcType, typename OpType::dstType >
{
 public:
   StrongFreeSlipWrapper( std::shared_ptr< OpType > op, std::shared_ptr< ProjOpType > projOp, DoFType projFlag )
   : Operator< typename OpType::srcType, typename OpType::dstType >( op->getStorage(), op->getMinLevel(), op->getMaxLevel() )
   , op_( op )
   , projOp_( projOp )
   , projFlag_( projFlag )
   , tmp_( PreProjection ?
               std::make_shared< typename OpType::srcType >( "tmp", op->getStorage(), op->getMinLevel(), op->getMaxLevel() ) :
               nullptr )
   {}

   void apply( const typename OpType::srcType& src,
               const typename OpType::dstType& dst,
               size_t                          level,
               DoFType                         flag,
               UpdateType                      updateType = Replace ) const
   {
      WALBERLA_CHECK( updateType == Replace, "Operator concatenation only supported for updateType Replace" );

      if ( PreProjection )
      {
         tmp_->assign( { 1 }, { src }, level, All );
         projOp_->apply( *tmp_, level, projFlag_ );
         op_->apply( *tmp_, dst, level, flag );
      }
      else
      {
         op_->apply( src, dst, level, flag );
      }

      projOp_->apply( dst, level, projFlag_ );
   }

 private:
   std::shared_ptr< OpType >     op_;
   std::shared_ptr< ProjOpType > projOp_;

   DoFType projFlag_;

   std::shared_ptr< typename OpType::dstType > tmp_;
};

} // namespace hyteg
