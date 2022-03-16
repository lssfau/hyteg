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
#include <vector>

#include "hyteg/composites/P1BlendingStokesOperator.hpp"
#include "hyteg/composites/P1P1StokesOperator.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p1functionspace/P1Petsc.hpp"
#include "hyteg/p1functionspace/P1ProjectNormalOperator.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"

namespace hyteg {

/// Concatenates a HyTeG operators with the strong free-slip boundary operator.
/// This allows us to enforce a strong free-slip boundary condition, without modifying the operators directly.
///
/// Usage Example:
/// \code
///     using StokesOperator = hyteg::StrongFreeSlipWrapper< hyteg::P1P1StokesOperator, hyteg::P1ProjectNormalOperator, true >;
///     auto stokes = std::make_shared< hyteg::P1P1StokesOperator > ( storage, minLevel, maxLevel );
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
   typedef typename OpType::BlockPreconditioner_T BlockPreconditioner_T;

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

      if ( PreProjection && projFlag_ != None )
      {
         tmp_->assign( {1}, {src}, level, All );
         projOp_->project( *tmp_, level, projFlag_ );
         op_->apply( *tmp_, dst, level, flag );
      }
      else
      {
         op_->apply( src, dst, level, flag );
      }

      if ( projFlag_ != None )
         projOp_->project( dst, level, projFlag_ );
   }

   /// Assemble operator as sparse matrix
   ///
   /// \param mat   a sparse matrix proxy
   /// \param numeratorSrc numerator for determining col indices
   /// \param numeratorDst numerator for determining row indices
   /// \param level level in mesh hierarchy for which local operator is to be assembled
   /// \param flag  determines on which primitives this operator is assembled
   ///
   void toMatrix( const std::shared_ptr< SparseMatrixProxy >&                     mat,
                  const typename OpType::srcType::template FunctionType< idx_t >& numeratorSrc,
                  const typename OpType::dstType::template FunctionType< idx_t >& numeratorDst,
                  uint_t                                                          level,
                  DoFType                                                         flag ) const
   {
      auto matProxyOp = mat->createCopy();
      op_->toMatrix( matProxyOp, numeratorSrc, numeratorDst, level, flag );

      auto matProxyProjectionPost = mat->createCopy();

      projOp_->toMatrix( matProxyProjectionPost, numeratorSrc.uvw(), numeratorDst.uvw(), level, projFlag_ );

      // we need the Id also in the pressure block
      saveIdentityOperator( numeratorDst.p(), matProxyProjectionPost, level, All );

      std::vector< std::shared_ptr< SparseMatrixProxy > > matrices;
      matrices.push_back( matProxyProjectionPost );
      matrices.push_back( matProxyOp );

      if ( PreProjection )
      {
         auto matProxyProjectionPre = mat->createCopy();

         projOp_->toMatrix( matProxyProjectionPre, numeratorSrc.uvw(), numeratorDst.uvw(), level, projFlag_ );

         saveIdentityOperator( numeratorDst.p(), matProxyProjectionPre, level, All );

         matrices.push_back( matProxyProjectionPre );
      }

      mat->createFromMatrixProduct( matrices );
   }

 private:
   std::shared_ptr< OpType >     op_;
   std::shared_ptr< ProjOpType > projOp_;

   DoFType projFlag_;

   std::shared_ptr< typename OpType::dstType > tmp_;
};

} // namespace hyteg
