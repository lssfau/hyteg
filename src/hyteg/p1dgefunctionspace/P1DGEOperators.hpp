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

#include "hyteg/dgfunctionspace/DGVectorLaplaceForm.hpp"
#include "hyteg/dgfunctionspace/DGVectorMassForm.hpp"
#include "hyteg/mixedoperators/P0ScalarToP1VectorOperator.hpp"
#include "hyteg/mixedoperators/P1VectorToP0ScalarOperator.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/operators/VectorLaplaceOperator.hpp"
#include "hyteg/operators/VectorMassOperator.hpp"
#include "hyteg/p0functionspace/P0Operator.hpp"
#include "hyteg/p1dgefunctionspace/P1DGEFunction.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"

namespace hyteg {

template < typename ValueType >
class P1DGEMassOperator final : public Operator< P1DGEFunction< real_t >, P1DGEFunction< real_t > >
{
 public:
   P1DGEMassOperator( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : Operator< P1DGEFunction< real_t >, P1DGEFunction< real_t > >( storage, minLevel, maxLevel )
   , cg_to_cg_coupling_( storage, minLevel, maxLevel )
   , eg_to_cg_coupling_( storage, minLevel, maxLevel )
   , cg_to_eg_coupling_( storage, minLevel, maxLevel )
   , eg_to_eg_coupling_( storage, minLevel, maxLevel, std::make_shared< dg::DGVectorMassFormEDGEDG >() )
   {}

   void apply( const P1DGEFunction< real_t >& src,
               const P1DGEFunction< real_t >& dst,
               size_t                         level,
               DoFType                        flag,
               UpdateType                     updateType ) const override
   {
      eg_to_cg_coupling_.apply( *src.getDiscontinuousPart(), *dst.getConformingPart(), level, flag, updateType );
      cg_to_cg_coupling_.apply( *src.getConformingPart(), *dst.getConformingPart(), level, flag, Add );

      cg_to_eg_coupling_.apply( *src.getConformingPart(), *dst.getDiscontinuousPart(), level, flag, updateType );
      eg_to_eg_coupling_.apply( *src.getDiscontinuousPart(), *dst.getDiscontinuousPart(), level, flag, Add );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P1DGEFunction< idx_t >&               src,
                  const P1DGEFunction< idx_t >&               dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      communication::syncVectorFunctionBetweenPrimitives( *src.getConformingPart(), level );
      communication::syncVectorFunctionBetweenPrimitives( *dst.getConformingPart(), level );
      src.getDiscontinuousPart()->communicate( level );
      dst.getDiscontinuousPart()->communicate( level );

      cg_to_cg_coupling_.toMatrix( mat, *src.getConformingPart(), *dst.getConformingPart(), level, flag );
      eg_to_cg_coupling_.toMatrix( mat, *src.getDiscontinuousPart(), *dst.getConformingPart(), level, flag );
      cg_to_eg_coupling_.toMatrix( mat, *src.getConformingPart(), *dst.getDiscontinuousPart(), level, flag );
      eg_to_eg_coupling_.toMatrix( mat, *src.getDiscontinuousPart(), *dst.getDiscontinuousPart(), level, flag );
   }

 private:
   P1ConstantVectorMassOperator                  cg_to_cg_coupling_;
   P1ToP0ConstantP1EDGVectorMassCouplingOperator cg_to_eg_coupling_;
   P0ToP1ConstantP1EDGVectorMassCouplingOperator eg_to_cg_coupling_;
   P0Operator< dg::DGVectorMassFormEDGEDG >      eg_to_eg_coupling_;
};

template < typename ValueType >
class P1DGELaplaceOperator final : public Operator< P1DGEFunction< real_t >, P1DGEFunction< real_t > >
{
 public:
   P1DGELaplaceOperator( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : Operator< P1DGEFunction< real_t >, P1DGEFunction< real_t > >( storage, minLevel, maxLevel )
   , cg_to_cg_coupling_( storage, minLevel, maxLevel )
   , eg_to_cg_coupling_( storage, minLevel, maxLevel )
   , cg_to_eg_coupling_( storage, minLevel, maxLevel )
   , eg_to_eg_coupling_( storage, minLevel, maxLevel, std::make_shared< dg::DGVectorLaplaceFormEDGEDG >() )
   {}

   void apply( const P1DGEFunction< real_t >& src,
               const P1DGEFunction< real_t >& dst,
               size_t                         level,
               DoFType                        flag,
               UpdateType                     updateType ) const override
   {
      WALBERLA_ABORT( "Not implemented." );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P1DGEFunction< idx_t >&               src,
                  const P1DGEFunction< idx_t >&               dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      communication::syncVectorFunctionBetweenPrimitives( *src.getConformingPart(), level );
      communication::syncVectorFunctionBetweenPrimitives( *dst.getConformingPart(), level );
      src.getDiscontinuousPart()->communicate( level );
      dst.getDiscontinuousPart()->communicate( level );

      cg_to_cg_coupling_.toMatrix( mat, *src.getConformingPart(), *dst.getConformingPart(), level, flag );
      eg_to_cg_coupling_.toMatrix( mat, *src.getDiscontinuousPart(), *dst.getConformingPart(), level, flag );
      cg_to_eg_coupling_.toMatrix( mat, *src.getConformingPart(), *dst.getDiscontinuousPart(), level, flag );
      eg_to_eg_coupling_.toMatrix( mat, *src.getDiscontinuousPart(), *dst.getDiscontinuousPart(), level, flag );
   }

 private:
   P1ConstantVectorLaplaceOperator                  cg_to_cg_coupling_;
   P1ToP0ConstantP1EDGVectorLaplaceCouplingOperator cg_to_eg_coupling_;
   P0ToP1ConstantP1EDGVectorLaplaceCouplingOperator eg_to_cg_coupling_;
   P0Operator< dg::DGVectorLaplaceFormEDGEDG >      eg_to_eg_coupling_;
};
} // namespace hyteg
