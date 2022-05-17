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

#include "hyteg/dgfunctionspace/DGDivForm.hpp"
#include "hyteg/dgfunctionspace/DGVectorLaplaceForm.hpp"
#include "hyteg/dgfunctionspace/DGVectorMassForm.hpp"
#include "hyteg/mixedoperators/P0ScalarToP1VectorOperator.hpp"
#include "hyteg/mixedoperators/P1VectorToP0ScalarOperator.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/operators/ScalarToVectorOperator.hpp"
#include "hyteg/operators/VectorLaplaceOperator.hpp"
#include "hyteg/operators/VectorMassOperator.hpp"
#include "hyteg/operators/VectorToScalarOperator.hpp"
#include "hyteg/p0functionspace/P0Operator.hpp"
#include "hyteg/p1dgefunctionspace/P1DGEFunction.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"

namespace hyteg {

class P1ToP1DGEDivTOperator final : public Operator< P1Function< real_t >, P1DGEFunction< real_t > >
{
 public:
   P1ToP1DGEDivTOperator( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : Operator< P1Function< real_t >, P1DGEFunction< real_t > >( storage, minLevel, maxLevel )
   , cg_to_cg_coupling_( storage, minLevel, maxLevel )
   , cg_to_eg_coupling_( storage, minLevel, maxLevel )
   {}

   void apply( const P1Function< real_t >&    src,
               const P1DGEFunction< real_t >& dst,
               size_t                         level,
               DoFType                        flag,
               UpdateType                     updateType ) const override
   {
      cg_to_eg_coupling_.apply( src, *dst.getDiscontinuousPart(), level, flag, updateType );
      cg_to_cg_coupling_.apply( src, *dst.getConformingPart(), level, flag, Add );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P1Function< idx_t >&                  src,
                  const P1DGEFunction< idx_t >&               dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      communication::syncFunctionBetweenPrimitives( src, level );
      communication::syncVectorFunctionBetweenPrimitives( dst, level );

      cg_to_cg_coupling_.toMatrix( mat, src, *dst.getConformingPart(), level, flag );
      cg_to_eg_coupling_.toMatrix( mat, src, *dst.getDiscontinuousPart(), level, flag );
   }

 private:
   P1ConstantDivTOperator                         cg_to_cg_coupling_;
   P1ToP0ConstantP1EDGVDivergenceCouplingOperator cg_to_eg_coupling_;
};

class P1DGEToP1DivOperator final : public Operator< P1DGEFunction< real_t >, P1Function< real_t > >
{
 public:
   P1DGEToP1DivOperator( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : Operator< P1DGEFunction< real_t >, P1Function< real_t > >( storage, minLevel, maxLevel )
   , cg_to_cg_coupling_( storage, minLevel, maxLevel )
   , eg_to_cg_coupling_( storage, minLevel, maxLevel )
   {}

   void apply( const P1DGEFunction< real_t >& src,
               const P1Function< real_t >&    dst,
               size_t                         level,
               DoFType                        flag,
               UpdateType                     updateType ) const override
   {
      eg_to_cg_coupling_.apply( *src.getDiscontinuousPart(), dst, level, flag, updateType );
      cg_to_cg_coupling_.apply( *src.getConformingPart(), dst, level, flag, Add );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P1DGEFunction< idx_t >&               src,
                  const P1Function< idx_t >&                  dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      communication::syncVectorFunctionBetweenPrimitives( *src.getConformingPart(), level );
      communication::syncFunctionBetweenPrimitives( dst, level );
      src.getDiscontinuousPart()->communicate( level );

      cg_to_cg_coupling_.toMatrix( mat, *src.getConformingPart(), dst, level, flag );
      eg_to_cg_coupling_.toMatrix( mat, *src.getDiscontinuousPart(), dst, level, flag );
   }

 private:
   P1ConstantDivOperator                         cg_to_cg_coupling_;
   P0ToP1ConstantP1EDGDivergenceCouplingOperator eg_to_cg_coupling_;
};

class P0ToP1DGEDivTOperator final : public Operator< P0Function< real_t >, P1DGEFunction< real_t > >
{
 public:
   P0ToP1DGEDivTOperator( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : Operator< P0Function< real_t >, P1DGEFunction< real_t > >( storage, minLevel, maxLevel )
   , p0_to_p1x( storage, minLevel, maxLevel )
   , p0_to_p1y( storage, minLevel, maxLevel )
   , p0_to_edg( storage, minLevel, maxLevel )
   {}

   void apply( const P0Function< real_t >&    src,
               const P1DGEFunction< real_t >& dst,
               size_t                         level,
               DoFType                        flag,
               UpdateType                     updateType ) const override
   {
      p0_to_p1x.apply( src, dst.getConformingPart()->component( 0 ), level, flag, updateType );
      p0_to_p1y.apply( src, dst.getConformingPart()->component( 1 ), level, flag, updateType );
      p0_to_edg.apply( src, *dst.getDiscontinuousPart(), level, flag, updateType );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P0Function< idx_t >&                  src,
                  const P1DGEFunction< idx_t >&               dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      communication::syncVectorFunctionBetweenPrimitives( dst, level );
      src.communicate( level );

      p0_to_p1x.toMatrix( mat, src, dst.getConformingPart()->component( 0 ), level, flag );
      p0_to_p1y.toMatrix( mat, src, dst.getConformingPart()->component( 1 ), level, flag );
      p0_to_edg.toMatrix( mat, src, *dst.getDiscontinuousPart(), level, flag );
   }

 private:
   P0ToP1Operator< dg::DGDivtFormP1P0_0 > p0_to_p1x;
   P0ToP1Operator< dg::DGDivtFormP1P0_1 > p0_to_p1y;
   P0Operator< dg::DGDivtFormEDGP0 >      p0_to_edg;
};

class P1DGEToP0DivOperator final : public Operator< P1DGEFunction< real_t >, P0Function< real_t > >
{
 public:
   P1DGEToP0DivOperator( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : Operator< P1DGEFunction< real_t >, P0Function< real_t > >( storage, minLevel, maxLevel )
   , p1x_to_p0( storage, minLevel, maxLevel )
   , p1y_to_p0( storage, minLevel, maxLevel )
   , p0_to_edg( storage, minLevel, maxLevel )
   {}

   void apply( const P1DGEFunction< real_t >& src,
               const P0Function< real_t >&    dst,
               size_t                         level,
               DoFType                        flag,
               UpdateType                     updateType ) const override
   {
      communication::syncVectorFunctionBetweenPrimitives( *src.getConformingPart(), level );
      dst.communicate( level );
      src.getDiscontinuousPart()->communicate( level );

      p1x_to_p0.apply( src.getConformingPart()->component( 0 ), dst, level, flag, updateType );
      p1y_to_p0.apply( src.getConformingPart()->component( 1 ), dst, level, flag, Add );
      p0_to_edg.apply( *src.getDiscontinuousPart(), dst, level, flag, Add );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P1DGEFunction< idx_t >&               src,
                  const P0Function< idx_t >&                  dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      communication::syncVectorFunctionBetweenPrimitives( *src.getConformingPart(), level );
      dst.communicate( level );
      src.getDiscontinuousPart()->communicate( level );

      p1x_to_p0.toMatrix( mat, src.getConformingPart()->component( 0 ), dst, level, flag );
      p1y_to_p0.toMatrix( mat, src.getConformingPart()->component( 1 ), dst, level, flag );
      p0_to_edg.toMatrix( mat, *src.getDiscontinuousPart(), dst, level, flag );
   }

 private:
   P1ToP0Operator< dg::DGDivFormP0P1_0 > p1x_to_p0;
   P1ToP0Operator< dg::DGDivFormP0P1_1 > p1y_to_p0;
   P0Operator< dg::DGDivFormP0EDG >      p0_to_edg;
};

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
