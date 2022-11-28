/*
* Copyright (c) 2022 Andreas Wagner.
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
#include "hyteg/dgfunctionspace/P1WithDGFormOperator.hpp"
#include "hyteg/egfunctionspace/EGDivForm.hpp"
#include "hyteg/egfunctionspace/EGDivtForm.hpp"
#include "hyteg/egfunctionspace/EGFunction.hpp"
#include "hyteg/egfunctionspace/EGVectorLaplaceFormNew.hpp"
#include "hyteg/mixedoperators/P0ScalarToP1VectorOperator.hpp"
#include "hyteg/mixedoperators/P1ToP0Operator.hpp"
#include "hyteg/mixedoperators/P1VectorToP0ScalarOperator.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/operators/ScalarToVectorOperator.hpp"
#include "hyteg/operators/VectorLaplaceOperator.hpp"
#include "hyteg/operators/VectorMassOperator.hpp"
#include "hyteg/operators/VectorToScalarOperator.hpp"
#include "hyteg/p0functionspace/P0Operator.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1EpsilonOperator.hpp"

namespace hyteg {
namespace dg {
namespace eg {

class EGLaplaceOperatorNew : public Operator< EGFunction< real_t >, EGFunction< real_t > >
{
 public:
   EGLaplaceOperatorNew( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : Operator< EGFunction< real_t >, EGFunction< real_t > >( storage, minLevel, maxLevel )
   // , cg_to_cg_coupling_( storage, minLevel, maxLevel, std::make_shared< EGVectorLaplaceFormP1P1_new_00 >() )
   , cg_(storage, minLevel, maxLevel, std::make_shared< EGVectorLaplaceFormP1P1_new_00 >() )
   // , cg_to_cg_coupling_( storage, minLevel, maxLevel )
   , cg_to_dg_coupling_( storage, minLevel, maxLevel )
   , dg_to_cg_coupling_( storage, minLevel, maxLevel )
   , dg_to_dg_coupling_( storage, minLevel, maxLevel, std::make_shared< EGVectorLaplaceFormEDGEDG_new >() )
   {
      // auto form = std::make_shared< EGVectorLaplaceFormP1P1_new_00 >();
      // auto op   = std::make_shared< P1WithDGFormOperator >( storage, minLevel, maxLevel, form );
      // cg_to_cg_coupling_.setSubOperator( 0, 0, op );
      // cg_to_cg_coupling_.setSubOperator( 1, 1, op );
   }

   void apply( const EGFunction< real_t >& src,
               const EGFunction< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType ) const override
   {
      dg_to_cg_coupling_.apply( *src.getDiscontinuousPart(), *dst.getConformingPart(), level, flag, updateType );
      cg_.apply( src.getConformingPart()->component(0), dst.getConformingPart()->component(0), level, flag, Add );
      cg_.apply( src.getConformingPart()->component(1), dst.getConformingPart()->component(1), level, flag, Add );
      // cg_to_cg_coupling_.apply( *src.getConformingPart(), *dst.getConformingPart(), level, flag, Add );

      cg_to_dg_coupling_.apply( *src.getConformingPart(), *dst.getDiscontinuousPart(), level, flag, updateType );
      dg_to_dg_coupling_.apply( *src.getDiscontinuousPart(), *dst.getDiscontinuousPart(), level, flag, Add );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const EGFunction< idx_t >&                  src,
                  const EGFunction< idx_t >&                  dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      communication::syncVectorFunctionBetweenPrimitives( *src.getConformingPart(), level );
      communication::syncVectorFunctionBetweenPrimitives( *dst.getConformingPart(), level );
      src.getDiscontinuousPart()->communicate( level );
      dst.getDiscontinuousPart()->communicate( level );

      // cg_to_cg_coupling_.toMatrix( mat, *src.getConformingPart(), *dst.getConformingPart(), level, flag );
      // cg_to_cg_coupling_.toMatrix( mat, *src.getConformingPart(), *dst.getConformingPart(), level, flag );
      cg_.toMatrix( mat, src.getConformingPart()->component(0), dst.getConformingPart()->component(0), level, flag );
      cg_.toMatrix( mat, src.getConformingPart()->component(1), dst.getConformingPart()->component(1), level, flag );
      dg_to_cg_coupling_.toMatrix( mat, *src.getDiscontinuousPart(), *dst.getConformingPart(), level, flag );
      cg_to_dg_coupling_.toMatrix( mat, *src.getConformingPart(), *dst.getDiscontinuousPart(), level, flag );
      dg_to_dg_coupling_.toMatrix( mat, *src.getDiscontinuousPart(), *dst.getDiscontinuousPart(), level, flag );
   }

 private:
   // using CGToCGOperatorType = VectorLaplaceOperator< real_t, P1VectorFunction, P1WithDGFormOperator >;
   // using CGToCGOperatorType = VectorToVectorOperator< real_t, P1VectorFunction, P1VectorFunction >;
   using CGToDGOperatorType = P1VectorToP0ScalarOperator< P1ToP0Operator< EGVectorLaplaceFormEDGP1_new_0 >,
                                                          P1ToP0Operator< EGVectorLaplaceFormEDGP1_new_1 >,
                                                          P1ToP0Operator< DGFormAbort > >;
   using DGToCGOperatorType = P0ScalarToP1VectorOperator< P0ToP1Operator< EGVectorLaplaceFormP1EDG_new_0 >,
                                                          P0ToP1Operator< EGVectorLaplaceFormP1EDG_new_1 >,
                                                          P0ToP1Operator< DGFormAbort > >;
   using DGToDGOperatorType = P0Operator< EGVectorLaplaceFormEDGEDG_new >;

   P1WithDGFormOperator cg_;
   // CGToCGOperatorType cg_to_cg_coupling_;
   CGToDGOperatorType cg_to_dg_coupling_;
   DGToCGOperatorType dg_to_cg_coupling_;
   DGToDGOperatorType dg_to_dg_coupling_;
};

} // namespace eg
} // namespace dg
} // namespace hyteg