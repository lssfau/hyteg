/*
 * Copyright (c) 2017-2020 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/elementwiseoperators/P1ToP2ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesBlockPreconditioner.hpp"
#include "hyteg/elementwiseoperators/P2ToP1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseBlendingFullViscousOperator.hpp"

#include "mixed_operator/P2P1TaylorHoodStokesBlockPreconditioner.hpp"
#include "mixed_operator/ScalarToVectorOperator.hpp"
#include "mixed_operator/VectorLaplaceOperator.hpp"
#include "mixed_operator/VectorToScalarOperator.hpp"

namespace hyteg {

class P2P1ElementwiseBlendingStokesOperator
: public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:
   typedef P2ElementwiseBlendingLaplaceOperator VelocityOperator_T;

   // This typedef could be a temporary solution to allow for using
   // the PETScBlockPreconditionedStokesSolver and for compiling
   // the StrongFreeSlipWrapper for the P2P1ElementwiseBlendingStokesOperator
   // typedef P2P1TaylorHoodStokesBlockPreconditioner BlockPreconditioner_T;

   typedef P2P1ElementwiseBlendingStokesBlockPreconditioner BlockPreconditioner_T;

   P2P1ElementwiseBlendingStokesOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , lapl( storage, minLevel, maxLevel )
   , div( storage, minLevel, maxLevel )
   , divT( storage, minLevel, maxLevel )
   , pspg_inv_diag_( storage, minLevel, maxLevel )
   , blockPrec( storage, minLevel, maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   void computeAndStoreLocalElementMatrices()
   {
      auto scalarA = dynamic_cast< P2ElementwiseBlendingLaplaceOperator& >( *lapl.getSubOperator( 0, 0 ) );
      scalarA.computeAndStoreLocalElementMatrices();

      div.getSubOperator< 0 >().computeAndStoreLocalElementMatrices();
      div.getSubOperator< 1 >().computeAndStoreLocalElementMatrices();

      divT.getSubOperator< 0 >().computeAndStoreLocalElementMatrices();
      divT.getSubOperator< 1 >().computeAndStoreLocalElementMatrices();

      div.getSubOperator< 0 >().computeAndStoreLocalElementMatrices();
      div.getSubOperator< 1 >().computeAndStoreLocalElementMatrices();

      divT.getSubOperator< 0 >().computeAndStoreLocalElementMatrices();
      divT.getSubOperator< 1 >().computeAndStoreLocalElementMatrices();

      if ( hasGlobalCells_ )
      {
         div.getSubOperator< 2 >().computeAndStoreLocalElementMatrices();
         divT.getSubOperator< 2 >().computeAndStoreLocalElementMatrices();

         div.getSubOperator< 2 >().computeAndStoreLocalElementMatrices();
         divT.getSubOperator< 2 >().computeAndStoreLocalElementMatrices();
      }
   }

   void apply( const P2P1TaylorHoodFunction< real_t >& src,
               const P2P1TaylorHoodFunction< real_t >& dst,
               const uint_t                            level,
               const DoFType                           flag,
               UpdateType                              updateType = Replace ) const override
   {
      lapl.apply( src.uvw(), dst.uvw(), level, flag, Replace );
      divT.apply( src.p(), dst.uvw(), level, flag, Add );
      div.apply( src.uvw(), dst.p(), level, flag, Replace );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2P1TaylorHoodFunction< idx_t >&      src,
                  const P2P1TaylorHoodFunction< idx_t >&      dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      lapl.toMatrix( mat, src.uvw(), dst.uvw(), level, flag );
      divT.toMatrix( mat, src.p(), dst.uvw(), level, flag );
      div.toMatrix( mat, src.uvw(), dst.p(), level, flag );
   }

   const P2ElementwiseBlendingLaplaceOperator& getA() const
   {
      auto ptr = lapl.getSubOperator( 0, 0 );
      return dynamic_cast< const P2ElementwiseBlendingLaplaceOperator& >( *ptr );
   }

   P2ElementwiseBlendingVectorLaplaceOperator lapl;
   P2ToP1ElementwiseBlendingDivOperator       div;
   P1ToP2ElementwiseBlendingDivTOperator      divT;

   /// this operator is need in the uzawa smoother
   P1PSPGInvDiagOperator pspg_inv_diag_;

    BlockPreconditioner_T blockPrec;
   bool hasGlobalCells_;
};

class P2P1ElementwiseBlendingFullViscousStokesOperator
: public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:
   typedef P2ElementwiseBlendingFullViscousOperator VelocityOperator_T;

   // This typedef could be a temporary solution to allow for using
   // the PETScBlockPreconditionedStokesSolver and for compiling
   // the StrongFreeSlipWrapper for the P2P1ElementwiseBlendingStokesOperator
   // typedef P2P1TaylorHoodStokesBlockPreconditioner BlockPreconditioner_T;

   typedef P2P1ElementwiseBlendingStokesBlockPreconditioner BlockPreconditioner_T;

   P2P1ElementwiseBlendingFullViscousStokesOperator(
       const std::shared_ptr< PrimitiveStorage >& storage,
       size_t                                     minLevel,
       size_t                                     maxLevel,
       std::function< real_t( const Point3D& ) >  viscosity = []( const Point3D& ) { return 1.0; } )
   : Operator( storage, minLevel, maxLevel )
   , lapl( storage, minLevel, maxLevel, viscosity )
   , div( storage, minLevel, maxLevel )
   , divT( storage, minLevel, maxLevel )
   , pspg_inv_diag_( storage, minLevel, maxLevel )
   , blockPrec( storage, minLevel, maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   void computeAndStoreLocalElementMatrices()
   {
      auto scalarA = dynamic_cast< P2ElementwiseBlendingLaplaceOperator& >( *lapl.getSubOperator( 0, 0 ) );
      scalarA.computeAndStoreLocalElementMatrices();

      div.getSubOperator< 0 >().computeAndStoreLocalElementMatrices();
      div.getSubOperator< 1 >().computeAndStoreLocalElementMatrices();

      divT.getSubOperator< 0 >().computeAndStoreLocalElementMatrices();
      divT.getSubOperator< 1 >().computeAndStoreLocalElementMatrices();

      div.getSubOperator< 0 >().computeAndStoreLocalElementMatrices();
      div.getSubOperator< 1 >().computeAndStoreLocalElementMatrices();

      divT.getSubOperator< 0 >().computeAndStoreLocalElementMatrices();
      divT.getSubOperator< 1 >().computeAndStoreLocalElementMatrices();

      if ( hasGlobalCells_ )
      {
         div.getSubOperator< 2 >().computeAndStoreLocalElementMatrices();
         divT.getSubOperator< 2 >().computeAndStoreLocalElementMatrices();

         div.getSubOperator< 2 >().computeAndStoreLocalElementMatrices();
         divT.getSubOperator< 2 >().computeAndStoreLocalElementMatrices();
      }
   }

   void apply( const P2P1TaylorHoodFunction< real_t >& src,
               const P2P1TaylorHoodFunction< real_t >& dst,
               const uint_t                            level,
               const DoFType                           flag,
               UpdateType                              updateType = Replace ) const override
   {
      lapl.apply( src.uvw(), dst.uvw(), level, flag, Replace );
      divT.apply( src.p(), dst.uvw(), level, flag, Add );
      div.apply( src.uvw(), dst.p(), level, flag, Replace );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2P1TaylorHoodFunction< idx_t >&      src,
                  const P2P1TaylorHoodFunction< idx_t >&      dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      lapl.toMatrix( mat, src.uvw(), dst.uvw(), level, flag );
      divT.toMatrix( mat, src.p(), dst.uvw(), level, flag );
      div.toMatrix( mat, src.uvw(), dst.p(), level, flag );
   }

   const VelocityOperator_T::visc_0_0& getA() const
   {
      auto ptr = lapl.getSubOperator( 0, 0 );
      return dynamic_cast< const VelocityOperator_T::visc_0_0& >( *ptr );
   }

   VelocityOperator_T                    lapl;
   P2ToP1ElementwiseBlendingDivOperator  div;
   P1ToP2ElementwiseBlendingDivTOperator divT;

   /// this operator is need in the uzawa smoother
   P1PSPGInvDiagOperator pspg_inv_diag_;

   BlockPreconditioner_T blockPrec;
   bool                  hasGlobalCells_;
};

} // namespace hyteg
