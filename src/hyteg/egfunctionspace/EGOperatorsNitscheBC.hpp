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

#include <hyteg/p0functionspace/P0P0MassForm.hpp>

#include "hyteg/dgfunctionspace/DGDivForm.hpp"
#include "hyteg/dgfunctionspace/DGVectorLaplaceForm.hpp"
#include "hyteg/dgfunctionspace/DGVectorMassForm.hpp"
#include "hyteg/dgfunctionspace/ScalarP1WithDGFormOperator.hpp"
#include "hyteg/egfunctionspace/EGDivForm.hpp"
#include "hyteg/egfunctionspace/EGDivFormNitscheBC.hpp"
#include "hyteg/egfunctionspace/EGDivtForm.hpp"
#include "hyteg/egfunctionspace/EGDivtFormNitscheBC.hpp"
#include "hyteg/egfunctionspace/EGFunction.hpp"
#include "hyteg/egfunctionspace/EGVectorLaplaceFormNitscheBC.hpp"
#include "hyteg/mixedoperators/P0ScalarToP1VectorOperator.hpp"
#include "hyteg/mixedoperators/P1ToP0Operator.hpp"
#include "hyteg/mixedoperators/P1VectorToP0ScalarOperator.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/operators/ScalarToVectorOperator.hpp"
#include "hyteg/operators/VectorLaplaceOperator.hpp"
#include "hyteg/operators/VectorMassOperator.hpp"
#include "hyteg/operators/VectorToScalarOperator.hpp"
#include "hyteg/p0functionspace/P0P0WeightedMassForm.hpp"
#include "hyteg/p0functionspace/P0Operator.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1EpsilonOperator.hpp"

namespace hyteg {
namespace dg {
namespace eg {

class EGLaplaceOperatorNitscheBC : public Operator< EGFunction< real_t >, EGFunction< real_t > >
{
 public:
   EGLaplaceOperatorNitscheBC( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : Operator< EGFunction< real_t >, EGFunction< real_t > >( storage, minLevel, maxLevel )
   , cg_( storage, minLevel, maxLevel, std::make_shared< EGVectorLaplaceFormNitscheBC_P1P1_00 >() )
   , cg_to_dg_coupling_( storage, minLevel, maxLevel )
   , dg_to_cg_coupling_( storage, minLevel, maxLevel )
   , dg_to_dg_coupling_( storage, minLevel, maxLevel, std::make_shared< EGVectorLaplaceFormNitscheBC_EE >() )
   {}

   void apply( const EGFunction< real_t >& src,
               const EGFunction< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType ) const override
   {
      dg_to_cg_coupling_.apply( *src.getDiscontinuousPart(), *dst.getConformingPart(), level, flag, updateType );
      cg_.apply( src.getConformingPart()->component( 0 ), dst.getConformingPart()->component( 0 ), level, flag, Add );
      cg_.apply( src.getConformingPart()->component( 1 ), dst.getConformingPart()->component( 1 ), level, flag, Add );
      if ( src.getDimension() == 3 )
      {
         cg_.apply( src.getConformingPart()->component( 2 ), dst.getConformingPart()->component( 2 ), level, flag, Add );
      }
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

      cg_.toMatrix( mat, src.getConformingPart()->component( 0 ), dst.getConformingPart()->component( 0 ), level, flag );
      cg_.toMatrix( mat, src.getConformingPart()->component( 1 ), dst.getConformingPart()->component( 1 ), level, flag );
      if ( src.getDimension() == 3 )
      {
         cg_.toMatrix( mat, src.getConformingPart()->component( 2 ), dst.getConformingPart()->component( 2 ), level, flag );
      }
      dg_to_cg_coupling_.toMatrix( mat, *src.getDiscontinuousPart(), *dst.getConformingPart(), level, flag );
      cg_to_dg_coupling_.toMatrix( mat, *src.getConformingPart(), *dst.getDiscontinuousPart(), level, flag );
      dg_to_dg_coupling_.toMatrix( mat, *src.getDiscontinuousPart(), *dst.getDiscontinuousPart(), level, flag );
   }

 private:
   using CGToDGOperatorType = P1VectorToP0ScalarOperator< P1ToP0Operator< EGVectorLaplaceFormNitscheBC_EP1_0 >,
                                                          P1ToP0Operator< EGVectorLaplaceFormNitscheBC_EP1_1 >,
                                                          P1ToP0Operator< EGVectorLaplaceFormNitscheBC_EP1_2 > >;
   using DGToCGOperatorType = P0ScalarToP1VectorOperator< P0ToP1Operator< EGVectorLaplaceFormNitscheBC_P1E_0 >,
                                                          P0ToP1Operator< EGVectorLaplaceFormNitscheBC_P1E_1 >,
                                                          P0ToP1Operator< EGVectorLaplaceFormNitscheBC_P1E_2 > >;
   using DGToDGOperatorType = P0Operator< EGVectorLaplaceFormNitscheBC_EE >;

   ScalarP1WithDGFormOperator cg_;
   // CGToCGOperatorType cg_to_cg_coupling_;
   CGToDGOperatorType cg_to_dg_coupling_;
   DGToCGOperatorType dg_to_cg_coupling_;
   DGToDGOperatorType dg_to_dg_coupling_;
};

// TODO: fix code duplication
class EGEpsilonOperatorNitscheBC : public Operator< EGFunction< real_t >, EGFunction< real_t > >
{
 public:
   EGEpsilonOperatorNitscheBC( const std::shared_ptr< PrimitiveStorage >& storage,
                               uint_t                                     minLevel,
                               uint_t                                     maxLevel,
                               std::function< real_t( const Point3D& ) >  viscosity )
   : Operator< EGFunction< real_t >, EGFunction< real_t > >( storage, minLevel, maxLevel )
   , viscosity_( viscosity )
   , cg_00( storage, minLevel, maxLevel, std::make_shared< eg::EGEpsilonFormNitscheBC_P1P1_00 >( viscosity ) )
   , cg_01( storage, minLevel, maxLevel, std::make_shared< eg::EGEpsilonFormNitscheBC_P1P1_01 >( viscosity ) )
   , cg_02( storage, minLevel, maxLevel, std::make_shared< eg::EGEpsilonFormNitscheBC_P1P1_02 >( viscosity ) )
   , cg_10( storage, minLevel, maxLevel, std::make_shared< eg::EGEpsilonFormNitscheBC_P1P1_10 >( viscosity ) )
   , cg_11( storage, minLevel, maxLevel, std::make_shared< eg::EGEpsilonFormNitscheBC_P1P1_11 >( viscosity ) )
   , cg_12( storage, minLevel, maxLevel, std::make_shared< eg::EGEpsilonFormNitscheBC_P1P1_12 >( viscosity ) )
   , cg_20( storage, minLevel, maxLevel, std::make_shared< eg::EGEpsilonFormNitscheBC_P1P1_20 >( viscosity ) )
   , cg_21( storage, minLevel, maxLevel, std::make_shared< eg::EGEpsilonFormNitscheBC_P1P1_21 >( viscosity ) )
   , cg_22( storage, minLevel, maxLevel, std::make_shared< eg::EGEpsilonFormNitscheBC_P1P1_22 >( viscosity ) )
   // , cg_(storage, minLevel, maxLevel, viscosity)
   , cg_to_dg_coupling_( storage,
                         minLevel,
                         maxLevel,
                         std::make_tuple( std::make_shared< typename CGToDGOperatorType::OperX_T::FormType >( viscosity ),
                                          std::make_shared< typename CGToDGOperatorType::OperY_T::FormType >( viscosity ),
                                          std::make_shared< typename CGToDGOperatorType::OperZ_T::FormType >( viscosity ) ) )
   , dg_to_cg_coupling_( storage,
                         minLevel,
                         maxLevel,
                         std::make_tuple( std::make_shared< typename DGToCGOperatorType::OperX_T::FormType >( viscosity ),
                                          std::make_shared< typename DGToCGOperatorType::OperY_T::FormType >( viscosity ),
                                          std::make_shared< typename DGToCGOperatorType::OperZ_T::FormType >( viscosity ) ) )
   , dg_to_dg_coupling_( storage, minLevel, maxLevel, std::make_shared< EGEpsilonFormNitscheBC_EE >( viscosity ) )
   {}

   void apply( const EGFunction< real_t >& src,
               const EGFunction< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType ) const override
   {
      dg_to_cg_coupling_.apply( *src.getDiscontinuousPart(), *dst.getConformingPart(), level, flag, updateType );

      cg_00.apply( src.getConformingPart()->component( 0 ), dst.getConformingPart()->component( 0 ), level, flag, Add );
      cg_01.apply( src.getConformingPart()->component( 1 ), dst.getConformingPart()->component( 0 ), level, flag, Add );
      cg_10.apply( src.getConformingPart()->component( 0 ), dst.getConformingPart()->component( 1 ), level, flag, Add );
      cg_11.apply( src.getConformingPart()->component( 1 ), dst.getConformingPart()->component( 1 ), level, flag, Add );

      if ( src.getDimension() == 3 )
      {
         WALBERLA_ABORT("not implemented.");
      }
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

      cg_00.toMatrix( mat, src.getConformingPart()->component( 0 ), dst.getConformingPart()->component( 0 ), level, flag );
      cg_01.toMatrix( mat, src.getConformingPart()->component( 1 ), dst.getConformingPart()->component( 0 ), level, flag );
      cg_10.toMatrix( mat, src.getConformingPart()->component( 0 ), dst.getConformingPart()->component( 1 ), level, flag );
      cg_11.toMatrix( mat, src.getConformingPart()->component( 1 ), dst.getConformingPart()->component( 1 ), level, flag );

      if ( src.getDimension() == 3 )
      {
         cg_02.toMatrix( mat, src.getConformingPart()->component( 2 ), dst.getConformingPart()->component( 0 ), level, flag );
         cg_12.toMatrix( mat, src.getConformingPart()->component( 2 ), dst.getConformingPart()->component( 1 ), level, flag );
         cg_20.toMatrix( mat, src.getConformingPart()->component( 0 ), dst.getConformingPart()->component( 2 ), level, flag );
         cg_21.toMatrix( mat, src.getConformingPart()->component( 1 ), dst.getConformingPart()->component( 2 ), level, flag );
         cg_22.toMatrix( mat, src.getConformingPart()->component( 2 ), dst.getConformingPart()->component( 2 ), level, flag );
      }
      //   cg_.toMatrix( mat, *src.getConformingPart(), *dst.getConformingPart(), level, flag );
      dg_to_cg_coupling_.toMatrix( mat, *src.getDiscontinuousPart(), *dst.getConformingPart(), level, flag );
      cg_to_dg_coupling_.toMatrix( mat, *src.getConformingPart(), *dst.getDiscontinuousPart(), level, flag );
      dg_to_dg_coupling_.toMatrix( mat, *src.getDiscontinuousPart(), *dst.getDiscontinuousPart(), level, flag );

   }

   using CGToDGOperatorType = P1VectorToP0ScalarOperator< P1ToP0Operator< EGEpsilonFormNitscheBC_EP1_0 >,
                                                          P1ToP0Operator< EGEpsilonFormNitscheBC_EP1_1 >,
                                                          P1ToP0Operator< EGEpsilonFormNitscheBC_EP1_2 > >;
   using DGToCGOperatorType = P0ScalarToP1VectorOperator< P0ToP1Operator< EGEpsilonFormNitscheBC_P1E_0 >,
                                                          P0ToP1Operator< EGEpsilonFormNitscheBC_P1E_1 >,
                                                          P0ToP1Operator< EGEpsilonFormNitscheBC_P1E_2 > >;
   using DGToDGOperatorType = P0Operator< EGEpsilonFormNitscheBC_EE >;

   //VectorialP1WithEpsilonFormOperator cg_;
   ScalarP1WithDGFormOperator cg_00;
   ScalarP1WithDGFormOperator cg_01;
   ScalarP1WithDGFormOperator cg_02;
   ScalarP1WithDGFormOperator cg_10;
   ScalarP1WithDGFormOperator cg_11;
   ScalarP1WithDGFormOperator cg_12;
   ScalarP1WithDGFormOperator cg_20;
   ScalarP1WithDGFormOperator cg_21;
   ScalarP1WithDGFormOperator cg_22;

   CGToDGOperatorType                        cg_to_dg_coupling_;
   DGToCGOperatorType                        dg_to_cg_coupling_;
   DGToDGOperatorType                        dg_to_dg_coupling_;
   std::function< real_t( const Point3D& ) > viscosity_;
};

// TODO: fix code duplication
class EGToP0DivOperatorNitscheBC final : public Operator< EGFunction< real_t >, P0Function< real_t > >
{
 public:
   EGToP0DivOperatorNitscheBC( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : Operator< EGFunction< real_t >, P0Function< real_t > >( storage, minLevel, maxLevel )
   , p1x_to_p0( storage, minLevel, maxLevel )
   , p1y_to_p0( storage, minLevel, maxLevel )
   , p1z_to_p0( storage, minLevel, maxLevel )

   // ,  p1_to_p0( storage, minLevel, maxLevel )
   , edg_to_p0( storage, minLevel, maxLevel )
   {}

   void apply( const EGFunction< real_t >& src,
               const P0Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType ) const override
   {
      communication::syncVectorFunctionBetweenPrimitives( *src.getConformingPart(), level );
      dst.communicate( level );
      src.getDiscontinuousPart()->communicate( level );

      p1x_to_p0.apply( src.getConformingPart()->component( 0 ), dst, level, flag, updateType );
      p1y_to_p0.apply( src.getConformingPart()->component( 1 ), dst, level, flag, Add );
      if ( src.getDimension() == 3 )
      {
         p1z_to_p0.apply( src.getConformingPart()->component( 2 ), dst, level, flag, Add );
      }

      // p1_to_p0.apply( *src.getConformingPart(), dst, level, flag, updateType );

      edg_to_p0.apply( *src.getDiscontinuousPart(), dst, level, flag, Add );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const EGFunction< idx_t >&                  src,
                  const P0Function< idx_t >&                  dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      communication::syncVectorFunctionBetweenPrimitives( *src.getConformingPart(), level );
      dst.communicate( level );
      src.getDiscontinuousPart()->communicate( level );

      p1x_to_p0.toMatrix( mat, src.getConformingPart()->component( 0 ), dst, level, flag );
      p1y_to_p0.toMatrix( mat, src.getConformingPart()->component( 1 ), dst, level, flag );
      if ( src.getDimension() == 3 )
      {
         p1z_to_p0.toMatrix( mat, src.getConformingPart()->component( 2 ), dst, level, flag );
      }

      //    p1_to_p0.toMatrix( mat, *src.getConformingPart(), dst, level, flag );
      edg_to_p0.toMatrix( mat, *src.getDiscontinuousPart(), dst, level, flag );
   }

 private:
   /*
    P1ToP0Operator< EGDivFormP0P1_new_0 > p1x_to_p0;
    P1ToP0Operator< EGDivFormP0P1_new_1 > p1y_to_p0;
    P1ToP0Operator< dg::DGFormAbort > p1z_to_p0;
    //P1ToP0ConstantDivOperator p1_to_p0;
    P0Operator< EGDivFormP0EDG_new > edg_to_p0;
*/

   P1ToP0Operator< dg::eg::EGDivFormNitscheBC_P0P1_0 > p1x_to_p0;
   P1ToP0Operator< dg::eg::EGDivFormNitscheBC_P0P1_1 > p1y_to_p0;
   P1ToP0Operator< dg::eg::EGDivFormNitscheBC_P0P1_2 > p1z_to_p0;

   //P1ToP0ConstantDivOperator p1_to_p0;
   P0Operator< dg::eg::EGDivFormNitscheBC_P0E > edg_to_p0;
};

// TODO: fix code duplication
class P0ToEGDivTOperatorNitscheBC final : public Operator< P0Function< real_t >, EGFunction< real_t > >
{
 public:
   P0ToEGDivTOperatorNitscheBC( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : Operator< P0Function< real_t >, EGFunction< real_t > >( storage, minLevel, maxLevel )

   , p0_to_p1x( storage, minLevel, maxLevel )
   , p0_to_p1y( storage, minLevel, maxLevel )
   , p0_to_p1z( storage, minLevel, maxLevel )
   // ,   p0_to_p1( storage, minLevel, maxLevel )
   , p0_to_edg( storage, minLevel, maxLevel )
   {}

   void apply( const P0Function< real_t >& src,
               const EGFunction< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType ) const override
   {
      //p0_to_p1.apply( src, *dst.getConformingPart(), level, flag, updateType );

      p0_to_p1x.apply( src, dst.getConformingPart()->component( 0 ), level, flag, updateType );
      p0_to_p1y.apply( src, dst.getConformingPart()->component( 1 ), level, flag, updateType );
      if ( src.getDimension() == 3 )
      {
         p0_to_p1z.apply( src, dst.getConformingPart()->component( 2 ), level, flag, updateType );
      }

      p0_to_edg.apply( src, *dst.getDiscontinuousPart(), level, flag, updateType );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P0Function< idx_t >&                  src,
                  const EGFunction< idx_t >&                  dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      communication::syncVectorFunctionBetweenPrimitives( *dst.getConformingPart(), level );
      dst.getDiscontinuousPart()->communicate( level );
      src.communicate( level );

      p0_to_p1x.toMatrix( mat, src, dst.getConformingPart()->component( 0 ), level, flag );
      p0_to_p1y.toMatrix( mat, src, dst.getConformingPart()->component( 1 ), level, flag );
      if ( src.getDimension() == 3 )
      {
         p0_to_p1z.toMatrix( mat, src, dst.getConformingPart()->component( 2 ), level, flag );
      }

      // p0_to_p1.toMatrix( mat, src, *dst.getConformingPart(), level, flag );
      p0_to_edg.toMatrix( mat, src, *dst.getDiscontinuousPart(), level, flag );
   }

 private:
   /*
    P0ToP1Operator< EGDivtFormP1P0_new_0 > p0_to_p1x;
    P0ToP1Operator< EGDivtFormP1P0_new_1 > p0_to_p1y;
    P0ToP1Operator< dg::DGFormAbort > p0_to_p1z;

    // P0ToP1ConstantDivTOperator p0_to_p1;
    P0Operator< EGDivtFormEDGP0_new > p0_to_edg;
*/

   P0ToP1Operator< dg::eg::EGDivtFormNitscheBC_P1P0_0 > p0_to_p1x;
   P0ToP1Operator< dg::eg::EGDivtFormNitscheBC_P1P0_1 > p0_to_p1y;
   P0ToP1Operator< dg::eg::EGDivtFormNitscheBC_P1P0_2 > p0_to_p1z;

   // P0ToP1ConstantDivTOperator p0_to_p1;
   P0Operator< dg::eg::EGDivtFormNitscheBC_EP0 > p0_to_edg;
};

class EGP0StokesPreconditionerNitscheBC : public Operator< EGP0StokesFunction< real_t >, EGP0StokesFunction< real_t > >
{
 public:
   EGP0StokesPreconditionerNitscheBC( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , viscOp( storage, minLevel, maxLevel )
   , P( storage, minLevel, maxLevel, std::make_shared< dg::P0P0MassForm >() )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const EGP0StokesFunction< idx_t >&          src,
                  const EGP0StokesFunction< idx_t >&          dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      viscOp.toMatrix( mat, src.uvw(), dst.uvw(), level, flag );
      P.toMatrix( mat, *src.p().getDGFunction(), *dst.p().getDGFunction(), level, flag );
   }

   eg::EGLaplaceOperatorNitscheBC viscOp;
   dg::DGOperator                 P;

   bool hasGlobalCells_;
};

class EGP0StokesOperatorNitscheBC : public Operator< EGP0StokesFunction< real_t >, EGP0StokesFunction< real_t > >
{
 public:
   typedef eg::EGLaplaceEnergyNormOperator   EnergyNormOperator_T;
   typedef EGP0StokesPreconditionerNitscheBC BlockPreconditioner_T;

   EGP0StokesOperatorNitscheBC( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , velocityBlockOp( storage, minLevel, maxLevel )
   , div( storage, minLevel, maxLevel )
   , divT( storage, minLevel, maxLevel )
   , energyNormOp( storage, minLevel, maxLevel )
   ,blockPrec( storage, minLevel, maxLevel )
   {}

   void apply( const EGP0StokesFunction< real_t >& src,
               const EGP0StokesFunction< real_t >& dst,
               const uint_t                        level,
               const DoFType                       flag ) const
   {
      velocityBlockOp.apply( src.uvw(), dst.uvw(), level, flag, Replace );
      divT.apply( src.p(), dst.uvw(), level, flag, Add );
      div.apply( src.uvw(), dst.p(), level, flag, Replace );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const EGP0StokesFunction< idx_t >&          src,
                  const EGP0StokesFunction< idx_t >&          dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      velocityBlockOp.toMatrix( mat, src.uvw(), dst.uvw(), level, flag );
      divT.toMatrix( mat, src.p(), dst.uvw(), level, flag );
      div.toMatrix( mat, src.uvw(), dst.p(), level, flag );
   }

   eg::EGLaplaceOperatorNitscheBC velocityBlockOp;
   EnergyNormOperator_T           energyNormOp;
   EGToP0DivOperatorNitscheBC     div;
   P0ToEGDivTOperatorNitscheBC    divT;
   BlockPreconditioner_T blockPrec;
};

class EGP0EpsilonStokesPreconditionerNitscheBC
: public Operator< EGP0StokesFunction< real_t >, EGP0StokesFunction< real_t > >
{
 public:
   EGP0EpsilonStokesPreconditionerNitscheBC( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel,
                                        std::function<real_t(const Point3D &)> viscosity )
   : Operator( storage, minLevel, maxLevel )
   , viscOp( storage, minLevel, maxLevel, viscosity )
   , P( storage, minLevel, maxLevel, std::make_shared< P0P0WeightedMassForm >(viscosity) )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const EGP0StokesFunction< idx_t >&          src,
                  const EGP0StokesFunction< idx_t >&          dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      viscOp.toMatrix( mat, src.uvw(), dst.uvw(), level, flag );
      P.toMatrix( mat, src.p(), dst.p(), level, flag );
   }
   eg::EGEpsilonOperatorNitscheBC viscOp;
   P0Operator< P0P0WeightedMassForm>               P;

   bool hasGlobalCells_;
};

class EGP0EpsilonOperatorStokesNitscheBC : public Operator< EGP0StokesFunction< real_t >, EGP0StokesFunction< real_t > >
{
 public:

   typedef  eg::EGEpsilonOperatorNitscheBC       VelocityOperator_T;
   typedef eg::EGEpsilonEnergyNormOperator EnergyNormOperator_T;
   typedef EGP0EpsilonStokesPreconditionerNitscheBC BlockPreconditioner_T;

   EGP0EpsilonOperatorStokesNitscheBC( const std::shared_ptr< PrimitiveStorage >& storage,
                                       size_t                                     minLevel,
                                       size_t                                     maxLevel,
                                       std::function< real_t( const Point3D& ) >  viscosity )
   : Operator( storage, minLevel, maxLevel )
   , velocityBlockOp( storage, minLevel, maxLevel, viscosity )
   , div( storage, minLevel, maxLevel )
   , divT( storage, minLevel, maxLevel )
   , energyNormOp( storage, minLevel, maxLevel, viscosity )
   , blockPrec( storage, minLevel, maxLevel, viscosity )
   {}

   void apply( const EGP0StokesFunction< real_t >& src,
               const EGP0StokesFunction< real_t >& dst,
               const uint_t                        level,
               const DoFType                       flag ) const
   {
      velocityBlockOp.apply( src.uvw(), dst.uvw(), level, flag, Replace );
      divT.apply( src.p(), dst.uvw(), level, flag, Add );
      div.apply( src.uvw(), dst.p(), level, flag, Replace );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const EGP0StokesFunction< idx_t >&          src,
                  const EGP0StokesFunction< idx_t >&          dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      velocityBlockOp.toMatrix( mat, src.uvw(), dst.uvw(), level, flag );
      divT.toMatrix( mat, src.p(), dst.uvw(), level, flag );
      div.toMatrix( mat, src.uvw(), dst.p(), level, flag );
   }

   const  eg::EGEpsilonOperatorNitscheBC& getA() const { return velocityBlockOp; }

   eg::EGEpsilonOperatorNitscheBC velocityBlockOp;
   EnergyNormOperator_T           energyNormOp;
   EGToP0DivOperatorNitscheBC     div;
   P0ToEGDivTOperatorNitscheBC    divT;
   BlockPreconditioner_T blockPrec;
};

} // namespace eg
} // namespace dg
} // namespace hyteg