/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include "hyteg/composites/P1DGEP0StokesFunction.hpp"
#include "hyteg/p0functionspace/P0P0WeightedMassForm.hpp"

#include "mixed_operator/EGOperators.hpp"

namespace hyteg {
namespace dg {
namespace eg {
template < typename VelocityBlockOperator >
class EGP0StokesPreconditionerType : public Operator< EGP0StokesFunction< real_t >, EGP0StokesFunction< real_t > >
{
 public:
   EGP0StokesPreconditionerType( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
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

   VelocityBlockOperator viscOp;
   dg::DGOperator        P;

   bool hasGlobalCells_;
};

typedef EGP0StokesPreconditionerType< EGSIPGLaplaceOperator > EGP0StokesPreconditioner;
typedef EGP0StokesPreconditionerType< EGIIPGLaplaceOperator > EGP0IIPGStokesPreconditioner;

template < typename VelocityBlockOperator, typename EnergyNormOperator, typename BlockPreconditioner >
class EGP0StokesOperatorType : public Operator< EGP0StokesFunction< real_t >, EGP0StokesFunction< real_t > >
{
 public:
   typedef VelocityBlockOperator VelocityBlockOperator_T;
   typedef EnergyNormOperator    EnergyNormOperator_T;
   typedef BlockPreconditioner   BlockPreconditioner_T;

   EGP0StokesOperatorType( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , velocityBlockOp( storage, minLevel, maxLevel )
   , div( storage, minLevel, maxLevel )
   , divT( storage, minLevel, maxLevel )
   , energyNormOp( storage, minLevel, maxLevel )
   , blockPrec( storage, minLevel, maxLevel )
   {}

   void apply( const EGP0StokesFunction< real_t >& src,
               const EGP0StokesFunction< real_t >& dst,
               const uint_t                        level,
               const DoFType                       flag,
               UpdateType                          updateType = Replace ) const override
   {
      velocityBlockOp.apply( src.uvw(), dst.uvw(), level, flag, Replace );
      divT.apply( src.p(), dst.uvw(), level, flag, Add );
      div.apply( src.uvw(), dst.p(), level, flag, Replace );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const EGP0StokesFunction< idx_t >&          src,
                  const EGP0StokesFunction< idx_t >&          dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      velocityBlockOp.toMatrix( mat, src.uvw(), dst.uvw(), level, flag );
      divT.toMatrix( mat, src.p(), dst.uvw(), level, flag );
      div.toMatrix( mat, src.uvw(), dst.p(), level, flag );
   }

   VelocityBlockOperator_T velocityBlockOp;
   EGToP0DivOperator       div;
   P0ToEGDivTOperator      divT;
   EnergyNormOperator_T    energyNormOp;
   BlockPreconditioner_T   blockPrec;
};

typedef EGP0StokesOperatorType< EGSIPGLaplaceOperator, EGLaplaceEnergyNormOperator, EGP0StokesPreconditioner > EGP0StokesOperator;
typedef EGP0StokesOperatorType< EGIIPGLaplaceOperator, EGLaplaceEnergyNormOperator, EGP0IIPGStokesPreconditioner >
    EGP0IIPGStokesOperator;

template < typename VelocityBlockOperator >
class EGP0EpsilonStokesPreconditionerType : public Operator< EGP0StokesFunction< real_t >, EGP0StokesFunction< real_t > >
{
 public:
   EGP0EpsilonStokesPreconditionerType( const std::shared_ptr< PrimitiveStorage >& storage,
                                        uint_t                                     minLevel,
                                        uint_t                                     maxLevel,
                                        std::function< real_t( const Point3D& ) >  viscosity )
   : Operator( storage, minLevel, maxLevel )
   , viscOp( storage, minLevel, maxLevel, viscosity )
   , P( storage, minLevel, maxLevel, std::make_shared< P0P0WeightedMassForm >( viscosity ) )
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

   VelocityBlockOperator              viscOp;
   P0Operator< P0P0WeightedMassForm > P;

   bool hasGlobalCells_;
};

typedef EGP0EpsilonStokesPreconditionerType< EGEpsilonOperator > EGP0EpsilonStokesPreconditioner;

template < typename VelocityBlockOperator, typename BlockPreconditioner >
class EGP0EpsilonStokesOperatorType : public Operator< EGP0StokesFunction< real_t >, EGP0StokesFunction< real_t > >
{
 public:
   typedef VelocityBlockOperator       VelocityOperator_T;
   typedef EGEpsilonEnergyNormOperator EnergyNormOperator_T;
   typedef BlockPreconditioner         BlockPreconditioner_T;

   EGP0EpsilonStokesOperatorType( const std::shared_ptr< PrimitiveStorage >& storage,
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
               const DoFType                       flag,
               UpdateType                          updateType = Replace ) const override
   {
      velocityBlockOp.apply( src.uvw(), dst.uvw(), level, flag, Replace );
      divT.apply( src.p(), dst.uvw(), level, flag, Add );
      div.apply( src.uvw(), dst.p(), level, flag, Replace );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const EGP0StokesFunction< idx_t >&          src,
                  const EGP0StokesFunction< idx_t >&          dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      velocityBlockOp.toMatrix( mat, src.uvw(), dst.uvw(), level, flag );
      divT.toMatrix( mat, src.p(), dst.uvw(), level, flag );
      div.toMatrix( mat, src.uvw(), dst.p(), level, flag );
   }

   const VelocityOperator_T& getA() const { return velocityBlockOp; }

   VelocityOperator_T    velocityBlockOp;
   EGToP0DivOperator     div;
   P0ToEGDivTOperator    divT;
   EnergyNormOperator_T  energyNormOp;
   BlockPreconditioner_T blockPrec;
};

typedef EGP0EpsilonStokesOperatorType< EGEpsilonOperator, EGP0EpsilonStokesPreconditioner > EGP0EpsilonStokesOperator;

} // namespace eg
} // namespace dg
} // namespace hyteg
