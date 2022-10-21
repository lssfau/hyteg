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

#include "core/Abort.h"
#include "core/DataTypes.h"

#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/gridtransferoperators/RestrictionOperator.hpp"
#include "hyteg/types/types.hpp"

namespace hyteg {

class RestrictionForm
{
 public:
   using Point = Eigen::Matrix< real_t, 3, 1 >;

   virtual void integrate2D( const std::vector< Point >&                              dst,
                             const std::vector< Point >&                              src,
                             Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& localMat ) const
   {
      WALBERLA_UNUSED( dst );
      WALBERLA_UNUSED( src );
      WALBERLA_UNUSED( localMat );
      WALBERLA_LOG_INFO_ON_ROOT( "not implemented" )
   }
};

class RestrictionFormDG1 : public RestrictionForm
{
 public:
   void integrate2D( const std::vector< Point >&                              dst,
                     const std::vector< Point >&                              src,
                     Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& localMat ) const override;
};

class RestrictionFormDG0 : public RestrictionForm
{
 public:
   void integrate2D( const std::vector< Point >&                              dst,
                     const std::vector< Point >&                              src,
                     Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& localMat ) const override
   {
      WALBERLA_UNUSED( dst );
      WALBERLA_UNUSED( src );
      localMat.resize( 1, 1 );
      localMat( 0, 0 ) = 1;
   }
};

class DGRestriction : public RestrictionOperator< dg::DGFunction< real_t > >
{
 public:
   explicit DGRestriction( std::shared_ptr< RestrictionForm > form )
   : form_( form )
   {}

   void restrict( const dg::DGFunction< real_t >& function,
                  const walberla::uint_t&         fineLevel,
                  const DoFType&                  flag ) const override;

 protected:
   std::shared_ptr< RestrictionForm > form_;
};

class DG1Restriction : public DGRestriction
{
 public:
   DG1Restriction()
   : DGRestriction( std::make_shared< RestrictionFormDG1 >() )
   {}
};

class DG0Restriction : public DGRestriction
{
 public:
   DG0Restriction()
   : DGRestriction( std::make_shared< RestrictionFormDG0 >() )
   {}
};

} // namespace hyteg
