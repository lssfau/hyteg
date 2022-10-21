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
#include "hyteg/gridtransferoperators/ProlongationOperator.hpp"
#include "hyteg/types/types.hpp"

namespace hyteg {

class ProlongationForm
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

class ProlongationFormDG1 : public ProlongationForm
{
 public:
   void integrate2D( const std::vector< Point >&                              dst,
                     const std::vector< Point >&                              src,
                     Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& localMat ) const override;
};

class ProlongationFormDG0 : public ProlongationForm
{
 public:
   using Point3 = Eigen::Matrix< real_t, 3, 1 >;

   void integrate2D( const std::vector< Point3 >&                             dst,
                     const std::vector< Point3 >&                             src,
                     Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& localMat ) const override
   {
      WALBERLA_UNUSED( dst );
      WALBERLA_UNUSED( src );
      localMat.resize( 1, 1 );
      localMat( 0, 0 ) = 1;
   }
};

class DGProlongation : public ProlongationOperator< dg::DGFunction< real_t > >
{
 public:
   explicit DGProlongation( std::shared_ptr< ProlongationForm > form )
   : form_( form )
   {}

   void prolongate( const dg::DGFunction< real_t >& function,
                    const walberla::uint_t&         coarseLevel,
                    const DoFType&                  flag ) const override;

 protected:
   std::shared_ptr< ProlongationForm > form_;
};

class DG1Prolongation : public DGProlongation
{
 public:
   DG1Prolongation()
   : DGProlongation( std::make_shared< ProlongationFormDG1 >() )
   {}
};

class DG0Prolongation : public DGProlongation
{
 public:
   DG0Prolongation()
   : DGProlongation( std::make_shared< ProlongationFormDG0 >() )
   {}
};

} // namespace hyteg
