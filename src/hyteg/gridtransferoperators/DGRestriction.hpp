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

   /// \brief Evaluates the coupling of DoF between a coarse triangle and (and embedded) fine triangle.
   ///
   /// \param coarseTriangle Vertices of the coarse triangle, ordered w.r.t. the DoF.
   /// \param fineTriangle   Vertices of the fine triangle, ordered w.r.t. the DoF.
   /// \param localMat       Coupling between coarse triangle DoF (rows) and fine triangle DoF (cols).
   virtual void integrate2D( const std::vector< Point >&                              coarseTriangle,
                             const std::vector< Point >&                              fineTriangle,
                             Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& localMat ) const
   {
      WALBERLA_UNUSED( coarseTriangle );
      WALBERLA_UNUSED( fineTriangle );
      WALBERLA_UNUSED( localMat );
      WALBERLA_LOG_INFO_ON_ROOT( "not implemented" )
   }
};

/// \brief Restriction form for piecewise linear DG functions.
class RestrictionFormDG1 : public RestrictionForm
{
 public:
   void integrate2D( const std::vector< Point >&                              dst,
                     const std::vector< Point >&                              src,
                     Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& localMat ) const override;
};

/// \brief Restriction form for piecewise constant DG functions.
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

   /// \brief Restrictics the given function on the fine level to the coarser level.
   ///
   /// \param function   Function which we restrict to the coarse level.
   /// \param srcLevel   The fine level.
   /// \param flag       Is ignored, but needed to be consistent with the other restriction operators.
   void restrict( const dg::DGFunction< real_t >& function,
                  const walberla::uint_t&         srcLevel,
                  const DoFType&                  flag ) const override;

 protected:
   std::shared_ptr< RestrictionForm > form_;
};

/// \brief Restriction for piecewise linear DG functions.
class DG1Restriction : public DGRestriction
{
 public:
   DG1Restriction()
   : DGRestriction( std::make_shared< RestrictionFormDG1 >() )
   {}
};

/// \brief Restriction for piecewise constant DG functions.
class DG0Restriction : public DGRestriction
{
 public:
   DG0Restriction()
   : DGRestriction( std::make_shared< RestrictionFormDG0 >() )
   {}
};

} // namespace hyteg
