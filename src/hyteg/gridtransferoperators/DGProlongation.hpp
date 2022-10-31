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

   virtual ~ProlongationForm() = default;

   /// \brief Evaluates the coupling of DoF between a coarse triangle and (and embedded) fine triangle.
   ///
   /// \param fineTriangle   Vertices of the fine triangle, ordered w.r.t. the DoF.
   /// \param coarseTriangle Vertices of the coarse triangle, ordered w.r.t. the DoF.
   /// \param localMat       Coupling between fine triangle DoF (rows) and coarse triangle DoF (cols).
   virtual void integrate2D( const std::vector< Point >&                              fineTriangle,
                             const std::vector< Point >&                              coarseTriangle,
                             Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& localMat ) const
   {
      WALBERLA_UNUSED( fineTriangle );
      WALBERLA_UNUSED( coarseTriangle );
      WALBERLA_UNUSED( localMat );
      WALBERLA_LOG_INFO_ON_ROOT( "not implemented" )
   }
};

/// \brief Prolongation form for piecewise linear DG functions.
class ProlongationFormDG1 : public ProlongationForm
{
 public:
   void integrate2D( const std::vector< Point >&                              dst,
                     const std::vector< Point >&                              src,
                     Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& localMat ) const override;
};

/// \brief Prolongation form for piecewise constant DG functions.
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

   /// \brief Prolongates the given function on the coarse level to a finer level.
   ///
   /// \param function   Function which we prolongate to a finer level.
   /// \param srcLevel   The coarse level.
   /// \param flag       Is ignored, but needed to be consistent with the other prolongation operators.
   void prolongate( const dg::DGFunction< real_t >& function,
                    const walberla::uint_t&         srcLevel,
                    const DoFType&                  flag ) const override;

 protected:
   std::shared_ptr< ProlongationForm > form_;
};

/// \brief Prolongation for piecewise linear DG functions.
class DG1Prolongation : public DGProlongation
{
 public:
   DG1Prolongation()
   : DGProlongation( std::make_shared< ProlongationFormDG1 >() )
   {}
};

/// \brief Prolongation for piecewise constant DG functions.
class DG0Prolongation : public DGProlongation
{
 public:
   DG0Prolongation()
   : DGProlongation( std::make_shared< ProlongationFormDG0 >() )
   {}
};

} // namespace hyteg
