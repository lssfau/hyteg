/*
 * Copyright (c) 2024 Andreas Burkhart.
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

#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

#include "mixed_operator/VectorToVectorOperator.hpp"

namespace hyteg {

using walberla::real_t;

template < class MassOperatorType >
class P2VectorMassOperator : public VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >,
                             public OperatorWithInverseDiagonal< P2VectorFunction< real_t > >,
                             public WeightedJacobiSmoothable< P2VectorFunction< real_t > >
{
 public:
   P2VectorMassOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                         size_t                                     minLevel,
                         size_t                                     maxLevel,
                         const std::shared_ptr< MassOperatorType >& MassOperator )
   : VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >( storage, minLevel, maxLevel )
   , MassOperator_( MassOperator )
   {
      if ( this->dim_ == 3 )
      {
         this->subOper_[0][0] = MassOperator_;
         this->subOper_[0][1] = nullptr;
         this->subOper_[0][2] = nullptr;

         this->subOper_[1][0] = nullptr;
         this->subOper_[1][1] = MassOperator_;
         this->subOper_[1][2] = nullptr;

         this->subOper_[2][0] = nullptr;
         this->subOper_[2][1] = nullptr;
         this->subOper_[2][2] = MassOperator_;
      }
      else
      {
         this->subOper_[0][0] = MassOperator_;
         this->subOper_[0][1] = nullptr;

         this->subOper_[1][0] = nullptr;
         this->subOper_[1][1] = MassOperator_;
      }
   }

   void apply( const SrcVecFuncType& src,
               const DstVecFuncType& dst,
               size_t                level,
               DoFType               flag,
               UpdateType            updateType = Replace ) const override
   {
      VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >::apply( src, dst, level, flag, updateType );
   };

   void smooth_jac( const P2VectorFunction< real_t >& dst,
                    const P2VectorFunction< real_t >& rhs,
                    const P2VectorFunction< real_t >& src,
                    real_t                            relax,
                    size_t                            level,
                    DoFType                           flag ) const override
   {
      for ( uint_t k = 0; k < this->dim_; ++k )
      {
         if ( const auto subOp =
                  std::dynamic_pointer_cast< WeightedJacobiSmoothable< P2Function< real_t > > >( this->subOper_[k][k] ) )
         {
            subOp->smooth_jac( dst[k], rhs[k], src[k], relax, level, flag );
         }
         else
         {
            throw std::runtime_error(
                "Jacobi smoothing of a VectorDivKGradOperator requires its diagonal blocks to have the WeightedJacobiSmoothable interface." );
         }
      }
   }

   std::shared_ptr< P2VectorFunction< real_t > > getInverseDiagonalValues() const override final
   {
      return this->extractInverseDiagonal();
   }

   void computeInverseDiagonalOperatorValues() override final
   {
      this->VectorToVectorOperator< real_t, P2VectorFunction, P2VectorFunction >::computeInverseDiagonalOperatorValues();
   }

 private:
   std::shared_ptr< MassOperatorType > MassOperator_;
};

} // namespace hyteg
