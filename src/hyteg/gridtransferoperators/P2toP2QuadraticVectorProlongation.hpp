/*
 * Copyright (c) 2023 Andreas Burkhart.
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
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/gridtransferoperators/ProlongationOperator.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/p2functionspace/P2RotationOperator.hpp"

namespace hyteg {

class P2toP2QuadraticVectorProlongation : public ProlongationOperator< P2VectorFunction< real_t > >
{
 public:
   void prolongate( const P2VectorFunction< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      for ( uint_t k = 0; k < function.getDimension(); k++ )
      {
         quadraticProlongationOperator_.prolongate( function[k], sourceLevel, flag );
      }
   }

   void prolongateAndAdd( const P2VectorFunction< real_t >& function,
                          const uint_t&                     sourceLevel,
                          const DoFType&                    flag ) const override
   {
      for ( uint_t k = 0; k < function.getDimension(); k++ )
      {
         quadraticProlongationOperator_.prolongateAndAdd( function[k], sourceLevel, flag );
      }
   }

 private:
   P2toP2QuadraticProlongation quadraticProlongationOperator_;
};

/***************************************************************************
NOTE: This prolongates the FE function and calls the project function on it 
      so that the normal components are set to zero on the FreeslipBoundary
***************************************************************************/
class P2toP2QuadraticVectorProlongationWithFreeSlipProjection : public P2toP2QuadraticVectorProlongation
{
 public:
   P2toP2QuadraticVectorProlongationWithFreeSlipProjection( std::shared_ptr< P2VectorFunction< real_t > > temp,
                                                    std::shared_ptr< P2ProjectNormalOperator >    projection )
   : P2toP2QuadraticVectorProlongation()
   , temp_( temp )
   , projection_( projection )
   {}

   void prolongate( const P2VectorFunction< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      P2toP2QuadraticVectorProlongation::prolongate( function, sourceLevel, flag );
      projection_->project( function, sourceLevel + 1, FreeslipBoundary );
   }

   // prolongateAndAdd has a different implementation than prolongate and seemingly needs internal projections
   // as a quick fix we are using prolongate on a temporary function and add it manually
   void prolongateAndAdd( const P2VectorFunction< real_t >& function,
                          const uint_t&                     sourceLevel,
                          const DoFType&                    flag ) const override
   {
      temp_->assign( { 1.0 }, { function }, sourceLevel, All );
      prolongate( *temp_, sourceLevel, flag );
      function.assign( { 1.0, 1.0 }, { function, *temp_ }, sourceLevel + 1, flag );
   }

 private:
   std::shared_ptr< P2VectorFunction< real_t > > temp_;
   std::shared_ptr< P2ProjectNormalOperator >    projection_;
};

class P2toP2QuadraticVectorProlongationWithRotation : public P2toP2QuadraticVectorProlongation
{
 public:
   P2toP2QuadraticVectorProlongationWithRotation( std::shared_ptr< P2VectorFunction< real_t > > temp,
                                                  std::shared_ptr< P2VectorFunction< real_t > > temp1,
                                                  std::shared_ptr< P2RotationOperator >         rotation )
   : P2toP2QuadraticVectorProlongation()
   , temp_( temp )
   , temp1_( temp1 )
   , rotation_( rotation )
   {}

   void prolongate( const P2VectorFunction< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      temp_->assign( { 1.0 }, { function }, sourceLevel, All );
      rotation_->rotate( *temp_, sourceLevel, FreeslipBoundary, true );
      P2toP2QuadraticVectorProlongation::prolongate( *temp_, sourceLevel, flag );
      // removeRotationalModes( *temp_, sourceLevel + 1 );
      rotation_->rotate( *temp_, sourceLevel + 1, FreeslipBoundary, false );
      function.assign( { 1.0 }, { *temp_ }, sourceLevel + 1, All );
   }

   // prolongateAndAdd has a different implementation than prolongate and seemingly needs internal projections
   // as a quick fix we are using prolongate on a temporary function and add it manually
   void prolongateAndAdd( const P2VectorFunction< real_t >& function,
                          const uint_t&                     sourceLevel,
                          const DoFType&                    flag ) const override
   {
      temp1_->assign( { 1.0 }, { function }, sourceLevel, All );
      prolongate( *temp1_, sourceLevel, flag );
      function.assign( { 1.0, 1.0 }, { function, *temp1_ }, sourceLevel + 1, flag );
   }

 private:
   std::shared_ptr< P2VectorFunction< real_t > > temp_;
   std::shared_ptr< P2VectorFunction< real_t > > temp1_;
   std::shared_ptr< P2RotationOperator >         rotation_;
};
} // namespace hyteg
