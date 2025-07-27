/*
 * Copyright (c) 2017-2025 Dominik Thoennes, Nils Kohl, Andreas Burkhart.
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
#include "hyteg/functions/PressureMeanProjection.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticVectorProlongation.hpp"
#include "hyteg/gridtransferoperators/ProlongationOperator.hpp"
#include "hyteg/operators/NoOperator.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
namespace hyteg {

class P2P1StokesToP2P1StokesProlongation : public ProlongationOperator< P2P1TaylorHoodFunction< real_t > >
{
 public:
   typedef P2toP2QuadraticProlongation VelocityProlongation_T;
   typedef P1toP1LinearProlongation<>  PressureProlongation_T;

   P2P1StokesToP2P1StokesProlongation()
   : projectMeanAfterProlongation_( false )
   {}
   P2P1StokesToP2P1StokesProlongation( bool projectMeanAfterProlongation )
   : projectMeanAfterProlongation_( projectMeanAfterProlongation )
   {}

   void prolongate( const P2P1TaylorHoodFunction< real_t >& function,
                    const uint_t&                           sourceLevel,
                    const DoFType&                          flag ) const override
   {
      for ( uint_t k = 0; k < function.uvw().getDimension(); k++ )
      {
         quadraticProlongationOperator_.prolongate( function.uvw()[k], sourceLevel, flag );
      }
      linearProlongationOperator_.prolongate( function.p(), sourceLevel, flag );

      if ( projectMeanAfterProlongation_ )
      {
         vertexdof::projectMean( function.p(), sourceLevel + 1 );
      }
   }

   void prolongateAndAdd( const P2P1TaylorHoodFunction< real_t >& function,
                          const uint_t&                           sourceLevel,
                          const DoFType&                          flag ) const override
   {
      for ( uint_t k = 0; k < function.uvw().getDimension(); k++ )
      {
         quadraticProlongationOperator_.prolongateAndAdd( function.uvw()[k], sourceLevel, flag );
      }
      linearProlongationOperator_.prolongateAndAdd( function.p(), sourceLevel, flag );

      if ( projectMeanAfterProlongation_ )
      {
         vertexdof::projectMean( function.p(), sourceLevel + 1 );
      }
   }

 private:
   P2toP2QuadraticProlongation quadraticProlongationOperator_;
   P1toP1LinearProlongation<>  linearProlongationOperator_;

   bool projectMeanAfterProlongation_;
};

template < typename ProjectionOperatorType             = hyteg::NoOperator,
           bool preProjectVelocity                     = false,
           bool postProjectVelocity                    = true,
           bool allowPreProjectionToChangeVelocitySrc  = true,
           bool allowPostProjectionToChangeVelocityDst = true,
           bool preProjectPressure                     = false,
           bool postProjectPressure                    = true,
           bool allowPreProjectionToChangePressureSrc  = true,
           bool allowPostProjectionToChangePressureDst = true >
class P2P1StokesToP2P1StokesProlongationWithProjection : public ProlongationOperator< P2P1TaylorHoodFunction< real_t > >
{
 public:
   P2P1StokesToP2P1StokesProlongationWithProjection( const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
                                                     uint_t                                            minLevel,
                                                     uint_t                                            maxLevel,
                                                     std::shared_ptr< ProjectionOperatorType >         projection = nullptr,
                                                     DoFType projectionFlag = FreeslipBoundary,
                                                     bool    lowMemoryMode  = false )
   : quadraticProlongationOperator_( storage, minLevel, maxLevel, projection, projectionFlag, lowMemoryMode )
   , linearProlongationOperator_( storage, minLevel, maxLevel, lowMemoryMode )
   {}

   void prolongate( const P2P1TaylorHoodFunction< real_t >& function,
                    const uint_t&                           sourceLevel,
                    const DoFType&                          flag ) const override
   {
      quadraticProlongationOperator_.prolongate( function.uvw(), sourceLevel, flag );
      linearProlongationOperator_.prolongate( function.p(), sourceLevel, flag );
   }

   // prolongateAndAdd has a different implementation than prolongate and seemingly needs internal projections
   // as a quick fix we are using prolongate on a temporary function and add it manually
   void prolongateAndAdd( const P2P1TaylorHoodFunction< real_t >& function,
                          const uint_t&                           sourceLevel,
                          const DoFType&                          flag ) const override
   {
      quadraticProlongationOperator_.prolongateAndAdd( function.uvw(), sourceLevel, flag );
      linearProlongationOperator_.prolongateAndAdd( function.p(), sourceLevel, flag );
   }

 private:
   P2toP2QuadraticVectorProlongationWithProjection< ProjectionOperatorType,
                                                    preProjectVelocity,
                                                    postProjectVelocity,
                                                    allowPreProjectionToChangeVelocitySrc,
                                                    allowPostProjectionToChangeVelocityDst >
       quadraticProlongationOperator_;
   P1toP1LinearProlongationWithProjection< real_t,
                                           preProjectPressure,
                                           postProjectPressure,
                                           allowPreProjectionToChangePressureSrc,
                                           allowPostProjectionToChangePressureDst >
       linearProlongationOperator_;
};

/***************************************************************************
NOTE: This prolongates the FE function and calls the project function on it 
      so that the normal components are set to zero on the FreeslipBoundary
***************************************************************************/
class P2P1StokesToP2P1StokesProlongationWithFreeSlipProjection : public P2P1StokesToP2P1StokesProlongation
{
 public:
   P2P1StokesToP2P1StokesProlongationWithFreeSlipProjection( std::shared_ptr< P2P1TaylorHoodFunction< real_t > > temp,
                                                             std::shared_ptr< P2ProjectNormalOperator >          projection )
   : P2P1StokesToP2P1StokesProlongation()
   , temp_( temp )
   , projection_( projection )
   {}

   void prolongate( const P2P1TaylorHoodFunction< real_t >& function,
                    const uint_t&                           sourceLevel,
                    const DoFType&                          flag ) const override
   {
      P2P1StokesToP2P1StokesProlongation::prolongate( function, sourceLevel, flag );
      projection_->project( function, sourceLevel + 1, FreeslipBoundary );
      vertexdof::projectMean( function.p(), sourceLevel + 1 );
   }

   // prolongateAndAdd has a different implementation than prolongate and seemingly needs internal projections
   // as a quick fix we are using prolongate on a temporary function and add it manually
   void prolongateAndAdd( const P2P1TaylorHoodFunction< real_t >& function,
                          const uint_t&                           sourceLevel,
                          const DoFType&                          flag ) const override
   {
      temp_->assign( { 1.0 }, { function }, sourceLevel, All );
      prolongate( *temp_, sourceLevel, flag );
      function.assign( { 1.0, 1.0 }, { function, *temp_ }, sourceLevel + 1, flag );
   }

 private:
   std::shared_ptr< P2P1TaylorHoodFunction< real_t > > temp_;
   std::shared_ptr< P2ProjectNormalOperator >          projection_;
};

} // namespace hyteg
