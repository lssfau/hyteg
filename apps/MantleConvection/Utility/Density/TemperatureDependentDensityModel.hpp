/*
 * Copyright (c) 2024-2025 Andreas Burkhart.
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

#include "core/DataTypes.h"
#include "core/config/Config.h"

#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/types/PointND.hpp"

#include "../Parameters/NondimensionalisationParameters.hpp"

namespace MantleConvection {

template < class ValueType >
class TemperatureDependentDensityModel
{
 public:
   virtual ~TemperatureDependentDensityModel() = default;

   TemperatureDependentDensityModel( std::string prefix = "" )
   : prefix_( prefix )
   {}

   virtual ValueType evaluate( const hyteg::Point3D& x, ValueType temp )
   {
      WALBERLA_UNUSED( x );
      WALBERLA_UNUSED( temp );
      WALBERLA_ABORT( "evaluate not implemented for TemperatureDependentDensityModel abstract base class!" );
      return ValueType();
   }

   virtual std::function< real_t( const hyteg::Point3D&, const std::vector< real_t >& ) > getDensityFct()
   {
      return [this]( const hyteg::Point3D& x, const std::vector< real_t >& fields ) {
         if ( fields.size() > 0 )
         {
            return this->evaluate( x, std::reduce( fields.begin(), fields.end() ) );
         }
         return this->evaluate( x, real_c( 0.0 ) );
      };
   }

   virtual std::ostream& print( std::ostream& os, uint_t offset = 0 ) const
   {
      WALBERLA_UNUSED( offset );
      WALBERLA_ABORT( "print not implemented for DensityModel abstract base class!" );
      return os;
   }

   template < class FunctionType >
   void interpolate( FunctionType&                                                      f,
                     const std::vector< std::reference_wrapper< const FunctionType > >& srcFunctions,
                     uint_t                                                             level,
                     hyteg::DoFType                                                     flag )
   {
      f.interpolate( getDensityFct(), srcFunctions, level, flag );
   }

   template < class FunctionType >
   void interpolate( FunctionType& f, uint_t level, hyteg::DoFType flag )
   {
      interpolate( f, {}, level, flag );
   }

 protected:
   std::string prefix_;
};

} // namespace MantleConvection