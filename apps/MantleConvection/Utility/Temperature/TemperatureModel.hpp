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
#include "core/Environment.h"
#include "core/config/Config.h"

#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/types/PointND.hpp"

#include "../Parameters/NondimensionalisationParameters.hpp"

namespace MantleConvection {

template < class ValueType >
class TemperatureModel
{
 public:
   virtual ~TemperatureModel() = default;

   TemperatureModel( std::string prefix = "" )
   : prefix_( prefix )
   {}

   virtual ValueType evaluate( const hyteg::Point3D& x )
   {
      WALBERLA_UNUSED( x );
      WALBERLA_ABORT( "evaluate not implemented for TemperatureModel abstract base class!" );
      return ValueType();
   }

   virtual ValueType evaluateNoCheck( const hyteg::Point3D& x )
   {
      WALBERLA_UNUSED( x );
      WALBERLA_ABORT( "evaluateNoCheck not implemented for TemperatureModel abstract base class!" );
      return ValueType();
   }

   virtual std::function< real_t( const hyteg::Point3D& ) > getTemperatureFct()
   {
      return [this]( const hyteg::Point3D& x ) { return this->evaluate( x ); };
   }

   virtual std::function< real_t( const hyteg::Point3D& ) > getTemperatureFctNoCheck()
   {
      return [this]( const hyteg::Point3D& x ) { return this->evaluateNoCheck( x ); };
   }   

   virtual bool checkSurface( const hyteg::Point3D& x )
   {
      WALBERLA_UNUSED( x );
      WALBERLA_ABORT( "checkSurface not implemented for TemperatureModel abstract base class!" );
      return false;
   }

   virtual bool checkCMB( const hyteg::Point3D& x )
   {
      WALBERLA_UNUSED( x );
      WALBERLA_ABORT( "checkCMB not implemented for TemperatureModel abstract base class!" );
      return false;
   }

   virtual std::ostream& print( std::ostream& os, uint_t offset = 0 ) const
   {
      WALBERLA_UNUSED( offset );
      WALBERLA_ABORT( "print not implemented for TemperatureModel abstract base class!" );
      return os;
   }

   template < class FunctionType >
   void interpolate( FunctionType& f, uint_t level, hyteg::DoFType flag )
   {
      f.interpolate( getTemperatureFct(), level, flag );
   }

 protected:
   std::string prefix_;
};

} // namespace MantleConvection