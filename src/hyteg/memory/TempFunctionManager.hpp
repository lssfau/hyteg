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

#include <any>

#include "core/DataTypes.h"
#include "core/config/Config.h"
#include "core/math/Random.h"
#include "core/singleton/Singleton.h"

#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/types/types.hpp"

// clang-format off
#include "hyteg/dg1functionspace/DG1Function.hpp"
#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/egfunctionspace/EGFunction.hpp"
#include "hyteg/n1e1functionspace/N1E1VectorFunction.hpp"
#include "hyteg/p0functionspace/P0Function.hpp"
#include "hyteg/p1functionspace/VertexDoFFunction.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/volumedofspace/VolumeDoFFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P1DGEP0StokesFunction.hpp"
#include "hyteg/composites/P1P0StokesFunction.hpp"
#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/composites/P2P2StokesFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodBlockFunction.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
// clang-format on

using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

class TempFunctionManager : public walberla::singleton::Singleton< TempFunctionManager >
{
 public:
   WALBERLA_BEFRIEND_SINGLETON;

   ~TempFunctionManager() = default;

   template < class FunctionType >
   std::shared_ptr< FunctionType >
       getFunction( const std::shared_ptr< hyteg::PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   {
      // DG1 Functions
      if constexpr ( std::is_same_v< FunctionType, DG1Function< typename FunctionType::valueType > > )
      {
         return getFunctionInternal< FunctionType >( DG1Functions_, "DG1", storage, minLevel, maxLevel );
      }
      // EdgeDof Functions
      else if constexpr ( std::is_same_v< FunctionType, EdgeDoFFunction< typename FunctionType::valueType > > )
      {
         return getFunctionInternal< FunctionType >( EdgeDoFFunctions_, "EdgeDoF", storage, minLevel, maxLevel );
      }
      // EG Functions
      else if constexpr ( std::is_same_v< FunctionType, EGFunction< typename FunctionType::valueType > > )
      {
         return getFunctionInternal< FunctionType >( EGFunctions_, "EG", storage, minLevel, maxLevel );
      }
      // N1E1Vector Functions
      else if constexpr ( std::is_same_v< FunctionType, n1e1::N1E1VectorFunction< typename FunctionType::valueType > > )
      {
         return getFunctionInternal< FunctionType >( N1E1VectorFunctions_, "N1E1Vector", storage, minLevel, maxLevel );
      }
      // P0 Functions
      else if constexpr ( std::is_same_v< FunctionType, P0Function< typename FunctionType::valueType > > )
      {
         return getFunctionInternal< FunctionType >( P0Functions_, "P0", storage, minLevel, maxLevel );
      }
      // VertexDoF Functions
      else if constexpr ( std::is_same_v< FunctionType, vertexdof::VertexDoFFunction< typename FunctionType::valueType > > )
      {
         return getFunctionInternal< FunctionType >( VertexDoFFunctions_, "VertexDoF", storage, minLevel, maxLevel );
      }
      // P2 Functions
      else if constexpr ( std::is_same_v< FunctionType, P2Function< typename FunctionType::valueType > > )
      {
         return getFunctionInternal< FunctionType >( P2Functions_, "P2", storage, minLevel, maxLevel );
      }
      // VolumeDoF Functions
      else if constexpr ( std::is_same_v< FunctionType, VolumeDoFFunction< typename FunctionType::valueType > > )
      {
         return getFunctionInternal< FunctionType >( VolumeDoFFunctions_, "VolumeDoF", storage, minLevel, maxLevel );
      }
      // P2P1TaylorHood Functions
      else if constexpr ( std::is_same_v< FunctionType, P2P1TaylorHoodFunction< typename FunctionType::valueType > > )
      {
         return getFunctionInternal< FunctionType >( P2P1TaylorHoodFunctions_, "P2P1TaylorHood", storage, minLevel, maxLevel );
      }
      // P1DGEP0Stokes Functions
      else if constexpr ( std::is_same_v< FunctionType, EGP0StokesFunction< typename FunctionType::valueType > > )
      {
         return getFunctionInternal< FunctionType >( P1DGEP0StokesFunctions_, "P1DGEP0Stokes", storage, minLevel, maxLevel );
      }
      // P1P0Stokes Functions
      else if constexpr ( std::is_same_v< FunctionType, P1P0StokesFunction< typename FunctionType::valueType > > )
      {
         return getFunctionInternal< FunctionType >( P1P0StokesFunctions_, "P1P0Stokes", storage, minLevel, maxLevel );
      }
      // P1Stokes Functions
      else if constexpr ( std::is_same_v< FunctionType, P1StokesFunction< typename FunctionType::valueType > > )
      {
         return getFunctionInternal< FunctionType >( P1StokesFunctions_, "P1Stokes", storage, minLevel, maxLevel );
      }
      // P2P2Stokes Functions
      else if constexpr ( std::is_same_v< FunctionType, P2P2StokesFunction< typename FunctionType::valueType > > )
      {
         return getFunctionInternal< FunctionType >( P2P2StokesFunctions_, "P2P2Stokes", storage, minLevel, maxLevel );
      }
      // P2P1TaylorHoodBlock Functions
      else if constexpr ( std::is_same_v< FunctionType, P2P1TaylorHoodBlockFunction< typename FunctionType::valueType > > )
      {
         return getFunctionInternal< FunctionType >(
             P2P1TaylorHoodBlockFunctions_, "P2P1TaylorHoodBlock", storage, minLevel, maxLevel );
      }
      // P2Vector Functions
      else if constexpr ( std::is_same_v< FunctionType, P2VectorFunction< typename FunctionType::valueType > > )
      {
         return getFunctionInternal< FunctionType >( P2VectorFunctions_, "P2Vector", storage, minLevel, maxLevel );
      }
      // P1Vector Functions
      else if constexpr ( std::is_same_v< FunctionType, P1VectorFunction< typename FunctionType::valueType > > )
      {
         return getFunctionInternal< FunctionType >( P1VectorFunctions_, "P1Vector", storage, minLevel, maxLevel );
      }
      // Misc Functions
      else
      {
         return getFunctionInternal< FunctionType >( miscFunctions_, "Misc", storage, minLevel, maxLevel );
      }
   }

   /// if true the manager will always return return
   /// shared ptrs to newly generated functions
   /// that automatically free themselves after use
   void setAlwaysDestroy( bool alwaysDestroy ) { alwaysDestroy_ = alwaysDestroy; }

 private:
   TempFunctionManager()
   : alwaysDestroy_( false )
   {}

   template < class FunctionType >
   inline std::shared_ptr< FunctionType > getFunctionInternal( std::vector< std::any >&                          vec,
                                                               const std::string&                                name,
                                                               const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
                                                               uint_t                                            minLevel,
                                                               uint_t                                            maxLevel )
   {
      if (alwaysDestroy_)
      {
         return std::make_shared< FunctionType >( "temporary function", storage, minLevel, maxLevel );
      }

      // Try to find existing function that can be used
      for ( auto& f : vec )
      {
         if ( f.type() == typeid( std::shared_ptr< FunctionType > ) )
         {
            std::shared_ptr< FunctionType >& fptr = std::any_cast< std::shared_ptr< FunctionType >& >( f );

            // check if we already have a function that matches our criteria and that is not in use
            if ( fptr.use_count() <= 1 && storage.get() == fptr->getStorage().get() && fptr->getMinLevel() <= minLevel &&
                 fptr->getMaxLevel() >= maxLevel )
            {
               return fptr;
            }
         }
      }

      // Create a new function if necessary
      std::stringstream fname;
      fname << name << "_function" << vec.size();

      auto ptr = std::make_shared< FunctionType >( fname.str(), storage, minLevel, maxLevel );
      vec.push_back( ptr );

      return ptr;
   }

   // we are using std::any here because want containers that can hold a function of any value type
   std::vector< std::any > DG1Functions_;
   std::vector< std::any > EdgeDoFFunctions_;
   std::vector< std::any > EGFunctions_;
   std::vector< std::any > N1E1VectorFunctions_;
   std::vector< std::any > P0Functions_;
   std::vector< std::any > VertexDoFFunctions_;
   std::vector< std::any > P2Functions_;
   std::vector< std::any > VolumeDoFFunctions_;
   std::vector< std::any > P2P1TaylorHoodFunctions_;
   std::vector< std::any > P1DGEP0StokesFunctions_;
   std::vector< std::any > P1P0StokesFunctions_;
   std::vector< std::any > P1StokesFunctions_;
   std::vector< std::any > P2P2StokesFunctions_;
   std::vector< std::any > P2P1TaylorHoodBlockFunctions_;
   std::vector< std::any > P2VectorFunctions_;
   std::vector< std::any > P1VectorFunctions_;
   std::vector< std::any > miscFunctions_;

   bool alwaysDestroy_;
};

/// \brief Returns a temporary function of type FunctionType from a reusable pool of functions.
///
/// The returned function might have arbitrary boundary conditions and contain arbitrary information.
/// You are responsible yourself to make the function usable for your purpose (e.g. set boundary conditions and set the function to zero).
/// Once the shared_ptr obtained through this function is no longer referenced outside of the temporary function manager the function returns to the pool and can be reused.
///
/// \param storage  A PrimitiveStorage instance.
/// \param minLevel Minimum level of the function.
/// \param maxLevel Maximum level of the function.
/// \param destroyAfterUse If you want the function to be freed from memory no matter what after use, set this to true
///
/// \return Shared pointer to the function.
template < class FunctionType >
std::shared_ptr< FunctionType > getTemporaryFunction( const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
                                                      uint_t                                            minLevel,
                                                      uint_t                                            maxLevel,
                                                      bool                                              destroyAfterUse = false )
{
   if ( destroyAfterUse )
   {
      return std::make_shared< FunctionType >( "temporary function", storage, minLevel, maxLevel );
   }
   else
   {
      return TempFunctionManager::instance()->getFunction< FunctionType >( storage, minLevel, maxLevel );
   }
}

} // namespace hyteg