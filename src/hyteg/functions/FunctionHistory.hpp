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
#include <algorithm>
#include <fstream>
#include <sstream>

#include "core/DataTypes.h"
#include "core/config/Config.h"

#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/indexing/ConsistentEnumeration.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/sparseassembly/FileWritingVector.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"
#include "hyteg/types/types.hpp"

#ifdef HYTEG_BUILD_WITH_ADIOS2
#include "hyteg/checkpointrestore/ADIOS2/AdiosCheckpointExporter.hpp"
#include "hyteg/checkpointrestore/ADIOS2/AdiosCheckpointImporter.hpp"
#endif

using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

// Object containing multiple functions (e.g. states at certain time points) and associated data (e.g. the time step sizes between time points).
// The ability to get a shared pointer to a specific state offset and its data (even to a specific time step) and the fact that returned function references
// remain bound to a specific state offset is intentional.
// Holds a memoryCapacity many functions of a type FunctionType.
// When adding a new state, the oldest state information saved in a function history is lost.
// Usage:
// // Create a new function history with all state FEM functions already initialised as empty functions with boundary condition bc
// hyteg::FunctionHistory< FunctionType > > Hist("f", storage, minLevel, maxLevel, bc, memoryCapacity );
// // Introduce a new state
// real_t stepSize = ...;
// Hist.newState( stepSize );
// // Now you can acquire and manipulate the newest state via Hist.getState( 0 ), Hist.getStateStepSize( 0 ) and Hist.getStateAdditionalData( 0 )
// auto& state0 = Hist.getState( 0 );
// state0.assign(1, maxLevel, All);
// // After introducing a new state the old one gets shifted one offset back
// Hist.newState( stepSize );
// Hist.getState( 0 ).assign(2, maxLevel, All);
// // state0 still refers to the function at offset 0, i.e. the function now containing values of 2
// // If you want to access the past time step you need to use Hist.getState( 1 )
// Hist.getState( 1 ).assign(1, maxLevel, All);
template < class FunctionType, class AdditionalDataType = real_t >
class FunctionHistory
{
 public:
   // automatically create a vector of functions with identical boundary condition
   FunctionHistory( const std::string&                                name,
                    const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
                    uint_t                                            minLevel,
                    uint_t                                            maxLevel,
                    hyteg::BoundaryCondition                          bc,
                    uint_t                                            memoryCapacity )
   : memoryCapacity_( memoryCapacity )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , functions_( memoryCapacity )
   , stepSizes_( memoryCapacity )
   , additionalData_( memoryCapacity )
   , numberOfInitialisedStates_( 0 )
   {
      for ( uint_t i = 0; i < memoryCapacity; i++ )
      {
         std::stringstream fname;
         fname << name << " index " << i;
         std::shared_ptr< FunctionType > f = std::make_shared< FunctionType >( fname.str(), storage, minLevel, maxLevel, bc );
         functions_.at( i )                = f;
         std::shared_ptr< AdditionalDataType > additional = std::make_shared< AdditionalDataType >();
         additionalData_.at( i )                          = additional;
         std::shared_ptr< real_t > step                   = std::make_shared< real_t >();
         stepSizes_.at( i )                               = step;
      }
   }

   // add functions yourself
   // should be used if vector functions with different boundary conditions for each component are required
   FunctionHistory( uint_t memoryCapacity, uint_t minLevel, uint_t maxLevel )
   : memoryCapacity_( memoryCapacity )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , functions_( 0 )
   , stepSizes_( memoryCapacity )
   , additionalData_( memoryCapacity )
   , numberOfInitialisedStates_( 0 )
   {
      for ( uint_t i = 0; i < memoryCapacity; i++ )
      {
         std::shared_ptr< AdditionalDataType > additional = std::make_shared< AdditionalDataType >();
         additionalData_.at( i )                          = additional;
         std::shared_ptr< real_t > step                   = std::make_shared< real_t >();
         stepSizes_.at( i )                               = step;
      }
   }

   void addFunction( const std::shared_ptr< FunctionType >& f )
   {
      if ( ( f->getMinLevel() != minLevel_ ) || ( f->getMaxLevel() != maxLevel_ ) )
      {
         WALBERLA_ABORT( "minLevel or maxLevel mismatch in addFunction!" )
      }

      if ( functions_.size() < memoryCapacity_ )
      {
         functions_.push_back( f );
      }
      else
      {
         WALBERLA_LOG_WARNING(
             "Ignored addFunction since it tried to add more functions to a FunctionMemory than its inherit capacity of "
             << memoryCapacity_ << "!" );
      }
   }

   void saveFunctionsViaFileWritingVector( std::string dir = "./" )
   {
      if ( functions_.size() < memoryCapacity_ )
      {
         WALBERLA_ABORT( "Function History is not properly initialised." );
      }

      // temporarily create consistent enumerator
      auto enumerator = std::make_shared< typename FunctionType::template FunctionType< hyteg::idx_t > >(
          "enumerator", functions_.at( 0 )->getStorage(), minLevel_, maxLevel_ );

      for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
      {
         enumerateConsistently( *enumerator, level );
      }

      for ( uint_t i = 0; i < memoryCapacity_; i++ )
      {
         for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
         {
            auto vec = std::make_shared< hyteg::FileWritingVector< FunctionType > >(
                functions_.at( i )->getStorage(), level, *enumerator, functions_.at( i ) );

            std::stringstream fName;
            std::string       functionName = functions_.at( i )->getFunctionName();
            std::replace( functionName.begin(), functionName.end(), ' ', '_' );
            fName << dir << functionName << "_Level" << level << ".bin";

            vec->writeToFile( fName.str() );
         }
      }
   }

   void loadFunctionsViaFileWritingVector( std::string dir = "./" )
   {
      if ( functions_.size() < memoryCapacity_ )
      {
         WALBERLA_ABORT( "Function History is not properly initialised." );
      }

      // temporarily create consistent enumerator
      auto enumerator = std::make_shared< typename FunctionType::template FunctionType< hyteg::idx_t > >(
          "enumerator", functions_.at( 0 )->getStorage(), minLevel_, maxLevel_ );

      for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
      {
         enumerateConsistently( *enumerator, level );
      }

      for ( uint_t i = 0; i < memoryCapacity_; i++ )
      {
         for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
         {
            auto vec = std::make_shared< hyteg::FileWritingVector< FunctionType > >(
                functions_.at( i )->getStorage(), level, *enumerator, functions_.at( i ) );

            std::stringstream fName;
            std::string       functionName = functions_.at( i )->getFunctionName();
            std::replace( functionName.begin(), functionName.end(), ' ', '_' );
            fName << dir << functionName << "_Level" << level << ".bin";

            vec->readFromFile( fName.str() );

            functions_.at( i )->fromVector( *enumerator, vec, level, hyteg::All );
         }
      }
   }

   void addStateDataToStringVector( std::vector< std::string >& v )
   {
      if ( functions_.size() < memoryCapacity_ )
      {
         WALBERLA_ABORT( "Function History is not properly initialised." );
      }

      std::string historyName = functions_.at( 0 )->getFunctionName();
      std::replace( historyName.begin(), historyName.end(), ' ', '_' );

      // save number of initialsed states
      std::stringstream initState;
      initState << historyName << "_numberOfInitialisedStates " << numberOfInitialisedStates_ << ";";
      v.push_back( initState.str() );

      // save step sizes
      for ( uint_t i = 0; i < stepSizes_.size(); i++ )
      {
         std::stringstream stepSize;
         stepSize << historyName << "_stepSize_" << i << " " << std::setprecision( std::numeric_limits< real_t >::max_digits10 + 1 )
                  << *stepSizes_.at( i ) << ";";
         v.push_back( stepSize.str() );
      }
   }

   void loadStateDataFromBlockHandle( walberla::config::Config::BlockHandle& parameters )
   {
      if ( functions_.size() < memoryCapacity_ )
      {
         WALBERLA_ABORT( "Function History is not properly initialised." );
      }

      std::string historyName = functions_.at( 0 )->getFunctionName();
      std::replace( historyName.begin(), historyName.end(), ' ', '_' );

      // load number of initialsed states
      std::stringstream initState;
      initState << historyName << "_numberOfInitialisedStates";

      numberOfInitialisedStates_ = parameters.getParameter< uint_t >( initState.str() );

      // load step sizes
      for ( uint_t i = 0; i < stepSizes_.size(); i++ )
      {
         std::stringstream stepSize;
         stepSize << historyName << "_stepSize_" << i;

         *stepSizes_.at( i ) = parameters.getParameter< real_t >( stepSize.str() );
      }
   }

#ifdef HYTEG_BUILD_WITH_ADIOS2
   void registerFunctionsToAdiosCheckpointExporter( const std::shared_ptr< hyteg::AdiosCheckpointExporter >& exporter )
   {
      if ( functions_.size() < memoryCapacity_ )
      {
         WALBERLA_ABORT( "Function History is not properly initialised." );
      }
      for ( uint_t i = 0; i < memoryCapacity_; i++ )
      {
         exporter->registerFunction( *functions_.at( i ), minLevel_, maxLevel_ );
      }
   }

   void loadFunctionsFromAdiosCheckpointImporter( const std::shared_ptr< hyteg::AdiosCheckpointImporter >& importer,
                                                  uint_t                                                   adiosStep = 0 )
   {
      if ( functions_.size() < memoryCapacity_ )
      {
         WALBERLA_ABORT( "Function History is not properly initialised." );
      }
      for ( uint_t i = 0; i < memoryCapacity_; i++ )
      {
         if constexpr ( std::is_same_v< FunctionType, hyteg::P2P1TaylorHoodFunction< typename FunctionType::valueType > > )
         {
            bool successUVW = importer->restoreFunction( functions_.at( i )->uvw(), minLevel_, maxLevel_, adiosStep, false );
            if ( !successUVW )
            {
               WALBERLA_LOG_WARNING_ON_ROOT( "Warning: Function"
                                             << functions_.at( i )->uvw().getFunctionName()
                                             << " could not be loaded from checkpoint. This does not have to be a problem." );
            }

            bool successP = importer->restoreFunction( functions_.at( i )->p(), minLevel_, maxLevel_, adiosStep, false );
            if ( !successP )
            {
               WALBERLA_LOG_WARNING_ON_ROOT( "Warning: Function"
                                             << functions_.at( i )->p().getFunctionName()
                                             << " could not be loaded from checkpoint. This does not have to be a problem." );
            }
         }
         else
         {
            bool success = importer->restoreFunction( *functions_.at( i ), minLevel_, maxLevel_, adiosStep, false );
            if ( !success )
            {
               WALBERLA_LOG_WARNING_ON_ROOT( "Warning: Function"
                                             << functions_.at( i )->getFunctionName()
                                             << " could not be loaded from checkpoint. This does not have to be a problem." );
            }
         }
      }
   }

   void addStateDataToUserAttributes( std::vector< std::string >& names, std::vector< hyteg::adiosHelpers::adiostype_t >& values )
   {
      if ( functions_.size() < memoryCapacity_ )
      {
         WALBERLA_ABORT( "Function History is not properly initialised." );
      }

      // save number of initialsed states
      std::stringstream initStateName;
      initStateName << functions_.at( 0 )->getFunctionName() << "_numberOfInitialisedStates";

      names.push_back( initStateName.str() );
      values.push_back( numberOfInitialisedStates_ );

      // save step sizes
      std::stringstream stepSizesName;
      stepSizesName << functions_.at( 0 )->getFunctionName() << "_stepSizes";
      std::vector< real_t > steps;
      for ( uint_t i = 0; i < stepSizes_.size(); i++ )
      {
         steps.push_back( *stepSizes_.at( i ) );
      }

      names.push_back( stepSizesName.str() );
      values.push_back( steps );
   }

   void loadStateDataFromUserAttributes( const std::shared_ptr< hyteg::AdiosCheckpointImporter >& importer )
   {
      if ( functions_.size() < memoryCapacity_ )
      {
         WALBERLA_ABORT( "Function History is not properly initialised." );
      }

      // load number of initialsed states
      std::stringstream initStateName;
      initStateName << functions_.at( 0 )->getFunctionName() << "_numberOfInitialisedStates";

      numberOfInitialisedStates_ = importer->getUserAttributeValue< uint_t >( initStateName.str() );

      // load step sizes
      std::stringstream stepSizesName;
      stepSizesName << functions_.at( 0 )->getFunctionName() << "_stepSizes";
      std::vector< real_t > steps = importer->getUserAttributeValue< std::vector< real_t > >( stepSizesName.str() );

      uint_t maxStep = std::min( stepSizes_.size(), steps.size() );
      for ( uint_t i = 0; i < maxStep; i++ )
      {
         *stepSizes_.at( i ) = steps.at( i );
      }
   }
#else
   template < class T >
   void registerFunctionsToAdiosCheckpointExporter( const T& exporter )
   {
      WALBERLA_UNUSED( exporter );
      WALBERLA_ABORT( "ADIOS 2 unavailable. Try HYTEG_BUILD_WITH_ADIOS2=ON and defining ADIOS2_DIR." );
   }

   template < class T >
   void loadFunctionsFromAdiosCheckpointImporter( const T& importer, uint_t adiosStep = 0 )
   {
      WALBERLA_UNUSED( importer );
      WALBERLA_UNUSED( adiosStep );
      WALBERLA_ABORT( "ADIOS 2 unavailable. Try HYTEG_BUILD_WITH_ADIOS2=ON and defining ADIOS2_DIR." );
   }

   template < class T, class T2 >
   void addStateDataToUserAttributes( const T& names, const T2& values )
   {
      WALBERLA_UNUSED( names );
      WALBERLA_UNUSED( values );
      WALBERLA_ABORT( "ADIOS 2 unavailable. Try HYTEG_BUILD_WITH_ADIOS2=ON and defining ADIOS2_DIR." );
   }

   template < class T >
   void loadStateDataFromUserAttributes( const T& importer )
   {
      WALBERLA_UNUSED( importer );
      WALBERLA_ABORT( "ADIOS 2 unavailable. Try HYTEG_BUILD_WITH_ADIOS2=ON and defining ADIOS2_DIR." );
   }
#endif

   void setFunctionByIndex( uint_t index, const std::shared_ptr< FunctionType >& f )
   {
      functions_.at( index ) = f;
   }
   FunctionType& getFunctionByIndex( uint_t index )
   {
      return *functions_.at( index );
   }
   std::shared_ptr< FunctionType > getFunctionPtrByIndex( uint_t index )
   {
      return functions_.at( index );
   }
   std::shared_ptr< AdditionalDataType > getAdditionalDataPtrByIndex( uint_t index )
   {
      return additionalData_.at( index );
   }

   FunctionType&
       newState( real_t stepSize = real_c( 0 ), AdditionalDataType additionalData = AdditionalDataType(), bool useSwap = false )
   {
      if ( numberOfInitialisedStates_ < memoryCapacity_ )
      {
         numberOfInitialisedStates_++;
      }

      for ( uint_t i = memoryCapacity_ - 1; i > 0; i-- )
      {
         uint_t swapIndex = i - 1;
         if ( swapIndex >= memoryCapacity_ )
         {
            WALBERLA_ABORT( "Index error in function newState!" )
         }

         for ( uint_t l = minLevel_; l <= maxLevel_; l++ )
         {
            if ( useSwap )
            {
               functions_.at( i )->swap( *functions_.at( swapIndex ), l, hyteg::All );
            }
            else
            {
               functions_.at( i )->assign( { real_c( 1 ) }, { *functions_.at( swapIndex ) }, l, hyteg::All );
            }
         }

         *( additionalData_.at( i ) ) = *( additionalData_.at( swapIndex ) );
         *( stepSizes_.at( i ) )      = *( stepSizes_.at( swapIndex ) );
      }

      *additionalData_.at( 0 ) = additionalData;
      *stepSizes_.at( 0 )      = stepSize;

      return *functions_.at( 0 );
   }

   uint_t getNumberOfInitialisedStates() const
   {
      return numberOfInitialisedStates_;
   }
   uint_t getMemoryCapacity() const
   {
      return memoryCapacity_;
   }

   bool stateOffsetAvailable( uint_t offset ) const
   {
      return offset < numberOfInitialisedStates_;
   }

   FunctionType& getState( uint_t offset ) const
   {
      if ( offset >= memoryCapacity_ )
      {
         WALBERLA_ABORT( "Requested offset is not within the memory capacity limit!" )
      }
      if ( !stateOffsetAvailable( offset ) )
      {
         WALBERLA_ABORT( "Requested offset not initialised yet. Check with stateOffsetAvailable first!" )
      }

      return *functions_.at( offset );
   }

   std::shared_ptr< FunctionType > getStatePtr( uint_t offset ) const
   {
      if ( offset >= memoryCapacity_ )
      {
         WALBERLA_ABORT( "Requested offset is not within the memory capacity limit!" )
      }
      if ( !stateOffsetAvailable( offset ) )
      {
         WALBERLA_ABORT( "Requested offset not initialised yet. Check with stateOffsetAvailable first!" )
      }

      return functions_.at( offset );
   }

   AdditionalDataType& getStateAdditionalData( uint_t offset ) const
   {
      if ( offset >= memoryCapacity_ )
      {
         WALBERLA_ABORT( "Requested offset is not within the memory capacity limit!" )
      }
      if ( !stateOffsetAvailable( offset ) )
      {
         WALBERLA_ABORT( "Requested offset not initialised yet. Check with stateOffsetAvailable first!" )
      }

      return *additionalData_.at( offset );
   }

   std::shared_ptr< AdditionalDataType > getStateAdditionalDataPtr( uint_t offset ) const
   {
      if ( offset >= memoryCapacity_ )
      {
         WALBERLA_ABORT( "Requested offset is not within the memory capacity limit!" )
      }
      if ( !stateOffsetAvailable( offset ) )
      {
         WALBERLA_ABORT( "Requested offset not initialised yet. Check with stateOffsetAvailable first!" )
      }

      return additionalData_.at( offset );
   }

   real_t& getStateStepSize( uint_t offset ) const
   {
      if ( offset >= memoryCapacity_ )
      {
         WALBERLA_ABORT( "Requested offset is not within the memory capacity limit!" )
      }
      if ( !stateOffsetAvailable( offset ) )
      {
         WALBERLA_ABORT( "Requested offset not initialised yet. Check with stateOffsetAvailable first!" )
      }

      return *stepSizes_.at( offset );
   }

   std::shared_ptr< real_t > getStateStepSizePtr( uint_t offset ) const
   {
      if ( offset >= memoryCapacity_ )
      {
         WALBERLA_ABORT( "Requested offset is not within the memory capacity limit!" )
      }
      if ( !stateOffsetAvailable( offset ) )
      {
         WALBERLA_ABORT( "Requested offset not initialised yet. Check with stateOffsetAvailable first!" )
      }

      return stepSizes_.at( offset );
   }

   // Given Functions F0, F1, F2, ... at time points t0 >= t1 >= t2 saved as the
   // states of the FunctionHistory Hist this function expects the time step sizes
   // dt0 = t0 - t1, dt1 = t1 - t2, ... to be saved in the stepSizes vector of
   // the Function History (F0 holds dt0, F1 holds dt1, ...).
   // Sets dst to a backwards differentiation approximation of dF/dt at time t0
   // exact up to order order.
   // If not enough past states are available then the function automatically
   // returns an approximation of lower order or aborts if the history is empty.
   void derivativeApproximation( uint_t                  order,
                                 const FunctionType&     dst,
                                 const uint_t            level,
                                 const hyteg::DoFType    flag,
                                 const hyteg::UpdateType updateType = hyteg::UpdateType::Replace )
   {
      uint_t numberOfInitialisedStates = getNumberOfInitialisedStates();
      if ( numberOfInitialisedStates == 0 )
      {
         WALBERLA_ABORT( "Empty function history!" );
      }
      uint_t maximumOrder = numberOfInitialisedStates - 1;
      uint_t order_       = std::min( order, maximumOrder );

      switch ( order_ )
      {
      case 0:
         break;
      case 1: {
         const real_t dt0 = getStateStepSize( 0 );

         const real_t scaling0 = real_c( 1 );
         const real_t scaling1 = real_c( 1 );

         if ( updateType == hyteg::UpdateType::Replace )
         {
            dst.assign( { scaling0 / dt0, -scaling1 / dt0 }, { getState( 0 ), getState( 1 ) }, level, flag );
         }
         else
         {
            dst.assign( { real_c( 1 ), scaling0 / dt0, -scaling1 / dt0 }, { dst, getState( 0 ), getState( 1 ) }, level, flag );
         }
      }
      break;
      case 2: {
         const real_t dt0 = getStateStepSize( 0 );
         const real_t dt1 = getStateStepSize( 1 );

         const real_t scaling0 = real_c( 1 ) + dt0 / ( dt0 + dt1 );
         const real_t scaling1 = ( dt0 + dt1 ) / dt1;
         const real_t scaling2 = -( dt0 * dt0 ) / ( dt0 * dt1 + ( dt1 * dt1 ) );

         if ( updateType == hyteg::UpdateType::Replace )
         {
            dst.assign( { scaling0 / dt0, -scaling1 / dt0, -scaling2 / dt0 },
                        { getState( 0 ), getState( 1 ), getState( 2 ) },
                        level,
                        flag );
         }
         else
         {
            dst.assign( { real_c( 1 ), scaling0 / dt0, -scaling1 / dt0, -scaling2 / dt0 },
                        { dst, getState( 0 ), getState( 1 ), getState( 2 ) },
                        level,
                        flag );
         }
      }
      break;
      case 3: {
         const real_t dt0 = getStateStepSize( 0 );
         const real_t dt1 = getStateStepSize( 1 );
         const real_t dt2 = getStateStepSize( 2 );

         const real_t scaling0 = real_c( 1 ) + ( dt0 / ( dt0 + dt1 ) + dt0 / ( dt0 + dt1 + dt2 ) );
         const real_t scaling1 = ( ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) ) / ( dt1 * ( dt1 + dt2 ) );
         const real_t scaling2 = -( ( dt0 * dt0 * ( dt0 + dt1 + dt2 ) ) / ( dt1 * ( dt0 + dt1 ) * dt2 ) );
         const real_t scaling3 = ( dt0 * dt0 * ( dt0 + dt1 ) ) / ( dt2 * ( dt1 + dt2 ) * ( dt0 + dt1 + dt2 ) );

         if ( updateType == hyteg::UpdateType::Replace )
         {
            dst.assign( { scaling0 / dt0, -scaling1 / dt0, -scaling2 / dt0, -scaling3 / dt0 },
                        { getState( 0 ), getState( 1 ), getState( 2 ), getState( 3 ) },
                        level,
                        flag );
         }
         else
         {
            dst.assign( { real_c( 1 ), scaling0 / dt0, -scaling1 / dt0, -scaling2 / dt0, -scaling3 / dt0 },
                        { dst, getState( 0 ), getState( 1 ), getState( 2 ), getState( 3 ) },
                        level,
                        flag );
         }
      }
      break;
      case 4: {
         const real_t dt0 = getStateStepSize( 0 );
         const real_t dt1 = getStateStepSize( 1 );
         const real_t dt2 = getStateStepSize( 2 );
         const real_t dt3 = getStateStepSize( 3 );

         const real_t scaling0 =
             real_c( 1 ) + ( dt0 / ( dt0 + dt1 ) + dt0 / ( dt0 + dt1 + dt2 ) + dt0 / ( dt3 + dt0 + dt1 + dt2 ) );
         const real_t scaling1 =
             ( ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt3 + dt0 + dt1 + dt2 ) ) / ( dt1 * ( dt1 + dt2 ) * ( dt3 + dt1 + dt2 ) );
         const real_t scaling2 =
             -( ( dt0 * dt0 * ( dt0 + dt1 + dt2 ) * ( dt3 + dt0 + dt1 + dt2 ) ) / ( dt1 * ( dt0 + dt1 ) * dt2 * ( dt3 + dt2 ) ) );
         const real_t scaling3 =
             ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt3 + dt0 + dt1 + dt2 ) ) / ( dt3 * dt2 * ( dt1 + dt2 ) * ( dt0 + dt1 + dt2 ) );
         const real_t scaling4 = -( ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) ) /
                                    ( dt3 * ( dt3 + dt2 ) * ( dt3 + dt1 + dt2 ) * ( dt3 + dt0 + dt1 + dt2 ) ) );

         if ( updateType == hyteg::UpdateType::Replace )
         {
            dst.assign( { scaling0 / dt0, -scaling1 / dt0, -scaling2 / dt0, -scaling3 / dt0, -scaling4 / dt0 },
                        { getState( 0 ), getState( 1 ), getState( 2 ), getState( 3 ), getState( 4 ) },
                        level,
                        flag );
         }
         else
         {
            dst.assign( { real_c( 1 ), scaling0 / dt0, -scaling1 / dt0, -scaling2 / dt0, -scaling3 / dt0, -scaling4 / dt0 },
                        { dst, getState( 0 ), getState( 1 ), getState( 2 ), getState( 3 ), getState( 4 ) },
                        level,
                        flag );
         }
      }
      break;
      case 5: {
         const real_t dt0 = getStateStepSize( 0 );
         const real_t dt1 = getStateStepSize( 1 );
         const real_t dt2 = getStateStepSize( 2 );
         const real_t dt3 = getStateStepSize( 3 );
         const real_t dt4 = getStateStepSize( 4 );

         const real_t scaling0 = real_c( 1 ) + ( dt0 / ( dt0 + dt1 ) + dt0 / ( dt0 + dt1 + dt2 ) +
                                                 dt0 / ( dt0 + dt1 + dt2 + dt3 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 ) );
         const real_t scaling1 =
             ( ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) ) /
             ( dt1 * ( dt1 + dt2 ) * ( dt1 + dt2 + dt3 ) * ( dt1 + dt2 + dt3 + dt4 ) );
         const real_t scaling2 =
             -( ( dt0 * dt0 * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) ) /
                ( dt1 * ( dt0 + dt1 ) * dt2 * ( dt2 + dt3 ) * ( dt2 + dt3 + dt4 ) ) );
         const real_t scaling3 = ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) ) /
                                 ( dt2 * ( dt1 + dt2 ) * ( dt0 + dt1 + dt2 ) * dt3 * ( dt3 + dt4 ) );
         const real_t scaling4 = -( ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) ) /
                                    ( dt3 * ( dt2 + dt3 ) * ( dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 ) * dt4 ) );
         const real_t scaling5 =
             ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) ) /
             ( dt4 * ( dt3 + dt4 ) * ( dt2 + dt3 + dt4 ) * ( dt1 + dt2 + dt3 + dt4 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) );

         if ( updateType == hyteg::UpdateType::Replace )
         {
            dst.assign( { scaling0 / dt0, -scaling1 / dt0, -scaling2 / dt0, -scaling3 / dt0, -scaling4 / dt0, -scaling5 / dt0 },
                        { getState( 0 ), getState( 1 ), getState( 2 ), getState( 3 ), getState( 4 ), getState( 5 ) },
                        level,
                        flag );
         }
         else
         {
            dst.assign( { real_c( 1 ),
                          scaling0 / dt0,
                          -scaling1 / dt0,
                          -scaling2 / dt0,
                          -scaling3 / dt0,
                          -scaling4 / dt0,
                          -scaling5 / dt0 },
                        { dst, getState( 0 ), getState( 1 ), getState( 2 ), getState( 3 ), getState( 4 ), getState( 5 ) },
                        level,
                        flag );
         }
      }
      break;
      case 6: {
         const real_t dt0 = getStateStepSize( 0 );
         const real_t dt1 = getStateStepSize( 1 );
         const real_t dt2 = getStateStepSize( 2 );
         const real_t dt3 = getStateStepSize( 3 );
         const real_t dt4 = getStateStepSize( 4 );
         const real_t dt5 = getStateStepSize( 5 );

         const real_t scaling0 =
             real_c( 1 ) + ( dt0 / ( dt0 + dt1 ) + dt0 / ( dt0 + dt1 + dt2 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 ) +
                             dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) );
         const real_t scaling1 =
             ( ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) ) /
             ( dt1 * ( dt1 + dt2 ) * ( dt1 + dt2 + dt3 ) * ( dt1 + dt2 + dt3 + dt4 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 ) );
         const real_t scaling2 =
             -( ( dt0 * dt0 * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) ) /
                ( dt1 * ( dt0 + dt1 ) * dt2 * ( dt2 + dt3 ) * ( dt2 + dt3 + dt4 ) * ( dt2 + dt3 + dt4 + dt5 ) ) );
         const real_t scaling3 = ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                                   ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) ) /
                                 ( dt2 * ( dt1 + dt2 ) * ( dt0 + dt1 + dt2 ) * dt3 * ( dt3 + dt4 ) * ( dt3 + dt4 + dt5 ) );
         const real_t scaling4 =
             -( ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) ) /
                ( dt3 * ( dt2 + dt3 ) * ( dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 ) * dt4 * ( dt4 + dt5 ) ) );
         const real_t scaling5 =
             ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) ) /
             ( dt4 * ( dt3 + dt4 ) * ( dt2 + dt3 + dt4 ) * ( dt1 + dt2 + dt3 + dt4 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) * dt5 );
         const real_t scaling6 = -(
             ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) ) /
             ( dt5 * ( dt4 + dt5 ) * ( dt3 + dt4 + dt5 ) * ( dt2 + dt3 + dt4 + dt5 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) ) );

         if ( updateType == hyteg::UpdateType::Replace )
         {
            dst.assign(
                { scaling0 / dt0,
                  -scaling1 / dt0,
                  -scaling2 / dt0,
                  -scaling3 / dt0,
                  -scaling4 / dt0,
                  -scaling5 / dt0,
                  -scaling6 / dt0 },
                { getState( 0 ), getState( 1 ), getState( 2 ), getState( 3 ), getState( 4 ), getState( 5 ), getState( 6 ) },
                level,
                flag );
         }
         else
         {
            dst.assign(
                { real_c( 1 ),
                  scaling0 / dt0,
                  -scaling1 / dt0,
                  -scaling2 / dt0,
                  -scaling3 / dt0,
                  -scaling4 / dt0,
                  -scaling5 / dt0,
                  -scaling6 / dt0 },
                { dst, getState( 0 ), getState( 1 ), getState( 2 ), getState( 3 ), getState( 4 ), getState( 5 ), getState( 6 ) },
                level,
                flag );
         }
      }
      break;
      case 7: {
         const real_t dt0 = getStateStepSize( 0 );
         const real_t dt1 = getStateStepSize( 1 );
         const real_t dt2 = getStateStepSize( 2 );
         const real_t dt3 = getStateStepSize( 3 );
         const real_t dt4 = getStateStepSize( 4 );
         const real_t dt5 = getStateStepSize( 5 );
         const real_t dt6 = getStateStepSize( 6 );

         const real_t scaling0 =
             real_c( 1 ) + ( dt0 / ( dt0 + dt1 ) + dt0 / ( dt0 + dt1 + dt2 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 ) +
                             dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) +
                             dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) );
         const real_t scaling1 =
             ( ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) ) /
             ( dt1 * ( dt1 + dt2 ) * ( dt1 + dt2 + dt3 ) * ( dt1 + dt2 + dt3 + dt4 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 ) *
               ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) );
         const real_t scaling2 =
             -( ( dt0 * dt0 * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) ) /
                ( dt1 * ( dt0 + dt1 ) * dt2 * ( dt2 + dt3 ) * ( dt2 + dt3 + dt4 ) * ( dt2 + dt3 + dt4 + dt5 ) *
                  ( dt2 + dt3 + dt4 + dt5 + dt6 ) ) );
         const real_t scaling3 = ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                                   ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) ) /
                                 ( dt2 * ( dt1 + dt2 ) * ( dt0 + dt1 + dt2 ) * dt3 * ( dt3 + dt4 ) * ( dt3 + dt4 + dt5 ) *
                                   ( dt3 + dt4 + dt5 + dt6 ) );
         const real_t scaling4 = -( ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                                      ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) ) /
                                    ( dt3 * ( dt2 + dt3 ) * ( dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 ) * dt4 *
                                      ( dt4 + dt5 ) * ( dt4 + dt5 + dt6 ) ) );
         const real_t scaling5 = ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) *
                                   ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) ) /
                                 ( dt4 * ( dt3 + dt4 ) * ( dt2 + dt3 + dt4 ) * ( dt1 + dt2 + dt3 + dt4 ) *
                                   ( dt0 + dt1 + dt2 + dt3 + dt4 ) * dt5 * ( dt5 + dt6 ) );
         const real_t scaling6 = -( ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) *
                                      ( dt0 + dt1 + dt2 + dt3 + dt4 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) ) /
                                    ( dt5 * ( dt4 + dt5 ) * ( dt3 + dt4 + dt5 ) * ( dt2 + dt3 + dt4 + dt5 ) *
                                      ( dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * dt6 ) );
         const real_t scaling7 =
             ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) ) /
             ( dt6 * ( dt5 + dt6 ) * ( dt4 + dt5 + dt6 ) * ( dt3 + dt4 + dt5 + dt6 ) * ( dt2 + dt3 + dt4 + dt5 + dt6 ) *
               ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) );

         if ( updateType == hyteg::UpdateType::Replace )
         {
            dst.assign( { scaling0 / dt0,
                          -scaling1 / dt0,
                          -scaling2 / dt0,
                          -scaling3 / dt0,
                          -scaling4 / dt0,
                          -scaling5 / dt0,
                          -scaling6 / dt0,
                          -scaling7 / dt0 },
                        { getState( 0 ),
                          getState( 1 ),
                          getState( 2 ),
                          getState( 3 ),
                          getState( 4 ),
                          getState( 5 ),
                          getState( 6 ),
                          getState( 7 ) },
                        level,
                        flag );
         }
         else
         {
            dst.assign( { real_c( 1 ),
                          scaling0 / dt0,
                          -scaling1 / dt0,
                          -scaling2 / dt0,
                          -scaling3 / dt0,
                          -scaling4 / dt0,
                          -scaling5 / dt0,
                          -scaling6 / dt0,
                          -scaling7 / dt0 },
                        { dst,
                          getState( 0 ),
                          getState( 1 ),
                          getState( 2 ),
                          getState( 3 ),
                          getState( 4 ),
                          getState( 5 ),
                          getState( 6 ),
                          getState( 7 ) },
                        level,
                        flag );
         }
      }
      break;
      case 8: {
         const real_t dt0 = getStateStepSize( 0 );
         const real_t dt1 = getStateStepSize( 1 );
         const real_t dt2 = getStateStepSize( 2 );
         const real_t dt3 = getStateStepSize( 3 );
         const real_t dt4 = getStateStepSize( 4 );
         const real_t dt5 = getStateStepSize( 5 );
         const real_t dt6 = getStateStepSize( 6 );
         const real_t dt7 = getStateStepSize( 7 );

         const real_t scaling0 =
             real_c( 1 ) +
             ( dt0 / ( dt0 + dt1 ) + dt0 / ( dt0 + dt1 + dt2 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 ) +
               dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) +
               dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) );
         const real_t scaling1 =
             ( ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) ) /
             ( dt1 * ( dt1 + dt2 ) * ( dt1 + dt2 + dt3 ) * ( dt1 + dt2 + dt3 + dt4 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 ) *
               ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) );
         const real_t scaling2 =
             -( ( dt0 * dt0 * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) ) /
                ( dt1 * ( dt0 + dt1 ) * dt2 * ( dt2 + dt3 ) * ( dt2 + dt3 + dt4 ) * ( dt2 + dt3 + dt4 + dt5 ) *
                  ( dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) ) );
         const real_t scaling3 = ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                                   ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
                                   ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) ) /
                                 ( dt2 * ( dt1 + dt2 ) * ( dt0 + dt1 + dt2 ) * dt3 * ( dt3 + dt4 ) * ( dt3 + dt4 + dt5 ) *
                                   ( dt3 + dt4 + dt5 + dt6 ) * ( dt3 + dt4 + dt5 + dt6 + dt7 ) );
         const real_t scaling4 = -( ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                                      ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
                                      ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) ) /
                                    ( dt3 * ( dt2 + dt3 ) * ( dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 ) * dt4 *
                                      ( dt4 + dt5 ) * ( dt4 + dt5 + dt6 ) * ( dt4 + dt5 + dt6 + dt7 ) ) );
         const real_t scaling5 = ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) *
                                   ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
                                   ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) ) /
                                 ( dt4 * ( dt3 + dt4 ) * ( dt2 + dt3 + dt4 ) * ( dt1 + dt2 + dt3 + dt4 ) *
                                   ( dt0 + dt1 + dt2 + dt3 + dt4 ) * dt5 * ( dt5 + dt6 ) * ( dt5 + dt6 + dt7 ) );
         const real_t scaling6 =
             -( ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) ) /
                ( dt5 * ( dt4 + dt5 ) * ( dt3 + dt4 + dt5 ) * ( dt2 + dt3 + dt4 + dt5 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * dt6 * ( dt6 + dt7 ) ) );
         const real_t scaling7 =
             ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) ) /
             ( dt6 * ( dt5 + dt6 ) * ( dt4 + dt5 + dt6 ) * ( dt3 + dt4 + dt5 + dt6 ) * ( dt2 + dt3 + dt4 + dt5 + dt6 ) *
               ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * dt7 );
         const real_t scaling8 =
             -( ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) ) /
                ( dt7 * ( dt6 + dt7 ) * ( dt5 + dt6 + dt7 ) * ( dt4 + dt5 + dt6 + dt7 ) * ( dt3 + dt4 + dt5 + dt6 + dt7 ) *
                  ( dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) ) );

         if ( updateType == hyteg::UpdateType::Replace )
         {
            dst.assign( { scaling0 / dt0,
                          -scaling1 / dt0,
                          -scaling2 / dt0,
                          -scaling3 / dt0,
                          -scaling4 / dt0,
                          -scaling5 / dt0,
                          -scaling6 / dt0,
                          -scaling7 / dt0,
                          -scaling8 / dt0 },
                        { getState( 0 ),
                          getState( 1 ),
                          getState( 2 ),
                          getState( 3 ),
                          getState( 4 ),
                          getState( 5 ),
                          getState( 6 ),
                          getState( 7 ),
                          getState( 8 ) },
                        level,
                        flag );
         }
         else
         {
            dst.assign( { real_c( 1 ),
                          scaling0 / dt0,
                          -scaling1 / dt0,
                          -scaling2 / dt0,
                          -scaling3 / dt0,
                          -scaling4 / dt0,
                          -scaling5 / dt0,
                          -scaling6 / dt0,
                          -scaling7 / dt0,
                          -scaling8 / dt0 },
                        { dst,
                          getState( 0 ),
                          getState( 1 ),
                          getState( 2 ),
                          getState( 3 ),
                          getState( 4 ),
                          getState( 5 ),
                          getState( 6 ),
                          getState( 7 ),
                          getState( 8 ) },
                        level,
                        flag );
         }
      }
      break;
      case 9: {
         const real_t dt0 = getStateStepSize( 0 );
         const real_t dt1 = getStateStepSize( 1 );
         const real_t dt2 = getStateStepSize( 2 );
         const real_t dt3 = getStateStepSize( 3 );
         const real_t dt4 = getStateStepSize( 4 );
         const real_t dt5 = getStateStepSize( 5 );
         const real_t dt6 = getStateStepSize( 6 );
         const real_t dt7 = getStateStepSize( 7 );
         const real_t dt8 = getStateStepSize( 8 );

         const real_t scaling0 =
             real_c( 1 ) +
             ( dt0 / ( dt0 + dt1 ) + dt0 / ( dt0 + dt1 + dt2 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 ) +
               dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) +
               dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) +
               dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) );
         const real_t scaling1 =
             ( ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) /
             ( dt1 * ( dt1 + dt2 ) * ( dt1 + dt2 + dt3 ) * ( dt1 + dt2 + dt3 + dt4 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 ) *
               ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) *
               ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) );
         const real_t scaling2 =
             -( ( dt0 * dt0 * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) /
                ( dt1 * ( dt0 + dt1 ) * dt2 * ( dt2 + dt3 ) * ( dt2 + dt3 + dt4 ) * ( dt2 + dt3 + dt4 + dt5 ) *
                  ( dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) *
                  ( dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) );
         const real_t scaling3 =
             ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) /
             ( dt2 * ( dt1 + dt2 ) * ( dt0 + dt1 + dt2 ) * dt3 * ( dt3 + dt4 ) * ( dt3 + dt4 + dt5 ) * ( dt3 + dt4 + dt5 + dt6 ) *
               ( dt3 + dt4 + dt5 + dt6 + dt7 ) * ( dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) );
         const real_t scaling4 =
             -( ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) /
                ( dt3 * ( dt2 + dt3 ) * ( dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 ) * dt4 * ( dt4 + dt5 ) *
                  ( dt4 + dt5 + dt6 ) * ( dt4 + dt5 + dt6 + dt7 ) * ( dt4 + dt5 + dt6 + dt7 + dt8 ) ) );
         const real_t scaling5 =
             ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) /
             ( dt4 * ( dt3 + dt4 ) * ( dt2 + dt3 + dt4 ) * ( dt1 + dt2 + dt3 + dt4 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) * dt5 *
               ( dt5 + dt6 ) * ( dt5 + dt6 + dt7 ) * ( dt5 + dt6 + dt7 + dt8 ) );
         const real_t scaling6 =
             -( ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) /
                ( dt5 * ( dt4 + dt5 ) * ( dt3 + dt4 + dt5 ) * ( dt2 + dt3 + dt4 + dt5 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * dt6 * ( dt6 + dt7 ) * ( dt6 + dt7 + dt8 ) ) );
         const real_t scaling7 =
             ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) /
             ( dt6 * ( dt5 + dt6 ) * ( dt4 + dt5 + dt6 ) * ( dt3 + dt4 + dt5 + dt6 ) * ( dt2 + dt3 + dt4 + dt5 + dt6 ) *
               ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * dt7 * ( dt7 + dt8 ) );
         const real_t scaling8 =
             -( ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) /
                ( dt7 * ( dt6 + dt7 ) * ( dt5 + dt6 + dt7 ) * ( dt4 + dt5 + dt6 + dt7 ) * ( dt3 + dt4 + dt5 + dt6 + dt7 ) *
                  ( dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) * dt8 ) );
         const real_t scaling9 =
             ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) ) /
             ( dt8 * ( dt7 + dt8 ) * ( dt6 + dt7 + dt8 ) * ( dt5 + dt6 + dt7 + dt8 ) * ( dt4 + dt5 + dt6 + dt7 + dt8 ) *
               ( dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) * ( dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) *
               ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) );

         if ( updateType == hyteg::UpdateType::Replace )
         {
            dst.assign( { scaling0 / dt0,
                          -scaling1 / dt0,
                          -scaling2 / dt0,
                          -scaling3 / dt0,
                          -scaling4 / dt0,
                          -scaling5 / dt0,
                          -scaling6 / dt0,
                          -scaling7 / dt0,
                          -scaling8 / dt0,
                          -scaling9 / dt0 },
                        { getState( 0 ),
                          getState( 1 ),
                          getState( 2 ),
                          getState( 3 ),
                          getState( 4 ),
                          getState( 5 ),
                          getState( 6 ),
                          getState( 7 ),
                          getState( 8 ),
                          getState( 9 ) },
                        level,
                        flag );
         }
         else
         {
            dst.assign( { real_c( 1 ),
                          scaling0 / dt0,
                          -scaling1 / dt0,
                          -scaling2 / dt0,
                          -scaling3 / dt0,
                          -scaling4 / dt0,
                          -scaling5 / dt0,
                          -scaling6 / dt0,
                          -scaling7 / dt0,
                          -scaling8 / dt0,
                          -scaling9 / dt0 },
                        { dst,
                          getState( 0 ),
                          getState( 1 ),
                          getState( 2 ),
                          getState( 3 ),
                          getState( 4 ),
                          getState( 5 ),
                          getState( 6 ),
                          getState( 7 ),
                          getState( 8 ),
                          getState( 9 ) },
                        level,
                        flag );
         }
      }
      break;
      default:
         WALBERLA_ABORT( "Requested derivative approximation order not implemented!" );
         break;
      }
   }

   // Given Functions F0, F1, F2, ... at time points t0 >= t1 >= t2 saved as the
   // states of the FunctionHistory Hist this function expects the time step sizes
   // dt0 = t0 - t1, dt1 = t1 - t2, ... to be saved in the stepSizes vector of
   // the Function History (F0 holds dt0, F1 holds dt1, ...).
   // Sets dst to an extrapolation exact up to order order
   // If not enough past states are available then the function automatically
   // returns an approximation of lower order or aborts if the history is empty.
   void extrapolate( uint_t order, real_t dt, FunctionType& dst, uint_t level, hyteg::DoFType flag )
   {
      if ( numberOfInitialisedStates_ == 0 )
      {
         WALBERLA_ABORT( "Empty function history!" );
      }
      uint_t maximumOrder = numberOfInitialisedStates_ - 1;
      order               = std::min( order, maximumOrder );

      switch ( order )
      {
      case 0: {
         dst.assign( { real_c( 1.0 ) }, { getState( 0 ) }, level, flag );
      }
      break;
      case 1: {
         const real_t dt0      = getStateStepSize( 0 );
         const real_t scaling0 = real_c( 1.0 ) + dt / dt0;
         const real_t scaling1 = -dt / dt0;

         dst.assign( { scaling0, scaling1 }, { getState( 0 ), getState( 1 ) }, level, flag );
      }
      break;
      case 2: {
         const real_t dt0      = getStateStepSize( 0 );
         const real_t dt1      = getStateStepSize( 1 );
         const real_t scaling0 = ( dt + dt0 ) * ( dt + dt0 + dt1 ) / ( dt0 * ( dt0 + dt1 ) );
         const real_t scaling1 = -dt * ( dt + dt0 + dt1 ) / ( dt0 * dt1 );
         const real_t scaling2 = dt * ( dt0 + dt ) / ( dt1 * ( dt0 + dt1 ) );

         dst.assign( { scaling0, scaling1, scaling2 }, { getState( 0 ), getState( 1 ), getState( 2 ) }, level, flag );
      }
      break;
      case 3: {
         const real_t dt0 = getStateStepSize( 0 );
         const real_t dt1 = getStateStepSize( 1 );
         const real_t dt2 = getStateStepSize( 2 );
         const real_t scaling0 =
             ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) / ( dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) );
         const real_t scaling1 = -( dt * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) ) / ( dt0 * dt1 * ( dt1 + dt2 ) );
         const real_t scaling2 = ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 + dt2 ) ) / ( dt1 * dt2 * ( dt0 + dt1 ) );
         const real_t scaling3 = -( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 ) ) / ( dt2 * ( dt1 + dt2 ) * ( dt0 + dt1 + dt2 ) );

         dst.assign( { scaling0, scaling1, scaling2, scaling3 },
                     { getState( 0 ), getState( 1 ), getState( 2 ), getState( 3 ) },
                     level,
                     flag );
      }
      break;
      case 4: {
         const real_t dt0 = getStateStepSize( 0 );
         const real_t dt1 = getStateStepSize( 1 );
         const real_t dt2 = getStateStepSize( 2 );
         const real_t dt3 = getStateStepSize( 3 );

         const real_t scaling0 =
             ( ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) ) /
             ( dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) );
         const real_t scaling1 = -( ( dt * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) ) /
                                    ( dt0 * dt1 * ( dt1 + dt2 ) * ( dt1 + dt2 + dt3 ) ) );
         const real_t scaling2 = ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) ) /
                                 ( dt1 * ( dt0 + dt1 ) * dt2 * ( dt2 + dt3 ) );
         const real_t scaling3 = -( ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) ) /
                                    ( dt2 * ( dt1 + dt2 ) * ( dt0 + dt1 + dt2 ) * dt3 ) );
         const real_t scaling4 = ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) ) /
                                 ( dt3 * ( dt2 + dt3 ) * ( dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 ) );

         dst.assign( { scaling0, scaling1, scaling2, scaling3, scaling4 },
                     { getState( 0 ), getState( 1 ), getState( 2 ), getState( 3 ), getState( 4 ) },
                     level,
                     flag );
      }
      break;
      case 5: {
         const real_t dt0 = getStateStepSize( 0 );
         const real_t dt1 = getStateStepSize( 1 );
         const real_t dt2 = getStateStepSize( 2 );
         const real_t dt3 = getStateStepSize( 3 );
         const real_t dt4 = getStateStepSize( 4 );

         const real_t scaling0 =
             ( ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) *
               ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) ) /
             ( dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) );
         const real_t scaling1 = -( ( dt * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) *
                                      ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) ) /
                                    ( dt0 * dt1 * ( dt1 + dt2 ) * ( dt1 + dt2 + dt3 ) * ( dt1 + dt2 + dt3 + dt4 ) ) );
         const real_t scaling2 = ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) *
                                   ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) ) /
                                 ( dt1 * ( dt0 + dt1 ) * dt2 * ( dt2 + dt3 ) * ( dt2 + dt3 + dt4 ) );
         const real_t scaling3 = -(
             ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) ) /
             ( dt2 * ( dt1 + dt2 ) * ( dt0 + dt1 + dt2 ) * dt3 * ( dt3 + dt4 ) ) );
         const real_t scaling4 =
             ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) ) /
             ( dt3 * ( dt2 + dt3 ) * ( dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 ) * dt4 );
         const real_t scaling5 =
             -( ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) ) /
                ( dt4 * ( dt3 + dt4 ) * ( dt2 + dt3 + dt4 ) * ( dt1 + dt2 + dt3 + dt4 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) ) );

         dst.assign( { scaling0, scaling1, scaling2, scaling3, scaling4, scaling5 },
                     { getState( 0 ), getState( 1 ), getState( 2 ), getState( 3 ), getState( 4 ), getState( 5 ) },
                     level,
                     flag );
      }
      break;
      case 6: {
         const real_t dt0 = getStateStepSize( 0 );
         const real_t dt1 = getStateStepSize( 1 );
         const real_t dt2 = getStateStepSize( 2 );
         const real_t dt3 = getStateStepSize( 3 );
         const real_t dt4 = getStateStepSize( 4 );
         const real_t dt5 = getStateStepSize( 5 );

         const real_t scaling0 = ( ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) *
                                   ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) ) /
                                 ( dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) *
                                   ( dt0 + dt1 + dt2 + dt3 + dt4 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) );
         const real_t scaling1 = -(
             ( dt * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) *
               ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) ) /
             ( dt0 * dt1 * ( dt1 + dt2 ) * ( dt1 + dt2 + dt3 ) * ( dt1 + dt2 + dt3 + dt4 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 ) ) );
         const real_t scaling2 = ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) *
                                   ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) ) /
                                 ( dt1 * ( dt0 + dt1 ) * dt2 * ( dt2 + dt3 ) * ( dt2 + dt3 + dt4 ) * ( dt2 + dt3 + dt4 + dt5 ) );
         const real_t scaling3 = -( ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) *
                                      ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) ) /
                                    ( dt2 * ( dt1 + dt2 ) * ( dt0 + dt1 + dt2 ) * dt3 * ( dt3 + dt4 ) * ( dt3 + dt4 + dt5 ) ) );
         const real_t scaling4 = ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) *
                                   ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) ) /
                                 ( dt3 * ( dt2 + dt3 ) * ( dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 ) * dt4 * ( dt4 + dt5 ) );
         const real_t scaling5 = -(
             ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) *
               ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) ) /
             ( dt4 * ( dt3 + dt4 ) * ( dt2 + dt3 + dt4 ) * ( dt1 + dt2 + dt3 + dt4 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) * dt5 ) );
         const real_t scaling6 = ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) *
                                   ( dt + dt0 + dt1 + dt2 + dt3 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) ) /
                                 ( dt5 * ( dt4 + dt5 ) * ( dt3 + dt4 + dt5 ) * ( dt2 + dt3 + dt4 + dt5 ) *
                                   ( dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) );

         dst.assign( { scaling0, scaling1, scaling2, scaling3, scaling4, scaling5, scaling6 },
                     { getState( 0 ), getState( 1 ), getState( 2 ), getState( 3 ), getState( 4 ), getState( 5 ), getState( 6 ) },
                     level,
                     flag );
      }
      break;
      case 7: {
         const real_t dt0 = getStateStepSize( 0 );
         const real_t dt1 = getStateStepSize( 1 );
         const real_t dt2 = getStateStepSize( 2 );
         const real_t dt3 = getStateStepSize( 3 );
         const real_t dt4 = getStateStepSize( 4 );
         const real_t dt5 = getStateStepSize( 5 );
         const real_t dt6 = getStateStepSize( 6 );

         const real_t scaling0 =
             ( ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) *
               ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) *
               ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) ) /
             ( dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) );
         const real_t scaling1 = -( ( dt * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) *
                                      ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) *
                                      ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) ) /
                                    ( dt0 * dt1 * ( dt1 + dt2 ) * ( dt1 + dt2 + dt3 ) * ( dt1 + dt2 + dt3 + dt4 ) *
                                      ( dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) ) );
         const real_t scaling2 = ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) *
                                   ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) *
                                   ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) ) /
                                 ( dt1 * ( dt0 + dt1 ) * dt2 * ( dt2 + dt3 ) * ( dt2 + dt3 + dt4 ) * ( dt2 + dt3 + dt4 + dt5 ) *
                                   ( dt2 + dt3 + dt4 + dt5 + dt6 ) );
         const real_t scaling3 =
             -( ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) *
                  ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) ) /
                ( dt2 * ( dt1 + dt2 ) * ( dt0 + dt1 + dt2 ) * dt3 * ( dt3 + dt4 ) * ( dt3 + dt4 + dt5 ) *
                  ( dt3 + dt4 + dt5 + dt6 ) ) );
         const real_t scaling4 =
             ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) *
               ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) ) /
             ( dt3 * ( dt2 + dt3 ) * ( dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 ) * dt4 * ( dt4 + dt5 ) *
               ( dt4 + dt5 + dt6 ) );
         const real_t scaling5 =
             -( ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) *
                  ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) ) /
                ( dt4 * ( dt3 + dt4 ) * ( dt2 + dt3 + dt4 ) * ( dt1 + dt2 + dt3 + dt4 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) * dt5 *
                  ( dt5 + dt6 ) ) );
         const real_t scaling6 =
             ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) *
               ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) ) /
             ( dt5 * ( dt4 + dt5 ) * ( dt3 + dt4 + dt5 ) * ( dt2 + dt3 + dt4 + dt5 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * dt6 );
         const real_t scaling7 =
             -( ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) *
                  ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) ) /
                ( dt6 * ( dt5 + dt6 ) * ( dt4 + dt5 + dt6 ) * ( dt3 + dt4 + dt5 + dt6 ) * ( dt2 + dt3 + dt4 + dt5 + dt6 ) *
                  ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) ) );

         dst.assign( { scaling0, scaling1, scaling2, scaling3, scaling4, scaling5, scaling6, scaling7 },
                     { getState( 0 ),
                       getState( 1 ),
                       getState( 2 ),
                       getState( 3 ),
                       getState( 4 ),
                       getState( 5 ),
                       getState( 6 ),
                       getState( 7 ) },
                     level,
                     flag );
      }
      break;
      case 8: {
         const real_t dt0 = getStateStepSize( 0 );
         const real_t dt1 = getStateStepSize( 1 );
         const real_t dt2 = getStateStepSize( 2 );
         const real_t dt3 = getStateStepSize( 3 );
         const real_t dt4 = getStateStepSize( 4 );
         const real_t dt5 = getStateStepSize( 5 );
         const real_t dt6 = getStateStepSize( 6 );
         const real_t dt7 = getStateStepSize( 7 );

         const real_t scaling0 =
             ( ( dt0 + dt ) * ( dt0 + dt1 + dt ) * ( dt0 + dt1 + dt2 + dt ) * ( dt0 + dt1 + dt2 + dt3 + dt ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt ) ) /
             ( dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) );
         const real_t scaling1 =
             -( ( dt * ( dt0 + dt1 + dt ) * ( dt0 + dt1 + dt2 + dt ) * ( dt0 + dt1 + dt2 + dt3 + dt ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt ) ) /
                ( dt0 * dt1 * ( dt1 + dt2 ) * ( dt1 + dt2 + dt3 ) * ( dt1 + dt2 + dt3 + dt4 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 ) *
                  ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) ) );
         const real_t scaling2 =
             ( dt * ( dt0 + dt ) * ( dt0 + dt1 + dt2 + dt ) * ( dt0 + dt1 + dt2 + dt3 + dt ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt ) ) /
             ( dt1 * ( dt0 + dt1 ) * dt2 * ( dt2 + dt3 ) * ( dt2 + dt3 + dt4 ) * ( dt2 + dt3 + dt4 + dt5 ) *
               ( dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) );
         const real_t scaling3 =
             -( ( dt * ( dt0 + dt ) * ( dt0 + dt1 + dt ) * ( dt0 + dt1 + dt2 + dt3 + dt ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt ) ) /
                ( dt2 * ( dt1 + dt2 ) * ( dt0 + dt1 + dt2 ) * dt3 * ( dt3 + dt4 ) * ( dt3 + dt4 + dt5 ) *
                  ( dt3 + dt4 + dt5 + dt6 ) * ( dt3 + dt4 + dt5 + dt6 + dt7 ) ) );
         const real_t scaling4 =
             ( dt * ( dt0 + dt ) * ( dt0 + dt1 + dt ) * ( dt0 + dt1 + dt2 + dt ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt ) ) /
             ( dt3 * ( dt2 + dt3 ) * ( dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 ) * dt4 * ( dt4 + dt5 ) * ( dt4 + dt5 + dt6 ) *
               ( dt4 + dt5 + dt6 + dt7 ) );
         const real_t scaling5 =
             -( ( dt * ( dt0 + dt ) * ( dt0 + dt1 + dt ) * ( dt0 + dt1 + dt2 + dt ) * ( dt0 + dt1 + dt2 + dt3 + dt ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt ) ) /
                ( dt4 * ( dt3 + dt4 ) * ( dt2 + dt3 + dt4 ) * ( dt1 + dt2 + dt3 + dt4 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) * dt5 *
                  ( dt5 + dt6 ) * ( dt5 + dt6 + dt7 ) ) );
         const real_t scaling6 =
             ( dt * ( dt0 + dt ) * ( dt0 + dt1 + dt ) * ( dt0 + dt1 + dt2 + dt ) * ( dt0 + dt1 + dt2 + dt3 + dt ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt ) ) /
             ( dt5 * ( dt4 + dt5 ) * ( dt3 + dt4 + dt5 ) * ( dt2 + dt3 + dt4 + dt5 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * dt6 * ( dt6 + dt7 ) );
         const real_t scaling7 =
             -( ( dt * ( dt0 + dt ) * ( dt0 + dt1 + dt ) * ( dt0 + dt1 + dt2 + dt ) * ( dt0 + dt1 + dt2 + dt3 + dt ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt ) ) /
                ( dt6 * ( dt5 + dt6 ) * ( dt4 + dt5 + dt6 ) * ( dt3 + dt4 + dt5 + dt6 ) * ( dt2 + dt3 + dt4 + dt5 + dt6 ) *
                  ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * dt7 ) );
         const real_t scaling8 =
             ( dt * ( dt0 + dt ) * ( dt0 + dt1 + dt ) * ( dt0 + dt1 + dt2 + dt ) * ( dt0 + dt1 + dt2 + dt3 + dt ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt ) ) /
             ( dt7 * ( dt6 + dt7 ) * ( dt5 + dt6 + dt7 ) * ( dt4 + dt5 + dt6 + dt7 ) * ( dt3 + dt4 + dt5 + dt6 + dt7 ) *
               ( dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) );

         dst.assign( { scaling0, scaling1, scaling2, scaling3, scaling4, scaling5, scaling6, scaling7, scaling8 },
                     { getState( 0 ),
                       getState( 1 ),
                       getState( 2 ),
                       getState( 3 ),
                       getState( 4 ),
                       getState( 5 ),
                       getState( 6 ),
                       getState( 7 ),
                       getState( 8 ) },
                     level,
                     flag );
      }
      break;
      case 9: {
         const real_t dt0 = getStateStepSize( 0 );
         const real_t dt1 = getStateStepSize( 1 );
         const real_t dt2 = getStateStepSize( 2 );
         const real_t dt3 = getStateStepSize( 3 );
         const real_t dt4 = getStateStepSize( 4 );
         const real_t dt5 = getStateStepSize( 5 );
         const real_t dt6 = getStateStepSize( 6 );
         const real_t dt7 = getStateStepSize( 7 );
         const real_t dt8 = getStateStepSize( 8 );

         const real_t scaling0 =
             ( ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) *
               ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) *
               ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) *
               ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) /
             ( dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) );
         const real_t scaling1 =
             -( ( dt * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) *
                  ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) *
                  ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) *
                  ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) /
                ( dt0 * dt1 * ( dt1 + dt2 ) * ( dt1 + dt2 + dt3 ) * ( dt1 + dt2 + dt3 + dt4 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 ) *
                  ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) *
                  ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) );
         const real_t scaling2 =
             ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) *
               ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) *
               ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) *
               ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) /
             ( dt1 * ( dt0 + dt1 ) * dt2 * ( dt2 + dt3 ) * ( dt2 + dt3 + dt4 ) * ( dt2 + dt3 + dt4 + dt5 ) *
               ( dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) *
               ( dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) );
         const real_t scaling3 =
             -( ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) *
                  ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
                  ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) *
                  ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) /
                ( dt2 * ( dt1 + dt2 ) * ( dt0 + dt1 + dt2 ) * dt3 * ( dt3 + dt4 ) * ( dt3 + dt4 + dt5 ) *
                  ( dt3 + dt4 + dt5 + dt6 ) * ( dt3 + dt4 + dt5 + dt6 + dt7 ) * ( dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) );
         const real_t scaling4 =
             ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) *
               ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
               ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) *
               ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) /
             ( dt3 * ( dt2 + dt3 ) * ( dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 ) * dt4 * ( dt4 + dt5 ) * ( dt4 + dt5 + dt6 ) *
               ( dt4 + dt5 + dt6 + dt7 ) * ( dt4 + dt5 + dt6 + dt7 + dt8 ) );
         const real_t scaling5 =
             -( ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) *
                  ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
                  ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) *
                  ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) /
                ( dt4 * ( dt3 + dt4 ) * ( dt2 + dt3 + dt4 ) * ( dt1 + dt2 + dt3 + dt4 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) * dt5 *
                  ( dt5 + dt6 ) * ( dt5 + dt6 + dt7 ) * ( dt5 + dt6 + dt7 + dt8 ) ) );
         const real_t scaling6 =
             ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) *
               ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
               ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) *
               ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) /
             ( dt5 * ( dt4 + dt5 ) * ( dt3 + dt4 + dt5 ) * ( dt2 + dt3 + dt4 + dt5 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * dt6 * ( dt6 + dt7 ) * ( dt6 + dt7 + dt8 ) );
         const real_t scaling7 =
             -( ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) *
                  ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) *
                  ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) *
                  ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) /
                ( dt6 * ( dt5 + dt6 ) * ( dt4 + dt5 + dt6 ) * ( dt3 + dt4 + dt5 + dt6 ) * ( dt2 + dt3 + dt4 + dt5 + dt6 ) *
                  ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * dt7 * ( dt7 + dt8 ) ) );
         const real_t scaling8 =
             ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) *
               ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) *
               ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) /
             ( dt7 * ( dt6 + dt7 ) * ( dt5 + dt6 + dt7 ) * ( dt4 + dt5 + dt6 + dt7 ) * ( dt3 + dt4 + dt5 + dt6 + dt7 ) *
               ( dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) * dt8 );
         const real_t scaling9 =
             -( ( dt * ( dt + dt0 ) * ( dt + dt0 + dt1 ) * ( dt + dt0 + dt1 + dt2 ) * ( dt + dt0 + dt1 + dt2 + dt3 ) *
                  ( dt + dt0 + dt1 + dt2 + dt3 + dt4 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) *
                  ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt + dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) ) /
                ( dt8 * ( dt7 + dt8 ) * ( dt6 + dt7 + dt8 ) * ( dt5 + dt6 + dt7 + dt8 ) * ( dt4 + dt5 + dt6 + dt7 + dt8 ) *
                  ( dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) * ( dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) *
                  ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) );

         dst.assign( { scaling0, scaling1, scaling2, scaling3, scaling4, scaling5, scaling6, scaling7, scaling8, scaling9 },
                     { getState( 0 ),
                       getState( 1 ),
                       getState( 2 ),
                       getState( 3 ),
                       getState( 4 ),
                       getState( 5 ),
                       getState( 6 ),
                       getState( 7 ),
                       getState( 8 ),
                       getState( 9 ) },
                     level,
                     flag );
      }
      break;
      default:
         WALBERLA_ABORT( "Requested extrapolation order not implemented!" );
         break;
      }
   }

 private:
   uint_t                                         memoryCapacity_; // maximum number of saved functions
   uint_t                                         minLevel_;       // minimum function level
   uint_t                                         maxLevel_;       // maximum function level
   std::vector< std::shared_ptr< FunctionType > > functions_;      // vector containing all function handles
   std::vector< std::shared_ptr< real_t > >       stepSizes_;      // vector containing the step sizes between states
   std::vector< std::shared_ptr< AdditionalDataType > >
          additionalData_;            // vector containing additional data that you saved for the individual states
   uint_t numberOfInitialisedStates_; // current number of initialised states (maximum memoryCapacity_)
};

} // namespace hyteg