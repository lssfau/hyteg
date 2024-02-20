/*
 * Copyright (c) 2023 Hamish Brown.
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

#include <cmath>
#include <vector>

#include "core/DataTypes.h"

using namespace hyteg;

namespace terraneo {
//! Class providing methods for working with radial profiles (e.g. mean / max / min )
//! Note: Currently only implemented for P2 Scalar functions on spherical meshes with evenly distributed radial layers
template < typename FunctionType >
class RadialProfileTool
{
 public:
   RadialProfileTool( const FunctionType& u, FunctionType& tmp, real_t rMax, real_t rMin, uint_t nRad, uint_t level )
   : u_( u )
   , tmp_( tmp )
   , rMax_( rMax )
   , rMin_( rMin )
   , nRad_( nRad )
   , level_( level )
   {
      numberOfLayers_ = 2 * ( nRad_ - 1 ) * ( levelinfo::num_microvertices_per_edge( level_ ) - 1 );
      //set radius and counter values, constant during lifetime of object
      setRadii_();
      setCounter_();
   }

   std::vector< real_t > getRadii() { return shellRadii_; };

   std::vector< real_t > getMeanProfile()
   {
      setMeanProfile_();
      return meanProfile_;
   };

   std::vector< real_t > getMaxProfile()
   {
      setMaxProfile_();
      return maxProfile_;
   };

   std::vector< real_t > getMinProfile()
   {
      setMinProfile_();
      return minProfile_;
   };

   void logProfilesToFile( std::string fileName, std::string fieldName )
   {
      //update all profiles prior to logging
      setMeanProfile_();
      setMaxProfile_();
      setMinProfile_();

      WALBERLA_ROOT_SECTION()
      {
         std::ofstream outFile( fileName );

         if ( outFile.fail() )
         {
            WALBERLA_ABORT( "Failed to open file \"" << fileName << "\" for logging radial profiles. " )
         }

         outFile << std::string( "Radius  " ) << fieldName << std::string( "_Mean " ) << fieldName << std::string( "_Max " )
                 << fieldName << std::string( "_Min \n" );
         for ( uint_t shell = 0; shell < numberOfLayers_ + 1; ++shell )
         {
            outFile << walberla::format( "%6.4f  %7.4f  %6.4f  %6.4f \n",
                                         shellRadii_.at( shell ),
                                         meanProfile_.at( shell ),
                                         maxProfile_.at( shell ),
                                         minProfile_.at( shell ) );
         }
         outFile.close();

         if ( outFile.is_open() )
         {
            WALBERLA_ABORT( "Failed to close file \"" << fileName << "\" when logging radial profiles. " )
         }
      }
   };

 private:
   uint_t                numberOfLayers_;
   real_t                rMin_, rMax_;
   uint_t                nRad_;
   uint_t                level_;
   std::vector< real_t > shellRadii_;
   std::vector< real_t > meanProfile_;
   std::vector< real_t > maxProfile_;
   std::vector< real_t > minProfile_;
   std::vector< uint_t > counter_;

   const FunctionType& u_;
   FunctionType&       tmp_;

   //get radius based on shell number
   real_t getRadius_( const uint_t& shell ) { return rMin_ + real_c( shell ) * ( rMax_ - rMin_ ) / real_c( numberOfLayers_ ); };

   //get shell number based on radius
   uint_t getShell_( const real_t& radius )
   {
      return static_cast< uint_t >( std::round( real_c( numberOfLayers_ ) * ( ( radius - rMin_ ) / ( rMax_ - rMin_ ) ) ) );
   };

   //fill shellRadii
   void setRadii_()
   {
      shellRadii_.reserve( numberOfLayers_ + 1 );

      for ( uint_t shell = 0; shell < numberOfLayers_ + 1; ++shell )
      {
         shellRadii_.push_back( getRadius_( shell ) );
      }
   };

   void setCounter_()
   {
      counter_.reserve( numberOfLayers_ + 1 );

      for ( uint_t shell = 0; shell < numberOfLayers_ + 1; ++shell )
      {
         counter_.push_back( 0 );
      }

      std::function< real_t( const Point3D& ) > gatherDoFs = [&]( const Point3D& x ) {
         real_t radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );

         uint_t shell = getShell_( radius );

         //manual bounds checking
         WALBERLA_ASSERT( shell < counter_.size() );

         //increment counter
         ++counter_.at( shell );
         // something to interpolate:
         return real_c( 0 );
      };

      tmp_.interpolate( gatherDoFs, level_, All );
      walberla::mpi::allReduceInplace( counter_, walberla::mpi::SUM );
   }

   void setMeanProfile_()
   {
      if ( meanProfile_.empty() )
      {
         meanProfile_ = std::vector< real_t >( numberOfLayers_ + 1, real_c( 0 ) );         
      }

      else
      {
         std::fill( meanProfile_.begin(), meanProfile_.end(), real_c( 0 ) );
      }

      std::function< real_t( const Point3D&, const std::vector< real_t >& ) > gatherValues =
          [&]( const Point3D& x, const std::vector< real_t >& Values ) {
             real_t radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );

             uint_t shell = getShell_( radius );

             //manual bounds checking
             WALBERLA_ASSERT( shell < meanProfile_.size() );

             //add function value to corresponding row in profile (using shell number)
             meanProfile_.at( shell ) += Values[0];

             // something to interpolate:
             return real_c( 0 );
          };

      //interpolate is used to cycle through all DoFs on a process and fill relevant parts of profile with total temperature and number of DoFs

      tmp_.interpolate( gatherValues, { u_ }, level_, All );

      //sum total values on each shell over all processes
      walberla::mpi::allReduceInplace( meanProfile_, walberla::mpi::SUM );

      //now find mean with total / counter
      for ( uint_t shell = 0; shell < numberOfLayers_ + 1; ++shell )
      {
         meanProfile_.at( shell ) /= real_c( counter_.at( shell ) );
      }
   };

   void setMaxProfile_()
   {
      real_t absMin = -std::numeric_limits< real_t >::max();

      //we start with the minimum possible value of a real_t and check whether the values on a shell are larger

      if ( maxProfile_.empty() )
      {
         maxProfile_ = std::vector< real_t >( numberOfLayers_ + 1, absMin );
      }

      else
      {
         std::fill( maxProfile_.begin(), maxProfile_.end(), absMin );
      }

      std::function< real_t( const Point3D&, const std::vector< real_t >& ) > gatherMaxValues =
          [&]( const Point3D& x, const std::vector< real_t >& Values ) {
             real_t radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );

             uint_t shell = getShell_( radius );

             //manual bounds checking
             WALBERLA_ASSERT( shell < maxProfile_.size() );

             // check whether value is larger than current value on shell, and update if so
             if ( Values[0] > maxProfile_.at( shell ) )
             {
                maxProfile_.at( shell ) = Values[0];
             }

             //need to return something to interpolate:
             return real_c( 0 );
          };

      //interpolate is used to cycle through all DoFs on a process and fill relevant parts of profile with max value on a process

      tmp_.interpolate( gatherMaxValues, { u_ }, level_, All );

      //get max profile over all processes
      walberla::mpi::allReduceInplace( maxProfile_, walberla::mpi::MAX );
   };

   void setMinProfile_()
   {
      real_t absMax = std::numeric_limits< real_t >::max();

      //we start with the maximum possible value and check whether the values on a shell are smaller
      if ( minProfile_.empty() )
      {
         minProfile_ = std::vector< real_t >( numberOfLayers_ + 1, absMax );
      }

      else
      {
         std::fill( minProfile_.begin(), minProfile_.end(), absMax );
      }

      std::function< real_t( const Point3D&, const std::vector< real_t >& ) > gatherMinValues =
          [&]( const Point3D& x, const std::vector< real_t >& Values ) {
             real_t radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );

             uint_t shell = getShell_( radius );

             //manual bounds checking
             WALBERLA_ASSERT( shell < minProfile_.size() );

             // check whether value is smaller than current value on shell, and update if so
             if ( Values[0] < minProfile_.at( shell ) )
             {
                minProfile_.at( shell ) = Values[0];
             }

             //need to return something to interpolate:
             return real_c( 0 );
          };

      //interpolate is used to cycle through all DoFs on a process and fill relevant parts of profile with min value on a process

      tmp_.interpolate( gatherMinValues, { u_ }, level_, All );

      //get min profile over all processes
      walberla::mpi::allReduceInplace( minProfile_, walberla::mpi::MIN );
   };
};
} //namespace terraneo