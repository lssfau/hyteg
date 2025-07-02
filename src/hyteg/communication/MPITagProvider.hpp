/*
 * Copyright (c) 2024-2025 Marcus Mohr, Andreas Burkhart.
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

#include <atomic>
#include <limits>

#include "core/logging/Logging.h"

namespace hyteg {
namespace communication {

/// Auxilliary class to provide MPI tag values for BufferedCommunicator and other MPI related objects/function calls
class MPITagProvider
{
   /// Maximal tag value supported by the MPI library implementation used
   ///
   /// The MPI 4.1 standard requires the largest tag to be at least 2^15-1 = 32,767.
   /// Larger values are possible, though. Since the value must be an int, it cannot
   /// exceed 2,147,483,647 for the standard 32-bit signed int setting.
   static int maxMPITag_;

   /// Stores the next tag that will be returned by getMPITag()
   static std::atomic_int nextMPITag_;

   /// Marks whether class can still provide tag values
   static bool poolExhausted_;

   /// Pool of returned MPI tags
   static std::vector< int > returnedTags_;

 public:
   /// Return the largest possible tag value supported by the MPI library in use
   static int getMaxMPITag()
   {
#ifdef HYTEG_BUILD_WITH_MPI
      void* maxTag;
      int   status;
      MPI_Comm_get_attr( walberla::mpi::MPIManager::instance()->comm(), MPI_TAG_UB, &maxTag, &status );
      if ( status == 0 )
      {
         WALBERLA_ABORT( "Failed to query maximal tag value from MPI implementation!" );
      }
      return *static_cast< int* >( maxTag );
#else
      return std::numeric_limits< int >::max();
#endif
   }

   /// Returns an MPI Tag back to the pool
   /// Make sure that the tag is no longer referenced before returning!
   static void returnMPITag( int returnedTag )
   {
      returnedTags_.push_back( returnedTag );
   }

   /// Return another MPI tag value
   ///
   /// The current implementation is very simple. In order to return unique tag values it starts with
   /// the smallest possbile value, i.e. 0, and then returns tags by incrementation until reaching the
   /// limit. Thus, there is no re-use of values that are no longer needed, and the pool of tags might
   /// get exhausted. In this case the class calls WALBERLA_ABORT().
   static int getMPITag()
   {
      // initialise largest available tag value (can only happen once MPI was activated)
      if ( maxMPITag_ == 0u )
      {
         maxMPITag_ = getMaxMPITag();
      }

      if ( poolExhausted_ )
      {
         WALBERLA_ABORT( "Your application exhausted the pool of available MPI tags.\n"
                         << "Your MPI implementation provides a maximum of " << maxMPITag_ << " tags." );
      }

      if ( returnedTags_.size() > 0 )
      {
         int recylcedTag = returnedTags_.back();
         returnedTags_.pop_back();
         return recylcedTag;
      }

      // this will store the return tag
      int freshTag{ nextMPITag_ };

      // check whether we can safely increase the tag
      if ( nextMPITag_ == maxMPITag_ )
      {
         poolExhausted_ = true;
      }
      else
      {
         ++nextMPITag_;
      }

      return freshTag;
   }
};

} // namespace communication
} // namespace hyteg
