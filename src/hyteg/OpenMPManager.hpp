/*
 * Copyright (c) 2017-2020 Nils Kohl.
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

#include "core/singleton/Singleton.h"
#include "core/OpenMP.h"

namespace hyteg {

/// \brief Class that may be useful to globally manage OpenMP functionality.
class OpenMPManager : public walberla::singleton::Singleton< OpenMPManager >
{
 public:

   WALBERLA_BEFRIEND_SINGLETON;

   /// \brief Forces the number of OpenMP threads to 1 after this call.
   void forceSerial()
   {
      #ifdef HYTEG_BUILD_WITH_OPENMP
      omp_set_num_threads(1);
      #endif
   }

   /// \brief Resets the (maximum) number of OpenMP threads to the number of threads that OpenMP was configured to
   ///        during construction of this manager.
   void resetToParallel()
   {
      #ifdef HYTEG_BUILD_WITH_OPENMP
      omp_set_num_threads(maxNumThreads_);
      #endif
   }

   /// \brief Returns the current (maximum) number of threads.
   int numThreads() const
   {
      #ifdef HYTEG_BUILD_WITH_OPENMP
      return omp_get_max_threads();
      #else
      return 1;
      #endif
   }

 private:

   OpenMPManager()
   : maxNumThreads_(
#ifdef HYTEG_BUILD_WITH_OPENMP
         omp_get_max_threads()
#else
         1
#endif
         ) {}

   int maxNumThreads_;

};

}