/*
 * Copyright (c) 2017-2019 Nils Kohl.
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

#include "core/logging/Logging.h"

namespace hyteg {

inline std::string gitSHA1()
{
   return "@GIT_COMMIT_HASH@";
}

inline std::string gitBranch()
{
   return "@GIT_BRANCH@";
}

inline std::string gitDiff()
{
   return R"(@GIT_DIFF@)";
}

inline void printGitInfo()
{
   WALBERLA_LOG_INFO_ON_ROOT( "Git info:" )
   WALBERLA_LOG_INFO_ON_ROOT( " - SHA1:   " << gitSHA1() );
   WALBERLA_LOG_INFO_ON_ROOT( " - branch: " << gitBranch() );
   WALBERLA_LOG_INFO_ON_ROOT( " - diff:  \n" <<  gitDiff() );
}

}
