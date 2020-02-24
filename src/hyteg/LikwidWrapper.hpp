/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
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
///
/// \brief Wrapper for the likwid tools
///
/// This file is used to include likwid if desired or take care of the macros if not
///

#ifndef PROJECT_LIKWIDWRAPPER_HPP_H
#define PROJECT_LIKWIDWRAPPER_HPP_H

#ifdef LIKWID_PERFMON
#include <likwid.h>
#else
#define LIKWID_MARKER_INIT (void)0
#define LIKWID_MARKER_THREADINIT (void)0
#define LIKWID_MARKER_SWITCH (void)0
#define LIKWID_MARKER_REGISTER(regionTag) (void)0
#define LIKWID_MARKER_START(regionTag) (void)0
#define LIKWID_MARKER_STOP(regionTag) (void)0
#define LIKWID_MARKER_CLOSE (void)0
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count) (void)0
#define LIKWID_MARKER_RESET(regionTag) (void)0
#endif // LIKWID_PERFMON


#endif //PROJECT_LIKWIDWRAPPER_HPP_H
