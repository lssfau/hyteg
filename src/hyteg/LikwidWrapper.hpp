///
/// \brief Wrapper for the likwid tools
///
/// \file likwidwrapper.hpp
///
/// This file is used to include likwid if desired or take care of the macros if not
///

#ifndef PROJECT_LIKWIDWRAPPER_HPP_H
#define PROJECT_LIKWIDWRAPPER_HPP_H

#ifdef LIKWID_PERFMON
#include <likwid.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#define LIKWID_MARKER_RESET(regionTag)
#endif // HHG_LIKWID


#endif //PROJECT_LIKWIDWRAPPER_HPP_H
