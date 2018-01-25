#pragma once

#pragma once

#include <string>
#include <cstdio>

namespace hhg {

/// format uses the printf syntax to format a given formatString and return an std::string
///\param formatString the format string
///\param all arguments which will be inserted into the string, these CANNOT be std::string but have to be converted using .c_str()
template<typename... Args>
std::string format(const char *formatString, Args &&... args) {
  ///this size is arbitrary
  size_t maxBufferSize = 4096;
  char buffer[maxBufferSize];
  int check = snprintf(buffer, maxBufferSize, formatString, args...);
  if (check <= 0 || check > maxBufferSize) {
    WALBERLA_ABORT("snprintf failed");
  }
  return std::string(buffer);
}

}/// namespace tinyhhg