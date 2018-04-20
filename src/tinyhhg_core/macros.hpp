#pragma once

#define SPECIALIZE_POLYNOMIAL(type, func_tmpl, func)                \
template< typename ValueType, typename... Args>                     \
inline type func(uint_t polyDegree, Args&&... args)                 \
{                                                                   \
  switch(polyDegree) {                                              \
    case 0:                                                         \
      return func_tmpl< ValueType, 0 >(args...);                    \
    case 1:                                                         \
      return func_tmpl< ValueType, 1 >(args...);                    \
    case 2:                                                         \
      return func_tmpl< ValueType, 2 >(args...);                    \
    case 3:                                                         \
      return func_tmpl< ValueType, 3 >(args...);                    \
    case 4:                                                         \
      return func_tmpl< ValueType, 4 >(args...);                    \
    case 5:                                                         \
      return func_tmpl< ValueType, 5 >(args...);                    \
    case 6:                                                         \
      return func_tmpl< ValueType, 6 >(args...);                    \
    case 7:                                                         \
      return func_tmpl< ValueType, 7 >(args...);                    \
    case 8:                                                         \
      return func_tmpl< ValueType, 8 >(args...);                    \
    case 9:                                                         \
      return func_tmpl< ValueType, 9 >(args...);                    \
    case 10:                                                        \
      return func_tmpl< ValueType, 10 >(args...);                   \
    case 11:                                                        \
      return func_tmpl< ValueType, 11 >(args...);                   \
    case 12:                                                        \
      return func_tmpl< ValueType, 12 >(args...);                   \
    default:                                                        \
      WALBERLA_ABORT("Degree " << polyDegree << " not supported")   \
  }                                                                 \
}