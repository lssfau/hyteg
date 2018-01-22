#ifndef TINYHHG_MACROS_HPP
#define TINYHHG_MACROS_HPP

#define SPECIALIZE_WITH_VALUETYPE(type, func_tmpl, func)                     \
template< typename ValueType, typename... Args>                                    \
inline type func(uint_t level, Args&&... args)                  \
{                                                             \
  switch(level)                                               \
  {                                                           \
    case 2:                                                   \
      return func_tmpl< ValueType,  2 >(args...);                           \
    case 3:                                                   \
      return func_tmpl< ValueType,  3 >(args...);                           \
    case 4:                                                   \
      return func_tmpl< ValueType,  4 >(args...);                           \
    case 5:                                                   \
      return func_tmpl< ValueType,  5 >(args...);                           \
    case 6:                                                   \
      return func_tmpl< ValueType,  6 >(args...);                           \
    case 7:                                                   \
      return func_tmpl< ValueType,  7 >(args...);                           \
    case 8:                                                   \
      return func_tmpl< ValueType,  8 >(args...);                           \
    case 9:                                                   \
      return func_tmpl< ValueType,  9 >(args...);                           \
    case 10:                                                  \
      return func_tmpl< ValueType, 10 >(args...);                          \
    case 11:                                                  \
      return func_tmpl< ValueType, 11 >(args...);                          \
    case 12:                                                  \
      return func_tmpl< ValueType, 12 >(args...);                          \
    case 13:                                                  \
      return func_tmpl< ValueType, 13 >(args...);                          \
    case 14:                                                  \
      return func_tmpl< ValueType, 14 >(args...);                          \
    default:                                                  \
      WALBERLA_ABORT("Level " << level << " not supported")   \
  }                                                           \
}

#define SPECIALIZE(type, func_tmpl, func)                     \
template<typename... Args>                                    \
inline type func(uint_t level, Args&&... args)                  \
{                                                             \
  switch(level)                                               \
  {                                                           \
    case 2:                                                   \
      return func_tmpl<2>(args...);                           \
    case 3:                                                   \
      return func_tmpl<3>(args...);                           \
    case 4:                                                   \
      return func_tmpl<4>(args...);                           \
    case 5:                                                   \
      return func_tmpl<5>(args...);                           \
    case 6:                                                   \
      return func_tmpl<6>(args...);                           \
    case 7:                                                   \
      return func_tmpl<7>(args...);                           \
    case 8:                                                   \
      return func_tmpl<8>(args...);                           \
    case 9:                                                   \
      return func_tmpl<9>(args...);                           \
    case 10:                                                  \
      return func_tmpl<10>(args...);                          \
    case 11:                                                  \
      return func_tmpl<11>(args...);                          \
    case 12:                                                  \
      return func_tmpl<12>(args...);                          \
    case 13:                                                  \
      return func_tmpl<13>(args...);                          \
    case 14:                                                  \
      return func_tmpl<14>(args...);                          \
    default:                                                  \
      WALBERLA_ABORT("Level " << level << " not supported")   \
  }                                                           \
}

#define SPECIALIZE_POLYNOMIAL(type, func_tmpl, func)                     \
template< typename ValueType, uint_t MaxPolyDegree, typename... Args>                                    \
inline type func(uint_t level, Args&&... args)                  \
{                                                             \
  switch(level)                                               \
  {                                                           \
    case 2:                                                   \
      return func_tmpl< ValueType, MaxPolyDegree,  2 >(args...);                           \
    case 3:                                                   \
      return func_tmpl< ValueType, MaxPolyDegree,  3 >(args...);                           \
    case 4:                                                   \
      return func_tmpl< ValueType, MaxPolyDegree,  4 >(args...);                           \
    case 5:                                                   \
      return func_tmpl< ValueType, MaxPolyDegree,  5 >(args...);                           \
    case 6:                                                   \
      return func_tmpl< ValueType, MaxPolyDegree,  6 >(args...);                           \
    case 7:                                                   \
      return func_tmpl< ValueType, MaxPolyDegree,  7 >(args...);                           \
    case 8:                                                   \
      return func_tmpl< ValueType, MaxPolyDegree,  8 >(args...);                           \
    case 9:                                                   \
      return func_tmpl< ValueType, MaxPolyDegree,  9 >(args...);                           \
    case 10:                                                  \
      return func_tmpl< ValueType, MaxPolyDegree, 10 >(args...);                          \
    case 11:                                                  \
      return func_tmpl< ValueType, MaxPolyDegree, 11 >(args...);                          \
    case 12:                                                  \
      return func_tmpl< ValueType, MaxPolyDegree, 12 >(args...);                          \
    case 13:                                                  \
      return func_tmpl< ValueType, MaxPolyDegree, 13 >(args...);                          \
    case 14:                                                  \
      return func_tmpl< ValueType, MaxPolyDegree, 14 >(args...);                          \
    default:                                                  \
      WALBERLA_ABORT("Level " << level << " not supported")   \
  }                                                           \
}

#endif //TINYHHG_MACROS_HPP
