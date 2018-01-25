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

#define SPECIALIZE_POLYNOMIAL(type, func_tmpl, func)                      \
template< typename ValueType, typename... Args>                           \
inline type func(uint_t level, uint_t polyDegree, Args&&... args)         \
{                                                                         \
  switch(level)                                                           \
  {                                                                       \
    case 2:                                                               \
      {                                                                   \
        switch(polyDegree) {                                              \
          case 0:                                                         \
            return func_tmpl< ValueType, 2, 0 >(args...);                 \
          case 1:                                                         \
            return func_tmpl< ValueType, 2, 1 >(args...);                 \
          case 2:                                                         \
            return func_tmpl< ValueType, 2, 2 >(args...);                 \
          case 3:                                                         \
            return func_tmpl< ValueType, 2, 3 >(args...);                 \
          case 4:                                                         \
            return func_tmpl< ValueType, 2, 4 >(args...);                 \
          case 5:                                                         \
            return func_tmpl< ValueType, 2, 5 >(args...);                 \
          case 6:                                                         \
            return func_tmpl< ValueType, 2, 6 >(args...);                 \
          case 7:                                                         \
            return func_tmpl< ValueType, 2, 7 >(args...);                 \
          default:                                                        \
            WALBERLA_ABORT("Degree " << polyDegree << " not supported")   \
        }                                                                 \
      }                                                                   \
    case 3:                                                               \
      {                                                                   \
        switch(polyDegree) {                                              \
          case 0:                                                         \
            return func_tmpl< ValueType, 3, 0 >(args...);                 \
          case 1:                                                         \
            return func_tmpl< ValueType, 3, 1 >(args...);                 \
          case 2:                                                         \
            return func_tmpl< ValueType, 3, 2 >(args...);                 \
          case 3:                                                         \
            return func_tmpl< ValueType, 3, 3 >(args...);                 \
          case 4:                                                         \
            return func_tmpl< ValueType, 3, 4 >(args...);                 \
          case 5:                                                         \
            return func_tmpl< ValueType, 3, 5 >(args...);                 \
          case 6:                                                         \
            return func_tmpl< ValueType, 3, 6 >(args...);                 \
          case 7:                                                         \
            return func_tmpl< ValueType, 3, 7 >(args...);                 \
          default:                                                        \
            WALBERLA_ABORT("Degree " << polyDegree << " not supported")   \
        }                                                                 \
      }                                                                   \
    case 4:                                                               \
      {                                                                   \
        switch(polyDegree) {                                              \
          case 0:                                                         \
            return func_tmpl< ValueType, 4, 0 >(args...);                 \
          case 1:                                                         \
            return func_tmpl< ValueType, 4, 1 >(args...);                 \
          case 2:                                                         \
            return func_tmpl< ValueType, 4, 2 >(args...);                 \
          case 3:                                                         \
            return func_tmpl< ValueType, 4, 3 >(args...);                 \
          case 4:                                                         \
            return func_tmpl< ValueType, 4, 4 >(args...);                 \
          case 5:                                                         \
            return func_tmpl< ValueType, 4, 5 >(args...);                 \
          case 6:                                                         \
            return func_tmpl< ValueType, 4, 6 >(args...);                 \
          case 7:                                                         \
            return func_tmpl< ValueType, 4, 7 >(args...);                 \
          default:                                                        \
            WALBERLA_ABORT("Degree " << polyDegree << " not supported")   \
        }                                                                 \
      }                                                                   \
    case 5:                                                               \
      {                                                                   \
        switch(polyDegree) {                                              \
          case 0:                                                         \
            return func_tmpl< ValueType, 5, 0 >(args...);                 \
          case 1:                                                         \
            return func_tmpl< ValueType, 5, 1 >(args...);                 \
          case 2:                                                         \
            return func_tmpl< ValueType, 5, 2 >(args...);                 \
          case 3:                                                         \
            return func_tmpl< ValueType, 5, 3 >(args...);                 \
          case 4:                                                         \
            return func_tmpl< ValueType, 5, 4 >(args...);                 \
          case 5:                                                         \
            return func_tmpl< ValueType, 5, 5 >(args...);                 \
          case 6:                                                         \
            return func_tmpl< ValueType, 5, 6 >(args...);                 \
          case 7:                                                         \
            return func_tmpl< ValueType, 5, 7 >(args...);                 \
          default:                                                        \
            WALBERLA_ABORT("Degree " << polyDegree << " not supported")   \
        }                                                                 \
      }                                                                   \
    case 6:                                                               \
      {                                                                   \
        switch(polyDegree) {                                              \
          case 0:                                                         \
            return func_tmpl< ValueType, 6, 0 >(args...);                 \
          case 1:                                                         \
            return func_tmpl< ValueType, 6, 1 >(args...);                 \
          case 2:                                                         \
            return func_tmpl< ValueType, 6, 2 >(args...);                 \
          case 3:                                                         \
            return func_tmpl< ValueType, 6, 3 >(args...);                 \
          case 4:                                                         \
            return func_tmpl< ValueType, 6, 4 >(args...);                 \
          case 5:                                                         \
            return func_tmpl< ValueType, 6, 5 >(args...);                 \
          case 6:                                                         \
            return func_tmpl< ValueType, 6, 6 >(args...);                 \
          case 7:                                                         \
            return func_tmpl< ValueType, 6, 7 >(args...);                 \
          default:                                                        \
            WALBERLA_ABORT("Degree " << polyDegree << " not supported")   \
        }                                                                 \
      }                                                                   \
    case 7:                                                               \
      {                                                                   \
        switch(polyDegree) {                                              \
          case 0:                                                         \
            return func_tmpl< ValueType, 7, 0 >(args...);                 \
          case 1:                                                         \
            return func_tmpl< ValueType, 7, 1 >(args...);                 \
          case 2:                                                         \
            return func_tmpl< ValueType, 7, 2 >(args...);                 \
          case 3:                                                         \
            return func_tmpl< ValueType, 7, 3 >(args...);                 \
          case 4:                                                         \
            return func_tmpl< ValueType, 7, 4 >(args...);                 \
          case 5:                                                         \
            return func_tmpl< ValueType, 7, 5 >(args...);                 \
          case 6:                                                         \
            return func_tmpl< ValueType, 7, 6 >(args...);                 \
          case 7:                                                         \
            return func_tmpl< ValueType, 7, 7 >(args...);                 \
          default:                                                        \
            WALBERLA_ABORT("Degree " << polyDegree << " not supported")   \
        }                                                                 \
      }                                                                   \
    case 8:                                                               \
      {                                                                   \
        switch(polyDegree) {                                              \
          case 0:                                                         \
            return func_tmpl< ValueType, 8, 0 >(args...);                 \
          case 1:                                                         \
            return func_tmpl< ValueType, 8, 1 >(args...);                 \
          case 2:                                                         \
            return func_tmpl< ValueType, 8, 2 >(args...);                 \
          case 3:                                                         \
            return func_tmpl< ValueType, 8, 3 >(args...);                 \
          case 4:                                                         \
            return func_tmpl< ValueType, 8, 4 >(args...);                 \
          case 5:                                                         \
            return func_tmpl< ValueType, 8, 5 >(args...);                 \
          case 6:                                                         \
            return func_tmpl< ValueType, 8, 6 >(args...);                 \
          case 7:                                                         \
            return func_tmpl< ValueType, 8, 7 >(args...);                 \
          default:                                                        \
            WALBERLA_ABORT("Degree " << polyDegree << " not supported")   \
        }                                                                 \
      }                                                                   \
    case 9:                                                               \
      {                                                                   \
        switch(polyDegree) {                                              \
          case 0:                                                         \
            return func_tmpl< ValueType, 9, 0 >(args...);                 \
          case 1:                                                         \
            return func_tmpl< ValueType, 9, 1 >(args...);                 \
          case 2:                                                         \
            return func_tmpl< ValueType, 9, 2 >(args...);                 \
          case 3:                                                         \
            return func_tmpl< ValueType, 9, 3 >(args...);                 \
          case 4:                                                         \
            return func_tmpl< ValueType, 9, 4 >(args...);                 \
          case 5:                                                         \
            return func_tmpl< ValueType, 9, 5 >(args...);                 \
          case 6:                                                         \
            return func_tmpl< ValueType, 9, 6 >(args...);                 \
          case 7:                                                         \
            return func_tmpl< ValueType, 9, 7 >(args...);                 \
          default:                                                        \
            WALBERLA_ABORT("Degree " << polyDegree << " not supported")   \
        }                                                                 \
      }                                                                   \
    case 10:                                                              \
      {                                                                   \
        switch(polyDegree) {                                              \
          case 0:                                                         \
            return func_tmpl< ValueType,10, 0 >(args...);                 \
          case 1:                                                         \
            return func_tmpl< ValueType,10, 1 >(args...);                 \
          case 2:                                                         \
            return func_tmpl< ValueType,10, 2 >(args...);                 \
          case 3:                                                         \
            return func_tmpl< ValueType,10, 3 >(args...);                 \
          case 4:                                                         \
            return func_tmpl< ValueType,10, 4 >(args...);                 \
          case 5:                                                         \
            return func_tmpl< ValueType,10, 5 >(args...);                 \
          case 6:                                                         \
            return func_tmpl< ValueType,10, 6 >(args...);                 \
          case 7:                                                         \
            return func_tmpl< ValueType,10, 7 >(args...);                 \
          default:                                                        \
            WALBERLA_ABORT("Degree " << polyDegree << " not supported")   \
        }                                                                 \
      }                                                                   \
    case 11:                                                              \
      {                                                                   \
        switch(polyDegree) {                                              \
          case 0:                                                         \
            return func_tmpl< ValueType,11, 0 >(args...);                 \
          case 1:                                                         \
            return func_tmpl< ValueType,11, 1 >(args...);                 \
          case 2:                                                         \
            return func_tmpl< ValueType,11, 2 >(args...);                 \
          case 3:                                                         \
            return func_tmpl< ValueType,11, 3 >(args...);                 \
          case 4:                                                         \
            return func_tmpl< ValueType,11, 4 >(args...);                 \
          case 5:                                                         \
            return func_tmpl< ValueType,11, 5 >(args...);                 \
          case 6:                                                         \
            return func_tmpl< ValueType,11, 6 >(args...);                 \
          case 7:                                                         \
            return func_tmpl< ValueType,11, 7 >(args...);                 \
          default:                                                        \
            WALBERLA_ABORT("Degree " << polyDegree << " not supported")   \
        }                                                                 \
      }                                                                   \
    case 12:                                                              \
      {                                                                   \
        switch(polyDegree) {                                              \
          case 0:                                                         \
            return func_tmpl< ValueType,12, 0 >(args...);                 \
          case 1:                                                         \
            return func_tmpl< ValueType,12, 1 >(args...);                 \
          case 2:                                                         \
            return func_tmpl< ValueType,12, 2 >(args...);                 \
          case 3:                                                         \
            return func_tmpl< ValueType,12, 3 >(args...);                 \
          case 4:                                                         \
            return func_tmpl< ValueType,12, 4 >(args...);                 \
          case 5:                                                         \
            return func_tmpl< ValueType,12, 5 >(args...);                 \
          case 6:                                                         \
            return func_tmpl< ValueType,12, 6 >(args...);                 \
          case 7:                                                         \
            return func_tmpl< ValueType,12, 7 >(args...);                 \
          default:                                                        \
            WALBERLA_ABORT("Degree " << polyDegree << " not supported")   \
        }                                                                 \
      }                                                                   \
    case 13:                                                               \
      {                                                                   \
        switch(polyDegree) {                                              \
          case 0:                                                         \
            return func_tmpl< ValueType,13, 0 >(args...);                 \
          case 1:                                                         \
            return func_tmpl< ValueType,13, 1 >(args...);                 \
          case 2:                                                         \
            return func_tmpl< ValueType,13, 2 >(args...);                 \
          case 3:                                                         \
            return func_tmpl< ValueType,13, 3 >(args...);                 \
          case 4:                                                         \
            return func_tmpl< ValueType,13, 4 >(args...);                 \
          case 5:                                                         \
            return func_tmpl< ValueType,13, 5 >(args...);                 \
          case 6:                                                         \
            return func_tmpl< ValueType,13, 6 >(args...);                 \
          case 7:                                                         \
            return func_tmpl< ValueType,13, 7 >(args...);                 \
          default:                                                        \
            WALBERLA_ABORT("Degree " << polyDegree << " not supported")   \
        }                                                                 \
      }                                                                   \
    case 14:                                                              \
      {                                                                   \
        switch(polyDegree) {                                              \
          case 0:                                                         \
            return func_tmpl< ValueType,14, 0 >(args...);                 \
          case 1:                                                         \
            return func_tmpl< ValueType,14, 1 >(args...);                 \
          case 2:                                                         \
            return func_tmpl< ValueType,14, 2 >(args...);                 \
          case 3:                                                         \
            return func_tmpl< ValueType,14, 3 >(args...);                 \
          case 4:                                                         \
            return func_tmpl< ValueType,14, 4 >(args...);                 \
          case 5:                                                         \
            return func_tmpl< ValueType,14, 5 >(args...);                 \
          case 6:                                                         \
            return func_tmpl< ValueType,14, 6 >(args...);                 \
          case 7:                                                         \
            return func_tmpl< ValueType,14, 7 >(args...);                 \
          default:                                                        \
            WALBERLA_ABORT("Degree " << polyDegree << " not supported")   \
        }                                                                 \
      }                                                                   \
    default:                                                              \
      WALBERLA_ABORT("Level " << level << " not supported")               \
  }                                                                       \
}

#endif //TINYHHG_MACROS_HPP
