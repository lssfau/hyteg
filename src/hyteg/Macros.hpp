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

#define SPECIALIZE_OPRTYPE_POLYNOMIAL(type, func_tmpl, func)          \
template< typename ValueType, OperatorType OprType, typename... Args> \
inline type func(uint_t polyDegree, Args&&... args)                   \
{                                                                     \
  switch(polyDegree) {                                                \
    case 0:                                                           \
      return func_tmpl< ValueType, OprType, 0 >(args...);             \
    case 1:                                                           \
      return func_tmpl< ValueType, OprType, 1 >(args...);             \
    case 2:                                                           \
      return func_tmpl< ValueType, OprType, 2 >(args...);             \
    case 3:                                                           \
      return func_tmpl< ValueType, OprType, 3 >(args...);             \
    case 4:                                                           \
      return func_tmpl< ValueType, OprType, 4 >(args...);             \
    case 5:                                                           \
      return func_tmpl< ValueType, OprType, 5 >(args...);             \
    case 6:                                                           \
      return func_tmpl< ValueType, OprType, 6 >(args...);             \
    case 7:                                                           \
      return func_tmpl< ValueType, OprType, 7 >(args...);             \
    case 8:                                                           \
      return func_tmpl< ValueType, OprType, 8 >(args...);             \
    case 9:                                                           \
      return func_tmpl< ValueType, OprType, 9 >(args...);             \
    case 10:                                                          \
      return func_tmpl< ValueType, OprType, 10 >(args...);            \
    case 11:                                                          \
      return func_tmpl< ValueType, OprType, 11 >(args...);            \
    case 12:                                                          \
      return func_tmpl< ValueType, OprType, 12 >(args...);            \
    default:                                                          \
      WALBERLA_ABORT("Degree " << polyDegree << " not supported")     \
  }                                                                   \
}

#define SPECIALIZE_R_POLYNOMIAL(type, func_tmpl, func)              \
template <typename... Args>                                         \
inline type func(uint_t polyDegree, Args&&... args)                 \
{                                                                   \
  switch(polyDegree) {                                              \
    case 0:                                                         \
      return func_tmpl< 0 >(args...);                               \
    case 1:                                                         \
      return func_tmpl< 1 >(args...);                               \
    case 2:                                                         \
      return func_tmpl< 2 >(args...);                               \
    case 3:                                                         \
      return func_tmpl< 3 >(args...);                               \
    case 4:                                                         \
      return func_tmpl< 4 >(args...);                               \
    case 5:                                                         \
      return func_tmpl< 5 >(args...);                               \
    case 6:                                                         \
      return func_tmpl< 6 >(args...);                               \
    case 7:                                                         \
      return func_tmpl< 7 >(args...);                               \
    case 8:                                                         \
      return func_tmpl< 8 >(args...);                               \
    case 9:                                                         \
      return func_tmpl< 9 >(args...);                               \
    case 10:                                                        \
      return func_tmpl< 10 >(args...);                              \
    case 11:                                                        \
      return func_tmpl< 11 >(args...);                              \
    case 12:                                                        \
      return func_tmpl< 12 >(args...);                              \
    default:                                                        \
      WALBERLA_ABORT("Degree " << polyDegree << " not supported")   \
  }                                                                 \
}
