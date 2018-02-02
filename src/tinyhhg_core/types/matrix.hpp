#pragma once

#include "core/DataTypes.h"
#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"
#include "core/debug/Debug.h"

#include <cmath>
#include <array>
#include <iomanip>

namespace hhg {

using walberla::real_t;
using walberla::uint_t;

/// \brief  NxM Matrix
/// \author Daniel Drzisga (drzisga@ma.tum.de)
/// \date   September, 2017
///
/// The Matrix class represents an MxN-dimensional matrix with basic support for algebraic operations
/// \tparam T Matrix value data type
/// \tparam M Number of rows
/// \tparam N Number of columns
template<typename T, uint_t M, uint_t N>
class Matrix {
public:
  static const uint_t Size = M * N;

  /// Default constructor setting all components to zero
  Matrix() {
    for (uint_t i = 0; i < Size; ++i) {
      x[i] = (T) 0;
    }
  }

  /// Constructs the matrix using values from M*N-dimensional array \p _x in row-major order
  /// \param _x Pointer to M*N-dimensional array
  Matrix(T _x[Size])
  {
    for (uint_t i = 0; i < Size; ++i)
    {
      x[i] = _x[i];
    }
  }

  /// Copy constructor
  /// \param b Reference to another instance of Matrix
  Matrix(const Matrix& b)
  {
    for (uint_t i = 0; i < Size; ++i)
    {
      x[i] = b.x[i];
    }
  }

  /// Get reference to a single matrix component
  /// \param row Row index
  /// \param col Column index
  /// \returns Reference to component at position [row,col] in matrix
  T& operator() (uint_t row, uint_t col)
  {
    WALBERLA_ASSERT(row < M, "Matrix row index out of bounds: row = " << row << " but M = " << M);
    WALBERLA_ASSERT(col < N, "Matrix column index out of bounds: col = " << col << " but N = " << N);
    return x[N*row + col];
  }

  /// Get const reference to a single matrix component
  /// \param row Row index
  /// \param col Column index
  /// \returns Const reference to component at position [row,col] in matrix
  const T& operator() (uint_t row, uint_t col) const
  {
    WALBERLA_ASSERT(row < M, "Matrix row index out of bounds: row = " << row << " but M = " << M);
    WALBERLA_ASSERT(col < N, "Matrix column index out of bounds: col = " << col << " but N = " << N);
    return x[N*row + col];
  }

  /// Get raw pointer to underlying matrix data
  /// \returns Pointer to first element of unerlying matrix data
  T* data()
  {
    return &x[0];
  }

  /// Get const raw pointer to underlying matrix data
  /// \returns Constant pointer to first element of unerlying matrix data
  const T* data() const
  {
    return &x[0];
  }

private:
  T x[Size];
};

template<typename T, uint_t M, uint_t N>
inline std::ostream& operator<<(std::ostream &os, const Matrix<T, M, N> &matrix)
{
  os << "[\n";

  for (uint_t i = 0; i < M; ++i) {
    os << "[";
    for (uint_t j = 0; j < N; ++j) {
      os << std::scientific << std::setw(13) << matrix(i,j);
      if (j != N - 1) {
        os << ", ";
      }
    }
    os << "],\n";
  }

  os << "]";

  return os;
}

typedef Matrix<real_t, 3, 3> Matrix3r;
typedef Matrix<real_t, 6, 6> Matrix6r;

}
