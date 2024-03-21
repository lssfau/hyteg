/*
* Copyright (c) 2017-2023 Nils Kohl, Daniel Bauer, Fabian Böhm.
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

/*
* The entire file was generated with the HyTeG form generator.
*
* Avoid modifying this file. If buggy, consider fixing the generator itself.
*/

#pragma once

#include "hyteg/geometry/GeometryMap.hpp"
#include "hyteg/forms/form_hyteg_base/P2FormHyTeG.hpp"

namespace hyteg {
namespace forms {

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p2_epsilonvar_0_0_affine_q4
/// - description:
/// - trial space: Lagrange, degree: 2
/// - test space:  Lagrange, degree: 2
///
class p2_epsilonvar_0_0_affine_q4 : public P2FormHyTeG
{

 public:

   p2_epsilonvar_0_0_affine_q4() { WALBERLA_ABORT("Not implemented."); }

   p2_epsilonvar_0_0_affine_q4( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_2D_mu, std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
   : callback_Scalar_Variable_Coefficient_2D_mu(_callback_Scalar_Variable_Coefficient_2D_mu)
   , callback_Scalar_Variable_Coefficient_3D_mu(_callback_Scalar_Variable_Coefficient_3D_mu)
   {}

 private:

   std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_mu;
   std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_mu;


 public:

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3
   /// - element matrix dimensions (rows, cols): (6, 6)
   /// - quadrature rule:                        Vioreanu-Rokhlin 2 | points: 6, degree: 4, test tolerance: 1.943e-16
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              393     837       2       0      1            418                 6
   ///
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 6, 6 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3
   /// - element matrix dimensions (rows, cols): (6, 6)
   /// - quadrature rule:                        Vioreanu-Rokhlin 2 | points: 6, degree: 4, test tolerance: 1.943e-16
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              192     336       1       0      1            170                 6
   ///
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 6 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (10, 10)
   /// - quadrature rule:                        Xiao-Gimbutas 4 | points: 11, degree: 4, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                             2752    4371       2       0      1           1819                11
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 10, 10 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (10, 10)
   /// - quadrature rule:                        Xiao-Gimbutas 4 | points: 11, degree: 4, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                             1015    1425       1       0      1            704                11
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 10 >& elMat ) const override;

 private:

   void Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t * out_0 ) const;

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p2_epsilonvar_0_1_affine_q4
/// - description:
/// - trial space: Lagrange, degree: 2
/// - test space:  Lagrange, degree: 2
///
class p2_epsilonvar_0_1_affine_q4 : public P2FormHyTeG
{

 public:

   p2_epsilonvar_0_1_affine_q4() { WALBERLA_ABORT("Not implemented."); }

   p2_epsilonvar_0_1_affine_q4( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_2D_mu, std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
   : callback_Scalar_Variable_Coefficient_2D_mu(_callback_Scalar_Variable_Coefficient_2D_mu)
   , callback_Scalar_Variable_Coefficient_3D_mu(_callback_Scalar_Variable_Coefficient_3D_mu)
   {}

 private:

   std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_mu;
   std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_mu;


 public:

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3
   /// - element matrix dimensions (rows, cols): (6, 6)
   /// - quadrature rule:                        Vioreanu-Rokhlin 2 | points: 6, degree: 4, test tolerance: 1.943e-16
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              300     503       2       0      1            341                 6
   ///
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 6, 6 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3
   /// - element matrix dimensions (rows, cols): (6, 6)
   /// - quadrature rule:                        Vioreanu-Rokhlin 2 | points: 6, degree: 4, test tolerance: 1.943e-16
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              132     173       1       0      1             83                 6
   ///
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 6 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (10, 10)
   /// - quadrature rule:                        Xiao-Gimbutas 4 | points: 11, degree: 4, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                             1527    1985       2       0      1           1070                11
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 10, 10 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (10, 10)
   /// - quadrature rule:                        Xiao-Gimbutas 4 | points: 11, degree: 4, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              528     550       1       0      1            258                11
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 10 >& elMat ) const override;

 private:

   void Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t * out_0 ) const;

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p2_epsilonvar_0_2_affine_q4
/// - description:
/// - trial space: Lagrange, degree: 2
/// - test space:  Lagrange, degree: 2
///
class p2_epsilonvar_0_2_affine_q4 : public P2FormHyTeG
{

 public:

   p2_epsilonvar_0_2_affine_q4() { WALBERLA_ABORT("Not implemented."); }

   p2_epsilonvar_0_2_affine_q4( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
   : callback_Scalar_Variable_Coefficient_3D_mu(_callback_Scalar_Variable_Coefficient_3D_mu)
   {}

 private:

   std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_mu;


 public:

   /// \brief Not implemented - does nothing.
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 6, 6 >& elMat ) const override;

   /// \brief Not implemented - does nothing.
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 6 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (10, 10)
   /// - quadrature rule:                        Xiao-Gimbutas 4 | points: 11, degree: 4, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                             1527    1985       2       0      1           1070                11
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 10, 10 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (10, 10)
   /// - quadrature rule:                        Xiao-Gimbutas 4 | points: 11, degree: 4, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              528     550       1       0      1            258                11
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 10 >& elMat ) const override;

 private:

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p2_epsilonvar_1_0_affine_q4
/// - description:
/// - trial space: Lagrange, degree: 2
/// - test space:  Lagrange, degree: 2
///
class p2_epsilonvar_1_0_affine_q4 : public P2FormHyTeG
{

 public:

   p2_epsilonvar_1_0_affine_q4() { WALBERLA_ABORT("Not implemented."); }

   p2_epsilonvar_1_0_affine_q4( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_2D_mu, std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
   : callback_Scalar_Variable_Coefficient_2D_mu(_callback_Scalar_Variable_Coefficient_2D_mu)
   , callback_Scalar_Variable_Coefficient_3D_mu(_callback_Scalar_Variable_Coefficient_3D_mu)
   {}

 private:

   std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_mu;
   std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_mu;


 public:

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3
   /// - element matrix dimensions (rows, cols): (6, 6)
   /// - quadrature rule:                        Vioreanu-Rokhlin 2 | points: 6, degree: 4, test tolerance: 1.943e-16
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              300     503       2       0      1            341                 6
   ///
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 6, 6 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3
   /// - element matrix dimensions (rows, cols): (6, 6)
   /// - quadrature rule:                        Vioreanu-Rokhlin 2 | points: 6, degree: 4, test tolerance: 1.943e-16
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              132     173       1       0      1             83                 6
   ///
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 6 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (10, 10)
   /// - quadrature rule:                        Xiao-Gimbutas 4 | points: 11, degree: 4, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                             1527    1985       2       0      1           1070                11
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 10, 10 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (10, 10)
   /// - quadrature rule:                        Xiao-Gimbutas 4 | points: 11, degree: 4, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              528     550       1       0      1            258                11
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 10 >& elMat ) const override;

 private:

   void Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t * out_0 ) const;

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p2_epsilonvar_1_1_affine_q4
/// - description:
/// - trial space: Lagrange, degree: 2
/// - test space:  Lagrange, degree: 2
///
class p2_epsilonvar_1_1_affine_q4 : public P2FormHyTeG
{

 public:

   p2_epsilonvar_1_1_affine_q4() { WALBERLA_ABORT("Not implemented."); }

   p2_epsilonvar_1_1_affine_q4( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_2D_mu, std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
   : callback_Scalar_Variable_Coefficient_2D_mu(_callback_Scalar_Variable_Coefficient_2D_mu)
   , callback_Scalar_Variable_Coefficient_3D_mu(_callback_Scalar_Variable_Coefficient_3D_mu)
   {}

 private:

   std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_mu;
   std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_mu;


 public:

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3
   /// - element matrix dimensions (rows, cols): (6, 6)
   /// - quadrature rule:                        Vioreanu-Rokhlin 2 | points: 6, degree: 4, test tolerance: 1.943e-16
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              393     835       2       0      1            418                 6
   ///
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 6, 6 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3
   /// - element matrix dimensions (rows, cols): (6, 6)
   /// - quadrature rule:                        Vioreanu-Rokhlin 2 | points: 6, degree: 4, test tolerance: 1.943e-16
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              192     336       1       0      1            170                 6
   ///
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 6 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (10, 10)
   /// - quadrature rule:                        Xiao-Gimbutas 4 | points: 11, degree: 4, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                             2752    4371       2       0      1           1819                11
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 10, 10 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (10, 10)
   /// - quadrature rule:                        Xiao-Gimbutas 4 | points: 11, degree: 4, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                             1015    1434       1       0      1            713                11
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 10 >& elMat ) const override;

 private:

   void Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t * out_0 ) const;

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p2_epsilonvar_1_2_affine_q4
/// - description:
/// - trial space: Lagrange, degree: 2
/// - test space:  Lagrange, degree: 2
///
class p2_epsilonvar_1_2_affine_q4 : public P2FormHyTeG
{

 public:

   p2_epsilonvar_1_2_affine_q4() { WALBERLA_ABORT("Not implemented."); }

   p2_epsilonvar_1_2_affine_q4( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
   : callback_Scalar_Variable_Coefficient_3D_mu(_callback_Scalar_Variable_Coefficient_3D_mu)
   {}

 private:

   std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_mu;


 public:

   /// \brief Not implemented - does nothing.
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 6, 6 >& elMat ) const override;

   /// \brief Not implemented - does nothing.
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 6 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (10, 10)
   /// - quadrature rule:                        Xiao-Gimbutas 4 | points: 11, degree: 4, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                             1527    1985       2       0      1           1070                11
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 10, 10 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (10, 10)
   /// - quadrature rule:                        Xiao-Gimbutas 4 | points: 11, degree: 4, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              528     550       1       0      1            258                11
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 10 >& elMat ) const override;

 private:

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p2_epsilonvar_2_0_affine_q4
/// - description:
/// - trial space: Lagrange, degree: 2
/// - test space:  Lagrange, degree: 2
///
class p2_epsilonvar_2_0_affine_q4 : public P2FormHyTeG
{

 public:

   p2_epsilonvar_2_0_affine_q4() { WALBERLA_ABORT("Not implemented."); }

   p2_epsilonvar_2_0_affine_q4( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
   : callback_Scalar_Variable_Coefficient_3D_mu(_callback_Scalar_Variable_Coefficient_3D_mu)
   {}

 private:

   std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_mu;


 public:

   /// \brief Not implemented - does nothing.
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 6, 6 >& elMat ) const override;

   /// \brief Not implemented - does nothing.
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 6 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (10, 10)
   /// - quadrature rule:                        Xiao-Gimbutas 4 | points: 11, degree: 4, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                             1527    1985       2       0      1           1070                11
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 10, 10 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (10, 10)
   /// - quadrature rule:                        Xiao-Gimbutas 4 | points: 11, degree: 4, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              528     550       1       0      1            258                11
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 10 >& elMat ) const override;

 private:

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p2_epsilonvar_2_1_affine_q4
/// - description:
/// - trial space: Lagrange, degree: 2
/// - test space:  Lagrange, degree: 2
///
class p2_epsilonvar_2_1_affine_q4 : public P2FormHyTeG
{

 public:

   p2_epsilonvar_2_1_affine_q4() { WALBERLA_ABORT("Not implemented."); }

   p2_epsilonvar_2_1_affine_q4( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
   : callback_Scalar_Variable_Coefficient_3D_mu(_callback_Scalar_Variable_Coefficient_3D_mu)
   {}

 private:

   std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_mu;


 public:

   /// \brief Not implemented - does nothing.
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 6, 6 >& elMat ) const override;

   /// \brief Not implemented - does nothing.
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 6 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (10, 10)
   /// - quadrature rule:                        Xiao-Gimbutas 4 | points: 11, degree: 4, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                             1527    1985       2       0      1           1070                11
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 10, 10 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (10, 10)
   /// - quadrature rule:                        Xiao-Gimbutas 4 | points: 11, degree: 4, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              528     550       1       0      1            258                11
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 10 >& elMat ) const override;

 private:

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p2_epsilonvar_2_2_affine_q4
/// - description:
/// - trial space: Lagrange, degree: 2
/// - test space:  Lagrange, degree: 2
///
class p2_epsilonvar_2_2_affine_q4 : public P2FormHyTeG
{

 public:

   p2_epsilonvar_2_2_affine_q4() { WALBERLA_ABORT("Not implemented."); }

   p2_epsilonvar_2_2_affine_q4( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
   : callback_Scalar_Variable_Coefficient_3D_mu(_callback_Scalar_Variable_Coefficient_3D_mu)
   {}

 private:

   std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_mu;


 public:

   /// \brief Not implemented - does nothing.
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 6, 6 >& elMat ) const override;

   /// \brief Not implemented - does nothing.
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 6 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (10, 10)
   /// - quadrature rule:                        Xiao-Gimbutas 4 | points: 11, degree: 4, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                             2752    4371       2       0      1           1819                11
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 10, 10 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (10, 10)
   /// - quadrature rule:                        Xiao-Gimbutas 4 | points: 11, degree: 4, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                             1015    1417       1       0      1            699                11
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 10 >& elMat ) const override;

 private:

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

} // namespace forms
} // namespace hyteg