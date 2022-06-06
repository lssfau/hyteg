/*
 * Copyright (c) 2017-2022 Nils Kohl.
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
 * Software:
 *
 * - quadpy version: 0.16.5
 *
 * Avoid modifying this file. If buggy, consider fixing the generator itself.
 */

#pragma once

#include "hyteg/geometry/GeometryMap.hpp"
#include "hyteg/forms/form_hyteg_base/P1FormHyTeG.hpp"

namespace hyteg {
namespace forms {

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p1_full_stokesvar_0_0_affine_q2
/// - description: 
/// - trial space: Lagrange, degree: 1
/// - test space:  Lagrange, degree: 1
///
class p1_full_stokesvar_0_0_affine_q2 : public P1FormHyTeG
{

 public:

   p1_full_stokesvar_0_0_affine_q2() { WALBERLA_ABORT("Not implemented."); }

   p1_full_stokesvar_0_0_affine_q2( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_2D_mu, std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
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
   /// - element matrix dimensions (rows, cols): (3, 3)
   /// - quadrature rule:                        Dunavant 2 | points: 3, degree: 2, test tolerance: 2.22e-16
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               67     129       2       0      1             74                 3
   ///
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3
   /// - element matrix dimensions (rows, cols): (3, 3)
   /// - quadrature rule:                        Dunavant 2 | points: 3, degree: 2, test tolerance: 2.22e-16
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               52      89       1       0      1             46                 3
   ///
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Vioreanu-Rokhlin 1 | points: 4, degree: 2, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              220     359       2       0      1            167                 4
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Vioreanu-Rokhlin 1 | points: 4, degree: 2, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              154     234       1       0      1            102                 4
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const override;

   bool assemble2D() const override { return true; }

   bool assembly2DDefined() const override { return true; }

   bool assemble3D() const override { return true; }

   bool assembly3DDefined() const override { return true; }

 private:

   void Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t * out_0 ) const;

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p1_full_stokesvar_0_1_affine_q2
/// - description: 
/// - trial space: Lagrange, degree: 1
/// - test space:  Lagrange, degree: 1
///
class p1_full_stokesvar_0_1_affine_q2 : public P1FormHyTeG
{

 public:

   p1_full_stokesvar_0_1_affine_q2() { WALBERLA_ABORT("Not implemented."); }

   p1_full_stokesvar_0_1_affine_q2( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_2D_mu, std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
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
   /// - element matrix dimensions (rows, cols): (3, 3)
   /// - quadrature rule:                        Dunavant 2 | points: 3, degree: 2, test tolerance: 2.22e-16
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               67     140       2       0      1             80                 3
   ///
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3
   /// - element matrix dimensions (rows, cols): (3, 3)
   /// - quadrature rule:                        Dunavant 2 | points: 3, degree: 2, test tolerance: 2.22e-16
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               43      75       1       0      1             40                 3
   ///
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Vioreanu-Rokhlin 1 | points: 4, degree: 2, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              187     340       2       0      1            158                 4
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Vioreanu-Rokhlin 1 | points: 4, degree: 2, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              115     181       1       0      1             81                 4
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const override;

   bool assemble2D() const override { return true; }

   bool assembly2DDefined() const override { return true; }

   bool assemble3D() const override { return true; }

   bool assembly3DDefined() const override { return true; }

 private:

   void Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t * out_0 ) const;

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p1_full_stokesvar_0_2_affine_q2
/// - description: 
/// - trial space: Lagrange, degree: 1
/// - test space:  Lagrange, degree: 1
///
class p1_full_stokesvar_0_2_affine_q2 : public P1FormHyTeG
{

 public:

   p1_full_stokesvar_0_2_affine_q2() { WALBERLA_ABORT("Not implemented."); }

   p1_full_stokesvar_0_2_affine_q2( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
   : callback_Scalar_Variable_Coefficient_3D_mu(_callback_Scalar_Variable_Coefficient_3D_mu)
   {}

 private:

   std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_mu;


 public:

   /// \brief Not implemented - does nothing.
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const override;

   /// \brief Not implemented - does nothing.
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Vioreanu-Rokhlin 1 | points: 4, degree: 2, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              187     340       2       0      1            158                 4
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Vioreanu-Rokhlin 1 | points: 4, degree: 2, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              115     181       1       0      1             81                 4
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const override;

   bool assemble2D() const override { return false; }

   bool assembly2DDefined() const override { return false; }

   bool assemble3D() const override { return true; }

   bool assembly3DDefined() const override { return true; }

 private:

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p1_full_stokesvar_1_0_affine_q2
/// - description: 
/// - trial space: Lagrange, degree: 1
/// - test space:  Lagrange, degree: 1
///
class p1_full_stokesvar_1_0_affine_q2 : public P1FormHyTeG
{

 public:

   p1_full_stokesvar_1_0_affine_q2() { WALBERLA_ABORT("Not implemented."); }

   p1_full_stokesvar_1_0_affine_q2( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_2D_mu, std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
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
   /// - element matrix dimensions (rows, cols): (3, 3)
   /// - quadrature rule:                        Dunavant 2 | points: 3, degree: 2, test tolerance: 2.22e-16
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               67     140       2       0      1             80                 3
   ///
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3
   /// - element matrix dimensions (rows, cols): (3, 3)
   /// - quadrature rule:                        Dunavant 2 | points: 3, degree: 2, test tolerance: 2.22e-16
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               43      75       1       0      1             40                 3
   ///
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Vioreanu-Rokhlin 1 | points: 4, degree: 2, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              187     340       2       0      1            158                 4
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Vioreanu-Rokhlin 1 | points: 4, degree: 2, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              115     181       1       0      1             81                 4
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const override;

   bool assemble2D() const override { return true; }

   bool assembly2DDefined() const override { return true; }

   bool assemble3D() const override { return true; }

   bool assembly3DDefined() const override { return true; }

 private:

   void Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t * out_0 ) const;

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p1_full_stokesvar_1_1_affine_q2
/// - description: 
/// - trial space: Lagrange, degree: 1
/// - test space:  Lagrange, degree: 1
///
class p1_full_stokesvar_1_1_affine_q2 : public P1FormHyTeG
{

 public:

   p1_full_stokesvar_1_1_affine_q2() { WALBERLA_ABORT("Not implemented."); }

   p1_full_stokesvar_1_1_affine_q2( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_2D_mu, std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
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
   /// - element matrix dimensions (rows, cols): (3, 3)
   /// - quadrature rule:                        Dunavant 2 | points: 3, degree: 2, test tolerance: 2.22e-16
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               67     124       2       0      1             72                 3
   ///
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3
   /// - element matrix dimensions (rows, cols): (3, 3)
   /// - quadrature rule:                        Dunavant 2 | points: 3, degree: 2, test tolerance: 2.22e-16
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               52      84       1       0      1             44                 3
   ///
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Vioreanu-Rokhlin 1 | points: 4, degree: 2, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              220     359       2       0      1            167                 4
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Vioreanu-Rokhlin 1 | points: 4, degree: 2, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              154     234       1       0      1            102                 4
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const override;

   bool assemble2D() const override { return true; }

   bool assembly2DDefined() const override { return true; }

   bool assemble3D() const override { return true; }

   bool assembly3DDefined() const override { return true; }

 private:

   void Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t * out_0 ) const;

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p1_full_stokesvar_1_2_affine_q2
/// - description: 
/// - trial space: Lagrange, degree: 1
/// - test space:  Lagrange, degree: 1
///
class p1_full_stokesvar_1_2_affine_q2 : public P1FormHyTeG
{

 public:

   p1_full_stokesvar_1_2_affine_q2() { WALBERLA_ABORT("Not implemented."); }

   p1_full_stokesvar_1_2_affine_q2( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
   : callback_Scalar_Variable_Coefficient_3D_mu(_callback_Scalar_Variable_Coefficient_3D_mu)
   {}

 private:

   std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_mu;


 public:

   /// \brief Not implemented - does nothing.
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const override;

   /// \brief Not implemented - does nothing.
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Vioreanu-Rokhlin 1 | points: 4, degree: 2, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              187     340       2       0      1            158                 4
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Vioreanu-Rokhlin 1 | points: 4, degree: 2, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              115     181       1       0      1             81                 4
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const override;

   bool assemble2D() const override { return false; }

   bool assembly2DDefined() const override { return false; }

   bool assemble3D() const override { return true; }

   bool assembly3DDefined() const override { return true; }

 private:

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p1_full_stokesvar_2_0_affine_q2
/// - description: 
/// - trial space: Lagrange, degree: 1
/// - test space:  Lagrange, degree: 1
///
class p1_full_stokesvar_2_0_affine_q2 : public P1FormHyTeG
{

 public:

   p1_full_stokesvar_2_0_affine_q2() { WALBERLA_ABORT("Not implemented."); }

   p1_full_stokesvar_2_0_affine_q2( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
   : callback_Scalar_Variable_Coefficient_3D_mu(_callback_Scalar_Variable_Coefficient_3D_mu)
   {}

 private:

   std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_mu;


 public:

   /// \brief Not implemented - does nothing.
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const override;

   /// \brief Not implemented - does nothing.
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Vioreanu-Rokhlin 1 | points: 4, degree: 2, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              187     340       2       0      1            158                 4
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Vioreanu-Rokhlin 1 | points: 4, degree: 2, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              115     181       1       0      1             81                 4
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const override;

   bool assemble2D() const override { return false; }

   bool assembly2DDefined() const override { return false; }

   bool assemble3D() const override { return true; }

   bool assembly3DDefined() const override { return true; }

 private:

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p1_full_stokesvar_2_1_affine_q2
/// - description: 
/// - trial space: Lagrange, degree: 1
/// - test space:  Lagrange, degree: 1
///
class p1_full_stokesvar_2_1_affine_q2 : public P1FormHyTeG
{

 public:

   p1_full_stokesvar_2_1_affine_q2() { WALBERLA_ABORT("Not implemented."); }

   p1_full_stokesvar_2_1_affine_q2( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
   : callback_Scalar_Variable_Coefficient_3D_mu(_callback_Scalar_Variable_Coefficient_3D_mu)
   {}

 private:

   std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_mu;


 public:

   /// \brief Not implemented - does nothing.
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const override;

   /// \brief Not implemented - does nothing.
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Vioreanu-Rokhlin 1 | points: 4, degree: 2, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              187     340       2       0      1            158                 4
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Vioreanu-Rokhlin 1 | points: 4, degree: 2, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              115     181       1       0      1             81                 4
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const override;

   bool assemble2D() const override { return false; }

   bool assembly2DDefined() const override { return false; }

   bool assemble3D() const override { return true; }

   bool assembly3DDefined() const override { return true; }

 private:

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p1_full_stokesvar_2_2_affine_q2
/// - description: 
/// - trial space: Lagrange, degree: 1
/// - test space:  Lagrange, degree: 1
///
class p1_full_stokesvar_2_2_affine_q2 : public P1FormHyTeG
{

 public:

   p1_full_stokesvar_2_2_affine_q2() { WALBERLA_ABORT("Not implemented."); }

   p1_full_stokesvar_2_2_affine_q2( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
   : callback_Scalar_Variable_Coefficient_3D_mu(_callback_Scalar_Variable_Coefficient_3D_mu)
   {}

 private:

   std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_mu;


 public:

   /// \brief Not implemented - does nothing.
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const override;

   /// \brief Not implemented - does nothing.
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Vioreanu-Rokhlin 1 | points: 4, degree: 2, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              220     362       2       0      1            170                 4
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Vioreanu-Rokhlin 1 | points: 4, degree: 2, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                              154     237       1       0      1            105                 4
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const override;

   bool assemble2D() const override { return false; }

   bool assembly2DDefined() const override { return false; }

   bool assemble3D() const override { return true; }

   bool assembly3DDefined() const override { return true; }

 private:

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

} // namespace forms
} // namespace hyteg
