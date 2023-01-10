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
 * Avoid modifying this file. If buggy, consider fixing the generator itself.
 */

#pragma once

#include "hyteg/geometry/GeometryMap.hpp"
#include "hyteg/forms/form_hyteg_base/P1FormHyTeG.hpp"

namespace hyteg {
namespace forms {

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p1_full_stokesvar_0_0_affine_q1
/// - description: 
/// - trial space: Lagrange, degree: 1
/// - test space:  Lagrange, degree: 1
///
class p1_full_stokesvar_0_0_affine_q1 : public P1FormHyTeG
{

 public:

   p1_full_stokesvar_0_0_affine_q1() { WALBERLA_ABORT("Not implemented."); }

   p1_full_stokesvar_0_0_affine_q1( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_2D_mu, std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
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
   /// - quadrature rule:                        Witherden-Vincent 1 | points: 1, degree: 1, test tolerance: 7.85e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               29      61       2       0      1             52                 1
   ///
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3
   /// - element matrix dimensions (rows, cols): (3, 3)
   /// - quadrature rule:                        Witherden-Vincent 1 | points: 1, degree: 1, test tolerance: 7.85e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               26      42       1       0      1             30                 1
   ///
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Xiao-Gimbutas 1 | points: 1, degree: 1, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               91     162       2       0      1            117                 1
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Xiao-Gimbutas 1 | points: 1, degree: 1, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               79     119       1       0      1             71                 1
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const override;

 private:

   void Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t * out_0 ) const;

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p1_full_stokesvar_0_1_affine_q1
/// - description: 
/// - trial space: Lagrange, degree: 1
/// - test space:  Lagrange, degree: 1
///
class p1_full_stokesvar_0_1_affine_q1 : public P1FormHyTeG
{

 public:

   p1_full_stokesvar_0_1_affine_q1() { WALBERLA_ABORT("Not implemented."); }

   p1_full_stokesvar_0_1_affine_q1( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_2D_mu, std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
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
   /// - quadrature rule:                        Witherden-Vincent 1 | points: 1, degree: 1, test tolerance: 7.85e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               27      63       2       0      1             52                 1
   ///
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3
   /// - element matrix dimensions (rows, cols): (3, 3)
   /// - quadrature rule:                        Witherden-Vincent 1 | points: 1, degree: 1, test tolerance: 7.85e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               23      37       1       0      1             26                 1
   ///
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Xiao-Gimbutas 1 | points: 1, degree: 1, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               73     143       2       0      1            109                 1
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Xiao-Gimbutas 1 | points: 1, degree: 1, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               64      96       1       0      1             62                 1
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const override;

 private:

   void Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t * out_0 ) const;

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p1_full_stokesvar_0_2_affine_q1
/// - description: 
/// - trial space: Lagrange, degree: 1
/// - test space:  Lagrange, degree: 1
///
class p1_full_stokesvar_0_2_affine_q1 : public P1FormHyTeG
{

 public:

   p1_full_stokesvar_0_2_affine_q1() { WALBERLA_ABORT("Not implemented."); }

   p1_full_stokesvar_0_2_affine_q1( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
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
   /// - quadrature rule:                        Xiao-Gimbutas 1 | points: 1, degree: 1, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               73     143       2       0      1            109                 1
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Xiao-Gimbutas 1 | points: 1, degree: 1, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               64      96       1       0      1             62                 1
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const override;

 private:

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p1_full_stokesvar_1_0_affine_q1
/// - description: 
/// - trial space: Lagrange, degree: 1
/// - test space:  Lagrange, degree: 1
///
class p1_full_stokesvar_1_0_affine_q1 : public P1FormHyTeG
{

 public:

   p1_full_stokesvar_1_0_affine_q1() { WALBERLA_ABORT("Not implemented."); }

   p1_full_stokesvar_1_0_affine_q1( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_2D_mu, std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
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
   /// - quadrature rule:                        Witherden-Vincent 1 | points: 1, degree: 1, test tolerance: 7.85e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               27      62       2       0      1             53                 1
   ///
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3
   /// - element matrix dimensions (rows, cols): (3, 3)
   /// - quadrature rule:                        Witherden-Vincent 1 | points: 1, degree: 1, test tolerance: 7.85e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               23      37       1       0      1             26                 1
   ///
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Xiao-Gimbutas 1 | points: 1, degree: 1, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               73     143       2       0      1            109                 1
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Xiao-Gimbutas 1 | points: 1, degree: 1, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               64      96       1       0      1             62                 1
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const override;

 private:

   void Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t * out_0 ) const;

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p1_full_stokesvar_1_1_affine_q1
/// - description: 
/// - trial space: Lagrange, degree: 1
/// - test space:  Lagrange, degree: 1
///
class p1_full_stokesvar_1_1_affine_q1 : public P1FormHyTeG
{

 public:

   p1_full_stokesvar_1_1_affine_q1() { WALBERLA_ABORT("Not implemented."); }

   p1_full_stokesvar_1_1_affine_q1( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_2D_mu, std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
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
   /// - quadrature rule:                        Witherden-Vincent 1 | points: 1, degree: 1, test tolerance: 7.85e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               29      61       2       0      1             52                 1
   ///
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       triangle, dim: 2, vertices: 3
   /// - element matrix dimensions (rows, cols): (3, 3)
   /// - quadrature rule:                        Witherden-Vincent 1 | points: 1, degree: 1, test tolerance: 7.85e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               26      42       1       0      1             30                 1
   ///
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Xiao-Gimbutas 1 | points: 1, degree: 1, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               91     162       2       0      1            117                 1
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Xiao-Gimbutas 1 | points: 1, degree: 1, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               79     119       1       0      1             71                 1
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const override;

 private:

   void Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t * out_0 ) const;

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p1_full_stokesvar_1_2_affine_q1
/// - description: 
/// - trial space: Lagrange, degree: 1
/// - test space:  Lagrange, degree: 1
///
class p1_full_stokesvar_1_2_affine_q1 : public P1FormHyTeG
{

 public:

   p1_full_stokesvar_1_2_affine_q1() { WALBERLA_ABORT("Not implemented."); }

   p1_full_stokesvar_1_2_affine_q1( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
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
   /// - quadrature rule:                        Xiao-Gimbutas 1 | points: 1, degree: 1, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               73     143       2       0      1            109                 1
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Xiao-Gimbutas 1 | points: 1, degree: 1, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               64      96       1       0      1             62                 1
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const override;

 private:

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p1_full_stokesvar_2_0_affine_q1
/// - description: 
/// - trial space: Lagrange, degree: 1
/// - test space:  Lagrange, degree: 1
///
class p1_full_stokesvar_2_0_affine_q1 : public P1FormHyTeG
{

 public:

   p1_full_stokesvar_2_0_affine_q1() { WALBERLA_ABORT("Not implemented."); }

   p1_full_stokesvar_2_0_affine_q1( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
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
   /// - quadrature rule:                        Xiao-Gimbutas 1 | points: 1, degree: 1, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               73     143       2       0      1            109                 1
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Xiao-Gimbutas 1 | points: 1, degree: 1, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               64      96       1       0      1             62                 1
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const override;

 private:

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p1_full_stokesvar_2_1_affine_q1
/// - description: 
/// - trial space: Lagrange, degree: 1
/// - test space:  Lagrange, degree: 1
///
class p1_full_stokesvar_2_1_affine_q1 : public P1FormHyTeG
{

 public:

   p1_full_stokesvar_2_1_affine_q1() { WALBERLA_ABORT("Not implemented."); }

   p1_full_stokesvar_2_1_affine_q1( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
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
   /// - quadrature rule:                        Xiao-Gimbutas 1 | points: 1, degree: 1, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               73     143       2       0      1            109                 1
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Xiao-Gimbutas 1 | points: 1, degree: 1, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               64      96       1       0      1             62                 1
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const override;

 private:

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

/// Implementation of the integration of a weak form over an element.
///
/// - name:        p1_full_stokesvar_2_2_affine_q1
/// - description: 
/// - trial space: Lagrange, degree: 1
/// - test space:  Lagrange, degree: 1
///
class p1_full_stokesvar_2_2_affine_q1 : public P1FormHyTeG
{

 public:

   p1_full_stokesvar_2_2_affine_q1() { WALBERLA_ABORT("Not implemented."); }

   p1_full_stokesvar_2_2_affine_q1( std::function< real_t ( const Point3D & ) > _callback_Scalar_Variable_Coefficient_3D_mu )
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
   /// - quadrature rule:                        Xiao-Gimbutas 1 | points: 1, degree: 1, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               91     162       2       0      1            117                 1
   ///
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const override;

   /// \brief Integrates the weak form over the passed element (vertices in computational space).
   ///
   /// - element geometry:                       tetrahedron, dim: 3, vertices: 4
   /// - element matrix dimensions (rows, cols): (4, 4)
   /// - quadrature rule:                        Xiao-Gimbutas 1 | points: 1, degree: 1, test tolerance: 2.379e-17
   /// - floating point operations:
   ///                                             adds    muls    divs    pows    abs    assignments    function_calls
   ///                                           ------  ------  ------  ------  -----  -------------  ----------------
   ///                                               79     119       1       0      1             71                 1
   ///
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const override;

 private:

   void Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const;

};

} // namespace forms
} // namespace hyteg
