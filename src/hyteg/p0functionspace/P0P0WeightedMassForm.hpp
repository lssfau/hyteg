
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

# pragma once

# include "core/DataTypes.h"

# include "hyteg/dgfunctionspace/DGBasisInfo.hpp"
# include "hyteg/dgfunctionspace/DGForm.hpp"
# include "hyteg/dgfunctionspace/DGForm2D.hpp"
# include "hyteg/types/matrix.hpp"
# include "hyteg/types/pointnd.hpp"

# include "Eigen/Eigen"

namespace hyteg {
namespace dg{
namespace eg{

class P0P0WeightedMassForm : public hyteg::dg::DGForm
{

 public:
    P0P0WeightedMassForm(std::function< real_t ( const Point3D & ) > mu)
: callback_Scalar_Variable_Coefficient_3D_k (mu)
, callback_Scalar_Variable_Coefficient_2D_k (mu)
    {}

void Scalar_Variable_Coefficient_2D_k( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_k( Point3D( {in_0, in_1, 0} ) );
}

void Scalar_Variable_Coefficient_3D_k( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_k( Point3D( {in_0, in_1, in_2} ) );
}

 protected:
  void integrateVolume2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      real_t Scalar_Variable_Coefficient_2D_k_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id5 = 0;
      Scalar_Variable_Coefficient_2D_k( 0.091576213509770743*p_affine_0_0 + 0.091576213509770743*p_affine_1_0 + 0.81684757298045851*p_affine_2_0, 0.091576213509770743*p_affine_0_1 + 0.091576213509770743*p_affine_1_1 + 0.81684757298045851*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id0 );
      Scalar_Variable_Coefficient_2D_k( 0.44594849091596489*p_affine_0_0 + 0.44594849091596489*p_affine_1_0 + 0.10810301816807022*p_affine_2_0, 0.44594849091596489*p_affine_0_1 + 0.44594849091596489*p_affine_1_1 + 0.10810301816807022*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id1 );
      Scalar_Variable_Coefficient_2D_k( 0.091576213509770743*p_affine_0_0 + 0.81684757298045851*p_affine_1_0 + 0.091576213509770743*p_affine_2_0, 0.091576213509770743*p_affine_0_1 + 0.81684757298045851*p_affine_1_1 + 0.091576213509770743*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id2 );
      Scalar_Variable_Coefficient_2D_k( 0.44594849091596489*p_affine_0_0 + 0.10810301816807022*p_affine_1_0 + 0.44594849091596489*p_affine_2_0, 0.44594849091596489*p_affine_0_1 + 0.10810301816807022*p_affine_1_1 + 0.44594849091596489*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id3 );
      Scalar_Variable_Coefficient_2D_k( 0.81684757298045851*p_affine_0_0 + 0.091576213509770743*p_affine_1_0 + 0.091576213509770743*p_affine_2_0, 0.81684757298045851*p_affine_0_1 + 0.091576213509770743*p_affine_1_1 + 0.091576213509770743*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id4 );
      Scalar_Variable_Coefficient_2D_k( 0.10810301816807022*p_affine_0_0 + 0.44594849091596489*p_affine_1_0 + 0.44594849091596489*p_affine_2_0, 0.10810301816807022*p_affine_0_1 + 0.44594849091596489*p_affine_1_1 + 0.44594849091596489*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id5 );
      real_t tmp_0 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t a_0_0 = 0.11169079483900572*tmp_0/std::pow(Scalar_Variable_Coefficient_2D_k_out0_id5, 1.0/2.0) + 0.054975871827660928*tmp_0/std::pow(Scalar_Variable_Coefficient_2D_k_out0_id4, 1.0/2.0) + 0.11169079483900572*tmp_0/std::pow(Scalar_Variable_Coefficient_2D_k_out0_id3, 1.0/2.0) + 0.054975871827660928*tmp_0/std::pow(Scalar_Variable_Coefficient_2D_k_out0_id2, 1.0/2.0) + 0.11169079483900572*tmp_0/std::pow(Scalar_Variable_Coefficient_2D_k_out0_id1, 1.0/2.0) + 0.054975871827660928*tmp_0/std::pow(Scalar_Variable_Coefficient_2D_k_out0_id0, 1.0/2.0);
      elMat( 0, 0) = a_0_0;
   }

   virtual void integrateFacetInner2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                       const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coordsElement[0]( 0 );
      const auto p_affine_0_1 = coordsElement[0]( 1 );

      const auto p_affine_1_0 = coordsElement[1]( 0 );
      const auto p_affine_1_1 = coordsElement[1]( 1 );

      const auto p_affine_2_0 = coordsElement[2]( 0 );
      const auto p_affine_2_1 = coordsElement[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_8_0 = oppositeVertex( 0 );
      const auto p_affine_8_1 = oppositeVertex( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      elMat( 0, 0) = 0;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementInner,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementOuter,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexInnerElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexOuterElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                          const DGBasisInfo&                                       trialBasis,
                                          const DGBasisInfo&                                       testBasis,
                                          int                                                      trialDegree,
                                          int                                                      testDegree,
                                          Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coordsElementInner[0]( 0 );
      const auto p_affine_0_1 = coordsElementInner[0]( 1 );

      const auto p_affine_1_0 = coordsElementInner[1]( 0 );
      const auto p_affine_1_1 = coordsElementInner[1]( 1 );

      const auto p_affine_2_0 = coordsElementInner[2]( 0 );
      const auto p_affine_2_1 = coordsElementInner[2]( 1 );

      const auto p_affine_3_0 = coordsElementOuter[0]( 0 );
      const auto p_affine_3_1 = coordsElementOuter[0]( 1 );

      const auto p_affine_4_0 = coordsElementOuter[1]( 0 );
      const auto p_affine_4_1 = coordsElementOuter[1]( 1 );

      const auto p_affine_5_0 = coordsElementOuter[2]( 0 );
      const auto p_affine_5_1 = coordsElementOuter[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_8_0 = oppositeVertexInnerElement( 0 );
      const auto p_affine_8_1 = oppositeVertexInnerElement( 1 );

      const auto p_affine_9_0 = oppositeVertexOuterElement( 0 );
      const auto p_affine_9_1 = oppositeVertexOuterElement( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      elMat( 0, 0) = 0;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                   const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coordsElement[0]( 0 );
      const auto p_affine_0_1 = coordsElement[0]( 1 );

      const auto p_affine_1_0 = coordsElement[1]( 0 );
      const auto p_affine_1_1 = coordsElement[1]( 1 );

      const auto p_affine_2_0 = coordsElement[2]( 0 );
      const auto p_affine_2_1 = coordsElement[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_8_0 = oppositeVertex( 0 );
      const auto p_affine_8_1 = oppositeVertex( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      elMat( 0, 0) = 0;
   }

    void integrateRHSDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                 const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
     elMat.resize( Eigen::Index( basis.numDoFsPerElement( 3, walberla::uint_c( degree ) ) ), 1 );

      const auto p_affine_0_0 = coordsElement[0]( 0 );
      const auto p_affine_0_1 = coordsElement[0]( 1 );

      const auto p_affine_1_0 = coordsElement[1]( 0 );
      const auto p_affine_1_1 = coordsElement[1]( 1 );

      const auto p_affine_2_0 = coordsElement[2]( 0 );
      const auto p_affine_2_1 = coordsElement[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      elMat( 0, 0) = 0;
   }
   void integrateRHSDirichletBoundary3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                 const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
     elMat.resize( Eigen::Index( basis.numDoFsPerElement( 3, walberla::uint_c( degree ) ) ), 1 );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );

      elMat( 0, 0) = 0;
   }
   void integrateVolume3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
{
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coords[0]( 0 );
   const auto p_affine_0_1 = coords[0]( 1 );
   const auto p_affine_0_2 = coords[0]( 2 );

   const auto p_affine_1_0 = coords[1]( 0 );
   const auto p_affine_1_1 = coords[1]( 1 );
   const auto p_affine_1_2 = coords[1]( 2 );

   const auto p_affine_2_0 = coords[2]( 0 );
   const auto p_affine_2_1 = coords[2]( 1 );
   const auto p_affine_2_2 = coords[2]( 2 );

   const auto p_affine_3_0 = coords[3]( 0 );
   const auto p_affine_3_1 = coords[3]( 1 );
   const auto p_affine_3_2 = coords[3]( 2 );

      real_t Scalar_Variable_Coefficient_3D_k_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id5 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id6 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id7 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id8 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id9 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id10 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id11 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id12 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id13 = 0;
      Scalar_Variable_Coefficient_3D_k( 0.3108859192633005*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.067342242210098213*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.067342242210098213*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.067342242210098213*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id0 );
      Scalar_Variable_Coefficient_3D_k( 0.092735250310891248*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.72179424906732625*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.72179424906732625*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.72179424906732625*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id1 );
      Scalar_Variable_Coefficient_3D_k( 0.45449629587435036*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id2 );
      Scalar_Variable_Coefficient_3D_k( 0.045503704125649629*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.045503704125649629*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.045503704125649629*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id3 );
      Scalar_Variable_Coefficient_3D_k( 0.45449629587435036*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id4 );
      Scalar_Variable_Coefficient_3D_k( 0.45449629587435036*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id5 );
      Scalar_Variable_Coefficient_3D_k( 0.3108859192633005*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.067342242210098213*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.067342242210098213*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.067342242210098213*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id6 );
      Scalar_Variable_Coefficient_3D_k( 0.092735250310891248*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.72179424906732625*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.72179424906732625*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.72179424906732625*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id7 );
      Scalar_Variable_Coefficient_3D_k( 0.3108859192633005*p_affine_0_0 + 0.067342242210098213*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.067342242210098213*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.067342242210098213*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id8 );
      Scalar_Variable_Coefficient_3D_k( 0.092735250310891248*p_affine_0_0 + 0.72179424906732625*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.72179424906732625*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.72179424906732625*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id9 );
      Scalar_Variable_Coefficient_3D_k( 0.067342242210098102*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.067342242210098102*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.067342242210098102*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id10 );
      Scalar_Variable_Coefficient_3D_k( 0.72179424906732625*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.72179424906732625*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.72179424906732625*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id11 );
      Scalar_Variable_Coefficient_3D_k( 0.045503704125649636*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.045503704125649636*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.045503704125649636*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id12 );
      Scalar_Variable_Coefficient_3D_k( 0.045503704125649636*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.045503704125649636*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.045503704125649636*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id13 );
      real_t tmp_0 = p_affine_0_0*p_affine_1_1;
      real_t tmp_1 = p_affine_0_0*p_affine_1_2;
      real_t tmp_2 = p_affine_2_1*p_affine_3_2;
      real_t tmp_3 = p_affine_0_1*p_affine_1_0;
      real_t tmp_4 = p_affine_0_1*p_affine_1_2;
      real_t tmp_5 = p_affine_2_2*p_affine_3_0;
      real_t tmp_6 = p_affine_0_2*p_affine_1_0;
      real_t tmp_7 = p_affine_0_2*p_affine_1_1;
      real_t tmp_8 = p_affine_2_0*p_affine_3_1;
      real_t tmp_9 = p_affine_2_2*p_affine_3_1;
      real_t tmp_10 = p_affine_2_0*p_affine_3_2;
      real_t tmp_11 = p_affine_2_1*p_affine_3_0;
      real_t tmp_12 = 0.5*std::abs(p_affine_0_0*tmp_2 - p_affine_0_0*tmp_9 - p_affine_0_1*tmp_10 + p_affine_0_1*tmp_5 - p_affine_0_2*tmp_11 + p_affine_0_2*tmp_8 - p_affine_1_0*tmp_2 + p_affine_1_0*tmp_9 + p_affine_1_1*tmp_10 - p_affine_1_1*tmp_5 + p_affine_1_2*tmp_11 - p_affine_1_2*tmp_8 + p_affine_2_0*tmp_4 - p_affine_2_0*tmp_7 - p_affine_2_1*tmp_1 + p_affine_2_1*tmp_6 + p_affine_2_2*tmp_0 - p_affine_2_2*tmp_3 - p_affine_3_0*tmp_4 + p_affine_3_0*tmp_7 + p_affine_3_1*tmp_1 - p_affine_3_1*tmp_6 - p_affine_3_2*tmp_0 + p_affine_3_2*tmp_3);
      real_t a_0_0 = 0.012248840519393657*tmp_12/std::pow(Scalar_Variable_Coefficient_3D_k_out0_id9, 1.0/2.0) + 0.018781320953002646*tmp_12/std::pow(Scalar_Variable_Coefficient_3D_k_out0_id8, 1.0/2.0) + 0.012248840519393657*tmp_12/std::pow(Scalar_Variable_Coefficient_3D_k_out0_id7, 1.0/2.0) + 0.018781320953002646*tmp_12/std::pow(Scalar_Variable_Coefficient_3D_k_out0_id6, 1.0/2.0) + 0.0070910034628469103*tmp_12/std::pow(Scalar_Variable_Coefficient_3D_k_out0_id5, 1.0/2.0) + 0.0070910034628469103*tmp_12/std::pow(Scalar_Variable_Coefficient_3D_k_out0_id4, 1.0/2.0) + 0.0070910034628469103*tmp_12/std::pow(Scalar_Variable_Coefficient_3D_k_out0_id3, 1.0/2.0) + 0.0070910034628469103*tmp_12/std::pow(Scalar_Variable_Coefficient_3D_k_out0_id2, 1.0/2.0) + 0.0070910034628469103*tmp_12/std::pow(Scalar_Variable_Coefficient_3D_k_out0_id13, 1.0/2.0) + 0.0070910034628469103*tmp_12/std::pow(Scalar_Variable_Coefficient_3D_k_out0_id12, 1.0/2.0) + 0.012248840519393657*tmp_12/std::pow(Scalar_Variable_Coefficient_3D_k_out0_id11, 1.0/2.0) + 0.018781320953002646*tmp_12/std::pow(Scalar_Variable_Coefficient_3D_k_out0_id10, 1.0/2.0) + 0.012248840519393657*tmp_12/std::pow(Scalar_Variable_Coefficient_3D_k_out0_id1, 1.0/2.0) + 0.018781320953002646*tmp_12/std::pow(Scalar_Variable_Coefficient_3D_k_out0_id0, 1.0/2.0);
      elMat( 0, 0) = a_0_0;
   }



   void integrateFacetInner3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElement,
                                                     const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
                                                     const Eigen::Matrix< real_t, 3, 1 >&,
                                                     const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                     const DGBasisInfo&                                       trialBasis,
                                                     const DGBasisInfo&                                       testBasis,
                                                     int                                                      trialDegree,
                                                     int                                                      testDegree,
                                                     Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
{
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );

         elMat( 0, 0) = 0;
   }




void integrateFacetCoupling3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElementInner,
                                                        const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElementOuter,
                                                        const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
                                                        const Eigen::Matrix< real_t, 3, 1 >&,
                                                        const Eigen::Matrix< real_t, 3, 1 >&,
                                                        const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                        const DGBasisInfo&                                       trialBasis,
                                                        const DGBasisInfo&                                       testBasis,
                                                        int                                                      trialDegree,
                                                        int                                                      testDegree,
                                                        Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
{
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElementInner[0]( 0 );
   const auto p_affine_0_1 = coordsElementInner[0]( 1 );
   const auto p_affine_0_2 = coordsElementInner[0]( 2 );

   const auto p_affine_1_0 = coordsElementInner[1]( 0 );
   const auto p_affine_1_1 = coordsElementInner[1]( 1 );
   const auto p_affine_1_2 = coordsElementInner[1]( 2 );

   const auto p_affine_2_0 = coordsElementInner[2]( 0 );
   const auto p_affine_2_1 = coordsElementInner[2]( 1 );
   const auto p_affine_2_2 = coordsElementInner[2]( 2 );

   const auto p_affine_3_0 = coordsElementInner[3]( 0 );
   const auto p_affine_3_1 = coordsElementInner[3]( 1 );
   const auto p_affine_3_2 = coordsElementInner[3]( 2 );

   const auto p_affine_4_0 = coordsElementOuter[0]( 0 );
   const auto p_affine_4_1 = coordsElementOuter[0]( 1 );
   const auto p_affine_4_2 = coordsElementOuter[0]( 2 );

   const auto p_affine_5_0 = coordsElementOuter[1]( 0 );
   const auto p_affine_5_1 = coordsElementOuter[1]( 1 );
   const auto p_affine_5_2 = coordsElementOuter[1]( 2 );

   const auto p_affine_6_0 = coordsElementOuter[2]( 0 );
   const auto p_affine_6_1 = coordsElementOuter[2]( 1 );
   const auto p_affine_6_2 = coordsElementOuter[2]( 2 );

   const auto p_affine_7_0 = coordsElementOuter[3]( 0 );
   const auto p_affine_7_1 = coordsElementOuter[3]( 1 );
   const auto p_affine_7_2 = coordsElementOuter[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );


      elMat( 0, 0) = 0;
}



void integrateFacetDirichletBoundary3D(
    const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElement,
    const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
    const Eigen::Matrix< real_t, 3, 1 >&,
    const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
    const DGBasisInfo&                                       trialBasis,
    const DGBasisInfo&                                       testBasis,
    int                                                      trialDegree,
    int                                                      testDegree,
    Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
{
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );


      elMat( 0, 0) = 0;
   }

public:

std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_k;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_k;

};


} //eg
} // dg
} // hyteg
