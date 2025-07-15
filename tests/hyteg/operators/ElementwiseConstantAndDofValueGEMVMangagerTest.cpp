/*
* Copyright (c) 2025 Andreas Burkhart
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
#include <iostream>

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P1ToP2ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ToP1ElementwiseOperator.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/geometry/PolarCoordsMap.hpp"
#include "hyteg/geometry/SphericalCoordsMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/operators/ScaledOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include "constant_stencil_operator/P1ConstantOperator.hpp"
#include "constant_stencil_operator/P2ConstantOperator.hpp"
#include "elementwise_dof_value_operator/generated/p1_k_mass_visc_centroid_blending_q4_ElementwiseOperator.hpp"
#include "elementwise_dof_value_operator/generated/p1_to_p2_div_blending_q6_ElementwiseOperator.hpp"
#include "elementwise_dof_value_operator/generated/p2_k_mass_visc_centroid_blending_q6_ElementwiseOperator.hpp"
#include "elementwise_dof_value_operator/generated/p2_to_p1_div_blending_q6_ElementwiseOperator.hpp"
#include "mixed_operator/P1ToP2ConstantOperator.hpp"
#include "mixed_operator/P2ToP1ConstantOperator.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;
using namespace hyteg;

// Dummy class to check if toMatrix works properly in a scaled fashion.
// The dummy class just adds up the absolute value of all matrix entries.
// This avoids having e.g. a PETSc dependency here.
class SumMatrixProxy : public SparseMatrixProxy
{
 public:
   SumMatrixProxy( real_t val )
   : val_( val )
   {}

   virtual ~SumMatrixProxy() = default;

   std::shared_ptr< SparseMatrixProxy > createCopy() const override { return std::make_shared< SumMatrixProxy >( val_ ); }
   std::shared_ptr< SparseMatrixProxy > createEmptyCopy() const override
   {
      return std::make_shared< SumMatrixProxy >( real_c( 0 ) );
   }
   std::shared_ptr< SparseMatrixProxy > createMatrix( uint_t          localRows,
                                                      uint_t          localCols,
                                                      uint_t          globalRows,
                                                      uint_t          globalCols,
                                                      const MPI_Comm& MpiCommunicator ) const override
   {
      return std::make_shared< SumMatrixProxy >( real_c( 0 ) );
   }

   void addValue( uint_t row, uint_t col, real_t value ) override { val_ += abs( value ); }

   void addValues( const std::vector< uint_t >& rows,
                   const std::vector< uint_t >& cols,
                   const std::vector< real_t >& values ) override
   {
      WALBERLA_ASSERT_EQUAL( values.size(), rows.size() * cols.size() );

      for ( uint_t i = 0; i < rows.size(); i++ )
      {
         for ( uint_t j = 0; j < cols.size(); j++ )
         {
            val_ += abs( values[i * cols.size() + j] );
         }
      }
   }

   void createFromMatrixProduct( const std::vector< std::shared_ptr< SparseMatrixProxy > >& matrices ) override
   {
      val_ = real_c( 1 );
      for ( auto& m : matrices )
      {
         auto proxy = std::dynamic_pointer_cast< SumMatrixProxy >( m );
         val_ *= proxy->getVal();
      }
   }

   void createFromMatrixLinComb( const std::vector< real_t >&                               scalars,
                                 const std::vector< std::shared_ptr< SparseMatrixProxy > >& matrices ) override
   {
      val_ = real_c( 0 );
      for ( uint_t i = 0; i < matrices.size(); i++ )
      {
         auto proxy = std::dynamic_pointer_cast< SumMatrixProxy >( matrices[i] );
         val_ += proxy->getVal() * scalars[i];
      }
   }

   real_t getVal() { return val_; }
   void   setVal( real_t val ) { val_ = val; }

 private:
   real_t val_;
};

int main( int argc, char** argv )
{
   // Create walberla & MPI environment
   walberla::Environment env( argc, argv );
   walberla::mpi::MPIManager::instance()->useWorldComm();

   if ( sizeof( real_t ) < 8 )
   {
      // The tolerances in this test are designed for at least double precision.
      WALBERLA_LOG_INFO( "Single precision or lower detected. Aborting test." );

      return EXIT_SUCCESS;
   }

   {
      // P1 2D Scaling Test

      const uint_t minLevel        = 3;
      const uint_t maxLevel        = 3;
      const real_t tolerance       = real_c( 9e-4 );
      const real_t toleranceDiag   = real_c( 1e-15 );
      const real_t toleranceMatrix = real_c( 1e-15 );

      // Init setup storage
      MeshInfo meshInfo = MeshInfo::meshRectangle(
          Point2D( real_c( 0 ), real_c( 0 ) ), Point2D( real_c( 1 ), real_c( 1 ) ), MeshInfo::CROSS, 3, 3 );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

      // Create storage
      std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

      // Create functions
      P2Function< real_t > T( "T", storage, minLevel, maxLevel );
      P1Function< real_t > u( "u", storage, minLevel, maxLevel );
      P1Function< real_t > resConst( "resConst", storage, minLevel, maxLevel );
      P1Function< real_t > resElementwise( "resElementwise", storage, minLevel, maxLevel );
      P1Function< real_t > resElementwiseDoFValue( "resElementwiseDoFValue", storage, minLevel, maxLevel );
      P1Function< real_t > O( "O", storage, minLevel, maxLevel );

      std::function< real_t( const hyteg::Point3D& ) > uFct = []( const hyteg::Point3D& x ) {
         WALBERLA_UNUSED( x );
         return std::sin( x[0] ) + x[1] * x[1];
      };

      u.interpolate( uFct, maxLevel, hyteg::All );
      O.interpolate( real_c( 1 ), maxLevel, hyteg::All );

      // operators
      std::function< real_t( const Point3D&, real_t, real_t, real_t, real_t, real_t ) > massScaling =
          [=]( const hyteg::Point3D& x, real_t temp, real_t tempCentroid, real_t centX, real_t centY, real_t centZ ) {
             return real_c( 1 );
          };

      auto CentroidMass = std::make_shared< p1_k_mass_visc_centroid_blending_q4_ElementwiseOperator >(
          storage, minLevel, maxLevel, T, massScaling );
      auto P1ElementwiseMass = std::make_shared< P1ElementwiseMassOperator >( storage, minLevel, maxLevel );
      auto P1ConstantMass    = std::make_shared< P1ConstantMassOperator >( storage, minLevel, maxLevel );

      CentroidMass->apply( u, resElementwiseDoFValue, maxLevel, All, Replace );
      P1ElementwiseMass->apply( u, resElementwise, maxLevel, All, Replace );
      P1ConstantMass->apply( u, resConst, maxLevel, All, Replace );

      real_t integValDoFValue = O.dotGlobal( resElementwiseDoFValue, maxLevel, All );
      real_t integValElem     = O.dotGlobal( resElementwise, maxLevel, All );
      real_t integValConst    = O.dotGlobal( resConst, maxLevel, All );
      real_t integValCtrl     = real_c( 4 ) / real_c( 3 ) - std::cos( 1 );

      WALBERLA_LOG_INFO( abs( integValDoFValue - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValDoFValue - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      auto CentroidMassScaled = std::make_shared<
          ScaledOperator< p1_k_mass_visc_centroid_blending_q4_ElementwiseOperator, P1Function< real_t >, P1Function< real_t > > >(
          storage, minLevel, maxLevel, CentroidMass, real_c( 2 ) );
      auto P1ElementwiseMassScaled =
          std::make_shared< ScaledOperator< P1ElementwiseMassOperator, P1Function< real_t >, P1Function< real_t > > >(
              storage, minLevel, maxLevel, P1ElementwiseMass, real_c( 2 ) );
      auto P1ConstantMassScaled =
          std::make_shared< ScaledOperator< P1ConstantMassOperator, P1Function< real_t >, P1Function< real_t > > >(
              storage, minLevel, maxLevel, P1ConstantMass, real_c( 2 ) );

      CentroidMassScaled->apply( u, resElementwiseDoFValue, maxLevel, All, Replace );
      P1ElementwiseMassScaled->apply( u, resElementwise, maxLevel, All, Replace );
      P1ConstantMassScaled->apply( u, resConst, maxLevel, All, Replace );

      integValDoFValue = O.dotGlobal( resElementwiseDoFValue, maxLevel, All );
      integValElem     = O.dotGlobal( resElementwise, maxLevel, All );
      integValConst    = O.dotGlobal( resConst, maxLevel, All );
      integValCtrl     = real_c( 8 ) / real_c( 3 ) - real_c( 2 ) * std::cos( 1 );

      WALBERLA_LOG_INFO( abs( integValDoFValue - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValDoFValue - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      auto elementwiseDoFType = hyteg::applyGEMV< p1_k_mass_visc_centroid_blending_q4_ElementwiseOperator >(
          *CentroidMass, real_c( 2 ), u, real_c( 1 ), resElementwiseDoFValue, maxLevel, All );
      auto elementwiseType = hyteg::applyGEMV< P1ElementwiseMassOperator >(
          *P1ElementwiseMass, real_c( 2 ), u, real_c( 1 ), resElementwise, maxLevel, All );
      auto constantType =
          hyteg::applyGEMV< P1ConstantMassOperator >( *P1ConstantMass, real_c( 2 ), u, real_c( 1 ), resConst, maxLevel, All );

      integValDoFValue = O.dotGlobal( resElementwiseDoFValue, maxLevel, All );
      integValElem     = O.dotGlobal( resElementwise, maxLevel, All );
      integValConst    = O.dotGlobal( resConst, maxLevel, All );
      integValCtrl     = real_c( 16 ) / real_c( 3 ) - real_c( 4 ) * std::cos( 1 );

      WALBERLA_LOG_INFO( abs( integValDoFValue - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValDoFValue - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled, smooth_jac_scaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseDoFType );
      WALBERLA_CHECK_EQUAL( elementwiseDoFType, GEMVType::GEMV );
      WALBERLA_LOG_INFO( elementwiseType );
      WALBERLA_CHECK_EQUAL( elementwiseType, GEMVType::GEMV );
      WALBERLA_LOG_INFO( constantType );
      WALBERLA_CHECK_EQUAL( constantType, GEMVType::MANUAL );

      // Now check the inverse diagonal
      P1Function< real_t > elementwiseDoFDiagStore( "elementwiseDoFDiagStore", storage, minLevel, maxLevel );
      P1Function< real_t > elementwiseDiagStore( "elementwiseDiagStore", storage, minLevel, maxLevel );
      P1Function< real_t > constDiagStore( "constDiagStore", storage, minLevel, maxLevel );

      CentroidMass->computeInverseDiagonalOperatorValues();
      P1ElementwiseMass->computeInverseDiagonalOperatorValues();
      P1ConstantMass->computeInverseDiagonalOperatorValues();

      elementwiseDoFDiagStore.assign( { real_c( 1 ) }, { *CentroidMass->getInverseDiagonalValues() }, maxLevel, All );
      elementwiseDiagStore.assign( { real_c( 1 ) }, { *P1ElementwiseMass->getInverseDiagonalValues() }, maxLevel, All );
      constDiagStore.assign( { real_c( 1 ) }, { *P1ConstantMass->getInverseDiagonalValues() }, maxLevel, All );

      CentroidMassScaled->computeInverseDiagonalOperatorValues();
      P1ElementwiseMassScaled->computeInverseDiagonalOperatorValues();
      P1ConstantMassScaled->computeInverseDiagonalOperatorValues();

      CentroidMass->getInverseDiagonalValues()->assign(
          { real_c( 1 ), real_c( -0.25 ), real_c( -0.25 ) },
          { *CentroidMass->getInverseDiagonalValues(), elementwiseDoFDiagStore, elementwiseDoFDiagStore },
          maxLevel,
          All );

      P1ElementwiseMass->getInverseDiagonalValues()->assign(
          { real_c( 1 ), real_c( -0.25 ), real_c( -0.25 ) },
          { *P1ElementwiseMass->getInverseDiagonalValues(), elementwiseDiagStore, elementwiseDiagStore },
          maxLevel,
          All );

      P1ConstantMassScaled->getInverseDiagonalValues()->assign(
          { real_c( 1 ), real_c( -0.25 ), real_c( -0.25 ) },
          { *P1ConstantMassScaled->getInverseDiagonalValues(), constDiagStore, constDiagStore },
          maxLevel,
          All );

      real_t diagValDoFElem =
          CentroidMass->getInverseDiagonalValues()->dotGlobal( *CentroidMass->getInverseDiagonalValues(), maxLevel, All );
      real_t diagValElem = P1ElementwiseMass->getInverseDiagonalValues()->dotGlobal(
          *P1ElementwiseMass->getInverseDiagonalValues(), maxLevel, All );
      real_t diagValConst = P1ConstantMassScaled->getInverseDiagonalValues()->dotGlobal(
          *P1ConstantMassScaled->getInverseDiagonalValues(), maxLevel, All );

      WALBERLA_LOG_INFO( abs( diagValDoFElem ) );
      WALBERLA_CHECK_LESS( abs( diagValDoFElem ), toleranceDiag );
      WALBERLA_LOG_INFO( abs( diagValElem ) );
      WALBERLA_CHECK_LESS( abs( diagValElem ), toleranceDiag );
      WALBERLA_LOG_INFO( abs( diagValConst ) );
      WALBERLA_CHECK_LESS( abs( diagValConst ), toleranceDiag );

      auto elementwiseDoFDiagType =
          hyteg::applyComputeInverseDiagonalOperatorValuesScaled< p1_k_mass_visc_centroid_blending_q4_ElementwiseOperator >(
              *CentroidMass, real_c( 2 ) );
      auto elementwiseDiagType =
          hyteg::applyComputeInverseDiagonalOperatorValuesScaled< P1ElementwiseMassOperator >( *P1ElementwiseMass, real_c( 2 ) );
      auto constantDiagType =
          hyteg::applyComputeInverseDiagonalOperatorValuesScaled< P1ConstantMassOperator >( *P1ConstantMass, real_c( 2 ) );

      CentroidMass->getInverseDiagonalValues()->assign(
          { real_c( 1 ), real_c( -0.25 ), real_c( -0.25 ) },
          { *CentroidMass->getInverseDiagonalValues(), elementwiseDoFDiagStore, elementwiseDoFDiagStore },
          maxLevel,
          All );

      P1ElementwiseMass->getInverseDiagonalValues()->assign(
          { real_c( 1 ), real_c( -0.25 ), real_c( -0.25 ) },
          { *P1ElementwiseMass->getInverseDiagonalValues(), elementwiseDiagStore, elementwiseDiagStore },
          maxLevel,
          All );

      P1ConstantMassScaled->getInverseDiagonalValues()->assign(
          { real_c( 1 ), real_c( -0.25 ), real_c( -0.25 ) },
          { *P1ConstantMassScaled->getInverseDiagonalValues(), constDiagStore, constDiagStore },
          maxLevel,
          All );

      diagValDoFElem =
          CentroidMass->getInverseDiagonalValues()->dotGlobal( *CentroidMass->getInverseDiagonalValues(), maxLevel, All );
      diagValElem = P1ElementwiseMass->getInverseDiagonalValues()->dotGlobal(
          *P1ElementwiseMass->getInverseDiagonalValues(), maxLevel, All );
      diagValConst = P1ConstantMassScaled->getInverseDiagonalValues()->dotGlobal(
          *P1ConstantMassScaled->getInverseDiagonalValues(), maxLevel, All );

      WALBERLA_LOG_INFO( abs( diagValDoFElem ) );
      WALBERLA_CHECK_LESS( abs( diagValDoFElem ), toleranceDiag );
      WALBERLA_LOG_INFO( abs( diagValElem ) );
      WALBERLA_CHECK_LESS( abs( diagValElem ), toleranceDiag );
      WALBERLA_LOG_INFO( abs( diagValConst ) );
      WALBERLA_CHECK_LESS( abs( diagValConst ), toleranceDiag );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled, smooth_jac_scaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseDoFDiagType );
      WALBERLA_CHECK_EQUAL( elementwiseDoFDiagType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( elementwiseDiagType );
      WALBERLA_CHECK_EQUAL( elementwiseDiagType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( constantDiagType );
      WALBERLA_CHECK_EQUAL( constantDiagType, GEMVType::MANUAL );

      // create dummy matrix proxy
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyDoFElem;
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyElem;
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyConst;

      dummyMatrixProxyDoFElem = std::make_shared< SumMatrixProxy >( real_c( 0 ) );
      dummyMatrixProxyElem    = std::make_shared< SumMatrixProxy >( real_c( 0 ) );
      dummyMatrixProxyConst   = std::make_shared< SumMatrixProxy >( real_c( 0 ) );

      P1Function< idx_t > enumerator( "enumerator", storage, minLevel, maxLevel );
      enumerator.enumerate( maxLevel );

      CentroidMass->toMatrix( dummyMatrixProxyDoFElem, enumerator, enumerator, maxLevel, All );
      P1ElementwiseMass->toMatrix( dummyMatrixProxyElem, enumerator, enumerator, maxLevel, All );
      P1ConstantMass->toMatrix( dummyMatrixProxyConst, enumerator, enumerator, maxLevel, All );

      real_t matrixSumValDoFElem = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal();
      real_t matrixSumValElem    = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal();
      real_t matrixSumValConst   = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal();

      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->setVal( real_c( 0 ) );

      CentroidMassScaled->toMatrix( dummyMatrixProxyDoFElem, enumerator, enumerator, maxLevel, All );
      P1ElementwiseMassScaled->toMatrix( dummyMatrixProxyElem, enumerator, enumerator, maxLevel, All );
      P1ConstantMassScaled->toMatrix( dummyMatrixProxyConst, enumerator, enumerator, maxLevel, All );

      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                              matrixSumValDoFElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                                matrixSumValDoFElem * real_c( 2 ) ),
                           toleranceMatrix );
      WALBERLA_LOG_INFO(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ),
          toleranceMatrix );
      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() -
                              matrixSumValConst * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() - matrixSumValConst * real_c( 2 ) ),
          toleranceMatrix );

      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->setVal( real_c( 0 ) );

      auto elementwiseDoFMatrixType = hyteg::applyToMatrixScaled(
          *CentroidMass, real_c( 2 ), dummyMatrixProxyDoFElem, enumerator, enumerator, maxLevel, All );
      auto elementwiseMatrixType = hyteg::applyToMatrixScaled(
          *P1ElementwiseMass, real_c( 2 ), dummyMatrixProxyElem, enumerator, enumerator, maxLevel, All );
      auto constMatrixType = hyteg::applyToMatrixScaled(
          *P1ConstantMass, real_c( 2 ), dummyMatrixProxyConst, enumerator, enumerator, maxLevel, All );

      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                              matrixSumValDoFElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                                matrixSumValDoFElem * real_c( 2 ) ),
                           toleranceMatrix );
      WALBERLA_LOG_INFO(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ),
          toleranceMatrix );
      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() -
                              matrixSumValConst * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() - matrixSumValConst * real_c( 2 ) ),
          toleranceMatrix );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled, smooth_jac_scaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseDoFMatrixType );
      WALBERLA_CHECK_EQUAL( elementwiseDoFMatrixType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( elementwiseMatrixType );
      WALBERLA_CHECK_EQUAL( elementwiseMatrixType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( constMatrixType );
      WALBERLA_CHECK_EQUAL( constMatrixType, GEMVType::MANUAL );
   }

   {
      // P1 3D Scaling Test

      const uint_t minLevel        = 2;
      const uint_t maxLevel        = 2;
      const real_t tolerance       = real_c( 9e-4 );
      const real_t toleranceDiag   = real_c( 1e-15 );
      const real_t toleranceMatrix = real_c( 1e-15 );

      // Init setup storage
      MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid(
          Point3D( real_c( 0 ), real_c( 0 ), real_c( 0 ) ), Point3D( real_c( 1 ), real_c( 1 ), real_c( 1 ) ), 3, 3, 3 );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

      // Create storage
      std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

      // Create functions
      P2Function< real_t > T( "T", storage, minLevel, maxLevel );
      P1Function< real_t > u( "u", storage, minLevel, maxLevel );
      P1Function< real_t > resConst( "resConst", storage, minLevel, maxLevel );
      P1Function< real_t > resElementwise( "resElementwise", storage, minLevel, maxLevel );
      P1Function< real_t > resElementwiseDoFValue( "resElementwiseDoFValue", storage, minLevel, maxLevel );
      P1Function< real_t > O( "O", storage, minLevel, maxLevel );

      std::function< real_t( const hyteg::Point3D& ) > uFct = []( const hyteg::Point3D& x ) {
         WALBERLA_UNUSED( x );
         return std::sin( x[0] ) + x[1] * x[1] + cos( x[2] );
      };

      u.interpolate( uFct, maxLevel, hyteg::All );
      O.interpolate( real_c( 1 ), maxLevel, hyteg::All );

      // operators
      std::function< real_t( const Point3D&, real_t, real_t, real_t, real_t, real_t ) > massScaling =
          [=]( const hyteg::Point3D& x, real_t temp, real_t tempCentroid, real_t centX, real_t centY, real_t centZ ) {
             return real_c( 1 );
          };

      auto CentroidMass = std::make_shared< p1_k_mass_visc_centroid_blending_q4_ElementwiseOperator >(
          storage, minLevel, maxLevel, T, massScaling );
      auto P1ElementwiseMass = std::make_shared< P1ElementwiseMassOperator >( storage, minLevel, maxLevel );
      auto P1ConstantMass    = std::make_shared< P1ConstantMassOperator >( storage, minLevel, maxLevel );

      CentroidMass->apply( u, resElementwiseDoFValue, maxLevel, All, Replace );
      P1ElementwiseMass->apply( u, resElementwise, maxLevel, All, Replace );
      P1ConstantMass->apply( u, resConst, maxLevel, All, Replace );

      real_t integValDoFValue = O.dotGlobal( resElementwiseDoFValue, maxLevel, All );
      real_t integValElem     = O.dotGlobal( resElementwise, maxLevel, All );
      real_t integValConst    = O.dotGlobal( resConst, maxLevel, All );
      real_t integValCtrl     = real_c( 4 ) / real_c( 3 ) + std::sin( 1 ) - std::cos( 1 );

      WALBERLA_LOG_INFO( abs( integValDoFValue - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValDoFValue - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      auto CentroidMassScaled = std::make_shared<
          ScaledOperator< p1_k_mass_visc_centroid_blending_q4_ElementwiseOperator, P1Function< real_t >, P1Function< real_t > > >(
          storage, minLevel, maxLevel, CentroidMass, real_c( 2 ) );
      auto P1ElementwiseMassScaled =
          std::make_shared< ScaledOperator< P1ElementwiseMassOperator, P1Function< real_t >, P1Function< real_t > > >(
              storage, minLevel, maxLevel, P1ElementwiseMass, real_c( 2 ) );
      auto P1ConstantMassScaled =
          std::make_shared< ScaledOperator< P1ConstantMassOperator, P1Function< real_t >, P1Function< real_t > > >(
              storage, minLevel, maxLevel, P1ConstantMass, real_c( 2 ) );

      CentroidMassScaled->apply( u, resElementwiseDoFValue, maxLevel, All, Replace );
      P1ElementwiseMassScaled->apply( u, resElementwise, maxLevel, All, Replace );
      P1ConstantMassScaled->apply( u, resConst, maxLevel, All, Replace );

      integValDoFValue = O.dotGlobal( resElementwiseDoFValue, maxLevel, All );
      integValElem     = O.dotGlobal( resElementwise, maxLevel, All );
      integValConst    = O.dotGlobal( resConst, maxLevel, All );
      integValCtrl     = real_c( 8 ) / real_c( 3 ) + real_c( 2 ) * std::sin( 1 ) - real_c( 2 ) * std::cos( 1 );

      WALBERLA_LOG_INFO( abs( integValDoFValue - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValDoFValue - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      auto elementwiseDoFType = hyteg::applyGEMV< p1_k_mass_visc_centroid_blending_q4_ElementwiseOperator >(
          *CentroidMass, real_c( 2 ), u, real_c( 1 ), resElementwiseDoFValue, maxLevel, All );
      auto elementwiseType = hyteg::applyGEMV< P1ElementwiseMassOperator >(
          *P1ElementwiseMass, real_c( 2 ), u, real_c( 1 ), resElementwise, maxLevel, All );
      auto constantType =
          hyteg::applyGEMV< P1ConstantMassOperator >( *P1ConstantMass, real_c( 2 ), u, real_c( 1 ), resConst, maxLevel, All );

      integValDoFValue = O.dotGlobal( resElementwiseDoFValue, maxLevel, All );
      integValElem     = O.dotGlobal( resElementwise, maxLevel, All );
      integValConst    = O.dotGlobal( resConst, maxLevel, All );
      integValCtrl     = real_c( 16 ) / real_c( 3 ) + real_c( 4 ) * std::sin( 1 ) - real_c( 4 ) * std::cos( 1 );

      WALBERLA_LOG_INFO( abs( integValDoFValue - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValDoFValue - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled, smooth_jac_scaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseDoFType );
      WALBERLA_CHECK_EQUAL( elementwiseDoFType, GEMVType::GEMV );
      WALBERLA_LOG_INFO( elementwiseType );
      WALBERLA_CHECK_EQUAL( elementwiseType, GEMVType::GEMV );
      WALBERLA_LOG_INFO( constantType );
      WALBERLA_CHECK_EQUAL( constantType, GEMVType::MANUAL );

      // Now check the inverse diagonal
      P1Function< real_t > elementwiseDoFDiagStore( "elementwiseDoFDiagStore", storage, minLevel, maxLevel );
      P1Function< real_t > elementwiseDiagStore( "elementwiseDiagStore", storage, minLevel, maxLevel );
      P1Function< real_t > constDiagStore( "constDiagStore", storage, minLevel, maxLevel );

      CentroidMass->computeInverseDiagonalOperatorValues();
      P1ElementwiseMass->computeInverseDiagonalOperatorValues();
      P1ConstantMass->computeInverseDiagonalOperatorValues();

      elementwiseDoFDiagStore.assign( { real_c( 1 ) }, { *CentroidMass->getInverseDiagonalValues() }, maxLevel, All );
      elementwiseDiagStore.assign( { real_c( 1 ) }, { *P1ElementwiseMass->getInverseDiagonalValues() }, maxLevel, All );
      constDiagStore.assign( { real_c( 1 ) }, { *P1ConstantMass->getInverseDiagonalValues() }, maxLevel, All );

      CentroidMassScaled->computeInverseDiagonalOperatorValues();
      P1ElementwiseMassScaled->computeInverseDiagonalOperatorValues();
      P1ConstantMassScaled->computeInverseDiagonalOperatorValues();

      CentroidMass->getInverseDiagonalValues()->assign(
          { real_c( 1 ), real_c( -0.25 ), real_c( -0.25 ) },
          { *CentroidMass->getInverseDiagonalValues(), elementwiseDoFDiagStore, elementwiseDoFDiagStore },
          maxLevel,
          All );

      P1ElementwiseMass->getInverseDiagonalValues()->assign(
          { real_c( 1 ), real_c( -0.25 ), real_c( -0.25 ) },
          { *P1ElementwiseMass->getInverseDiagonalValues(), elementwiseDiagStore, elementwiseDiagStore },
          maxLevel,
          All );

      P1ConstantMassScaled->getInverseDiagonalValues()->assign(
          { real_c( 1 ), real_c( -0.25 ), real_c( -0.25 ) },
          { *P1ConstantMassScaled->getInverseDiagonalValues(), constDiagStore, constDiagStore },
          maxLevel,
          All );

      real_t diagValDoFElem =
          CentroidMass->getInverseDiagonalValues()->dotGlobal( *CentroidMass->getInverseDiagonalValues(), maxLevel, All );
      real_t diagValElem = P1ElementwiseMass->getInverseDiagonalValues()->dotGlobal(
          *P1ElementwiseMass->getInverseDiagonalValues(), maxLevel, All );
      real_t diagValConst = P1ConstantMassScaled->getInverseDiagonalValues()->dotGlobal(
          *P1ConstantMassScaled->getInverseDiagonalValues(), maxLevel, All );

      WALBERLA_LOG_INFO( abs( diagValDoFElem ) );
      WALBERLA_CHECK_LESS( abs( diagValDoFElem ), toleranceDiag );
      WALBERLA_LOG_INFO( abs( diagValElem ) );
      WALBERLA_CHECK_LESS( abs( diagValElem ), toleranceDiag );
      WALBERLA_LOG_INFO( abs( diagValConst ) );
      WALBERLA_CHECK_LESS( abs( diagValConst ), toleranceDiag );

      auto elementwiseDoFDiagType =
          hyteg::applyComputeInverseDiagonalOperatorValuesScaled< p1_k_mass_visc_centroid_blending_q4_ElementwiseOperator >(
              *CentroidMass, real_c( 2 ) );
      auto elementwiseDiagType =
          hyteg::applyComputeInverseDiagonalOperatorValuesScaled< P1ElementwiseMassOperator >( *P1ElementwiseMass, real_c( 2 ) );
      auto constantDiagType =
          hyteg::applyComputeInverseDiagonalOperatorValuesScaled< P1ConstantMassOperator >( *P1ConstantMass, real_c( 2 ) );

      CentroidMass->getInverseDiagonalValues()->assign(
          { real_c( 1 ), real_c( -0.25 ), real_c( -0.25 ) },
          { *CentroidMass->getInverseDiagonalValues(), elementwiseDoFDiagStore, elementwiseDoFDiagStore },
          maxLevel,
          All );

      P1ElementwiseMass->getInverseDiagonalValues()->assign(
          { real_c( 1 ), real_c( -0.25 ), real_c( -0.25 ) },
          { *P1ElementwiseMass->getInverseDiagonalValues(), elementwiseDiagStore, elementwiseDiagStore },
          maxLevel,
          All );

      P1ConstantMassScaled->getInverseDiagonalValues()->assign(
          { real_c( 1 ), real_c( -0.25 ), real_c( -0.25 ) },
          { *P1ConstantMassScaled->getInverseDiagonalValues(), constDiagStore, constDiagStore },
          maxLevel,
          All );

      diagValDoFElem =
          CentroidMass->getInverseDiagonalValues()->dotGlobal( *CentroidMass->getInverseDiagonalValues(), maxLevel, All );
      diagValElem = P1ElementwiseMass->getInverseDiagonalValues()->dotGlobal(
          *P1ElementwiseMass->getInverseDiagonalValues(), maxLevel, All );
      diagValConst = P1ConstantMassScaled->getInverseDiagonalValues()->dotGlobal(
          *P1ConstantMassScaled->getInverseDiagonalValues(), maxLevel, All );

      WALBERLA_LOG_INFO( abs( diagValDoFElem ) );
      WALBERLA_CHECK_LESS( abs( diagValDoFElem ), toleranceDiag );
      WALBERLA_LOG_INFO( abs( diagValElem ) );
      WALBERLA_CHECK_LESS( abs( diagValElem ), toleranceDiag );
      WALBERLA_LOG_INFO( abs( diagValConst ) );
      WALBERLA_CHECK_LESS( abs( diagValConst ), toleranceDiag );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled, smooth_jac_scaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseDoFDiagType );
      WALBERLA_CHECK_EQUAL( elementwiseDoFDiagType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( elementwiseDiagType );
      WALBERLA_CHECK_EQUAL( elementwiseDiagType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( constantDiagType );
      WALBERLA_CHECK_EQUAL( constantDiagType, GEMVType::MANUAL );

      // create dummy matrix proxy
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyDoFElem;
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyElem;
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyConst;

      dummyMatrixProxyDoFElem = std::make_shared< SumMatrixProxy >( real_c( 0 ) );
      dummyMatrixProxyElem    = std::make_shared< SumMatrixProxy >( real_c( 0 ) );
      dummyMatrixProxyConst   = std::make_shared< SumMatrixProxy >( real_c( 0 ) );

      P1Function< idx_t > enumerator( "enumerator", storage, minLevel, maxLevel );
      enumerator.enumerate( maxLevel );

      CentroidMass->toMatrix( dummyMatrixProxyDoFElem, enumerator, enumerator, maxLevel, All );
      P1ElementwiseMass->toMatrix( dummyMatrixProxyElem, enumerator, enumerator, maxLevel, All );
      P1ConstantMass->toMatrix( dummyMatrixProxyConst, enumerator, enumerator, maxLevel, All );

      real_t matrixSumValDoFElem = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal();
      real_t matrixSumValElem    = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal();
      real_t matrixSumValConst   = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal();

      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->setVal( real_c( 0 ) );

      CentroidMassScaled->toMatrix( dummyMatrixProxyDoFElem, enumerator, enumerator, maxLevel, All );
      P1ElementwiseMassScaled->toMatrix( dummyMatrixProxyElem, enumerator, enumerator, maxLevel, All );
      P1ConstantMassScaled->toMatrix( dummyMatrixProxyConst, enumerator, enumerator, maxLevel, All );

      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                              matrixSumValDoFElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                                matrixSumValDoFElem * real_c( 2 ) ),
                           toleranceMatrix );
      WALBERLA_LOG_INFO(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ),
          toleranceMatrix );
      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() -
                              matrixSumValConst * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() - matrixSumValConst * real_c( 2 ) ),
          toleranceMatrix );

      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->setVal( real_c( 0 ) );

      auto elementwiseDoFMatrixType = hyteg::applyToMatrixScaled(
          *CentroidMass, real_c( 2 ), dummyMatrixProxyDoFElem, enumerator, enumerator, maxLevel, All );
      auto elementwiseMatrixType = hyteg::applyToMatrixScaled(
          *P1ElementwiseMass, real_c( 2 ), dummyMatrixProxyElem, enumerator, enumerator, maxLevel, All );
      auto constMatrixType = hyteg::applyToMatrixScaled(
          *P1ConstantMass, real_c( 2 ), dummyMatrixProxyConst, enumerator, enumerator, maxLevel, All );

      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                              matrixSumValDoFElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                                matrixSumValDoFElem * real_c( 2 ) ),
                           toleranceMatrix );
      WALBERLA_LOG_INFO(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ),
          toleranceMatrix );
      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() -
                              matrixSumValConst * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() - matrixSumValConst * real_c( 2 ) ),
          toleranceMatrix );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled, smooth_jac_scaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseDoFMatrixType );
      WALBERLA_CHECK_EQUAL( elementwiseDoFMatrixType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( elementwiseMatrixType );
      WALBERLA_CHECK_EQUAL( elementwiseMatrixType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( constMatrixType );
      WALBERLA_CHECK_EQUAL( constMatrixType, GEMVType::MANUAL );
   }

   {
      // P2 2D Scaling Test

      const uint_t minLevel        = 3;
      const uint_t maxLevel        = 3;
      const real_t tolerance       = real_c( 2e-9 );
      const real_t toleranceDiag   = real_c( 1e-15 );
      const real_t toleranceMatrix = real_c( 1e-15 );

      // Init setup storage
      MeshInfo meshInfo = MeshInfo::meshRectangle(
          Point2D( real_c( 0 ), real_c( 0 ) ), Point2D( real_c( 1 ), real_c( 1 ) ), MeshInfo::CROSS, 3, 3 );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

      // Create storage
      std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

      // Create functions
      P2Function< real_t > T( "T", storage, minLevel, maxLevel );
      P2Function< real_t > u( "u", storage, minLevel, maxLevel );
      P2Function< real_t > resConst( "resConst", storage, minLevel, maxLevel );
      P2Function< real_t > resElementwise( "resElementwise", storage, minLevel, maxLevel );
      P2Function< real_t > resElementwiseDoFValue( "resElementwiseDoFValue", storage, minLevel, maxLevel );
      P2Function< real_t > O( "O", storage, minLevel, maxLevel );

      std::function< real_t( const hyteg::Point3D& ) > uFct = []( const hyteg::Point3D& x ) {
         WALBERLA_UNUSED( x );
         return std::sin( x[0] ) + x[1] * x[1];
      };

      u.interpolate( uFct, maxLevel, hyteg::All );
      O.interpolate( real_c( 1 ), maxLevel, hyteg::All );

      // operators
      std::function< real_t( const Point3D&, real_t, real_t, real_t, real_t, real_t ) > massScaling =
          [=]( const hyteg::Point3D& x, real_t temp, real_t tempCentroid, real_t centX, real_t centY, real_t centZ ) {
             return real_c( 1 );
          };

      auto CentroidMass = std::make_shared< p2_k_mass_visc_centroid_blending_q6_ElementwiseOperator >(
          storage, minLevel, maxLevel, T, massScaling );
      auto P2ElementwiseMass = std::make_shared< P2ElementwiseMassOperator >( storage, minLevel, maxLevel );
      auto P2ConstantMass    = std::make_shared< P2ConstantMassOperator >( storage, minLevel, maxLevel );

      CentroidMass->apply( u, resElementwiseDoFValue, maxLevel, All, Replace );
      P2ElementwiseMass->apply( u, resElementwise, maxLevel, All, Replace );
      P2ConstantMass->apply( u, resConst, maxLevel, All, Replace );

      real_t integValDoFValue = O.dotGlobal( resElementwiseDoFValue, maxLevel, All );
      real_t integValElem     = O.dotGlobal( resElementwise, maxLevel, All );
      real_t integValConst    = O.dotGlobal( resConst, maxLevel, All );
      real_t integValCtrl     = real_c( 4 ) / real_c( 3 ) - std::cos( 1 );

      WALBERLA_LOG_INFO( abs( integValDoFValue - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValDoFValue - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      auto CentroidMassScaled = std::make_shared<
          ScaledOperator< p2_k_mass_visc_centroid_blending_q6_ElementwiseOperator, P2Function< real_t >, P2Function< real_t > > >(
          storage, minLevel, maxLevel, CentroidMass, real_c( 2 ) );
      auto P2ElementwiseMassScaled =
          std::make_shared< ScaledOperator< P2ElementwiseMassOperator, P2Function< real_t >, P2Function< real_t > > >(
              storage, minLevel, maxLevel, P2ElementwiseMass, real_c( 2 ) );
      auto P2ConstantMassScaled =
          std::make_shared< ScaledOperator< P2ConstantMassOperator, P2Function< real_t >, P2Function< real_t > > >(
              storage, minLevel, maxLevel, P2ConstantMass, real_c( 2 ) );

      CentroidMassScaled->apply( u, resElementwiseDoFValue, maxLevel, All, Replace );
      P2ElementwiseMassScaled->apply( u, resElementwise, maxLevel, All, Replace );
      P2ConstantMassScaled->apply( u, resConst, maxLevel, All, Replace );

      integValDoFValue = O.dotGlobal( resElementwiseDoFValue, maxLevel, All );
      integValElem     = O.dotGlobal( resElementwise, maxLevel, All );
      integValConst    = O.dotGlobal( resConst, maxLevel, All );
      integValCtrl     = real_c( 8 ) / real_c( 3 ) - real_c( 2 ) * std::cos( 1 );

      WALBERLA_LOG_INFO( abs( integValDoFValue - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValDoFValue - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      auto elementwiseDoFType = hyteg::applyGEMV< p2_k_mass_visc_centroid_blending_q6_ElementwiseOperator >(
          *CentroidMass, real_c( 2 ), u, real_c( 1 ), resElementwiseDoFValue, maxLevel, All );
      auto elementwiseType = hyteg::applyGEMV< P2ElementwiseMassOperator >(
          *P2ElementwiseMass, real_c( 2 ), u, real_c( 1 ), resElementwise, maxLevel, All );
      auto constantType =
          hyteg::applyGEMV< P2ConstantMassOperator >( *P2ConstantMass, real_c( 2 ), u, real_c( 1 ), resConst, maxLevel, All );

      integValDoFValue = O.dotGlobal( resElementwiseDoFValue, maxLevel, All );
      integValElem     = O.dotGlobal( resElementwise, maxLevel, All );
      integValConst    = O.dotGlobal( resConst, maxLevel, All );
      integValCtrl     = real_c( 16 ) / real_c( 3 ) - real_c( 4 ) * std::cos( 1 );

      WALBERLA_LOG_INFO( abs( integValDoFValue - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValDoFValue - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled, smooth_jac_scaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseDoFType );
      WALBERLA_CHECK_EQUAL( elementwiseDoFType, GEMVType::GEMV );
      WALBERLA_LOG_INFO( elementwiseType );
      WALBERLA_CHECK_EQUAL( elementwiseType, GEMVType::GEMV );
      WALBERLA_LOG_INFO( constantType );
      WALBERLA_CHECK_EQUAL( constantType, GEMVType::MANUAL );

      // Now check the inverse diagonal, P2 constant operators for some reason do not have a computeInverseDiagonalOperatorValues method.
      // Leaving the respective calls as a comment in case they can be reactivated later on.
      P2Function< real_t > elementwiseDoFDiagStore( "elementwiseDoFDiagStore", storage, minLevel, maxLevel );
      P2Function< real_t > elementwiseDiagStore( "elementwiseDiagStore", storage, minLevel, maxLevel );
      //   P2Function< real_t > constDiagStore( "constDiagStore", storage, minLevel, maxLevel );

      CentroidMass->computeInverseDiagonalOperatorValues();
      P2ElementwiseMass->computeInverseDiagonalOperatorValues();
      //   P2ConstantMass->computeInverseDiagonalOperatorValues();

      elementwiseDoFDiagStore.assign( { real_c( 1 ) }, { *CentroidMass->getInverseDiagonalValues() }, maxLevel, All );
      elementwiseDiagStore.assign( { real_c( 1 ) }, { *P2ElementwiseMass->getInverseDiagonalValues() }, maxLevel, All );
      //   constDiagStore.assign( { real_c( 1 ) }, { *P2ConstantMass->getInverseDiagonalValues() }, maxLevel, All );

      CentroidMassScaled->computeInverseDiagonalOperatorValues();
      P2ElementwiseMassScaled->computeInverseDiagonalOperatorValues();
      //   P2ConstantMassScaled->computeInverseDiagonalOperatorValues();

      CentroidMass->getInverseDiagonalValues()->assign(
          { real_c( 1 ), real_c( -0.25 ), real_c( -0.25 ) },
          { *CentroidMass->getInverseDiagonalValues(), elementwiseDoFDiagStore, elementwiseDoFDiagStore },
          maxLevel,
          All );

      P2ElementwiseMass->getInverseDiagonalValues()->assign(
          { real_c( 1 ), real_c( -0.25 ), real_c( -0.25 ) },
          { *P2ElementwiseMass->getInverseDiagonalValues(), elementwiseDiagStore, elementwiseDiagStore },
          maxLevel,
          All );

      //   P2ConstantMassScaled->getInverseDiagonalValues()->assign(
      //       { real_c( 1 ), real_c( -0.25 ), real_c( -0.25 ) },
      //       { *P2ConstantMassScaled->getInverseDiagonalValues(), constDiagStore, constDiagStore },
      //       maxLevel,
      //       All );

      real_t diagValDoFElem =
          CentroidMass->getInverseDiagonalValues()->dotGlobal( *CentroidMass->getInverseDiagonalValues(), maxLevel, All );
      real_t diagValElem = P2ElementwiseMass->getInverseDiagonalValues()->dotGlobal(
          *P2ElementwiseMass->getInverseDiagonalValues(), maxLevel, All );
      //   real_t diagValConst = P2ConstantMassScaled->getInverseDiagonalValues()->dotGlobal(
      //       *P2ConstantMassScaled->getInverseDiagonalValues(), maxLevel, All );

      WALBERLA_LOG_INFO( abs( diagValDoFElem ) );
      WALBERLA_CHECK_LESS( abs( diagValDoFElem ), toleranceDiag );
      WALBERLA_LOG_INFO( abs( diagValElem ) );
      WALBERLA_CHECK_LESS( abs( diagValElem ), toleranceDiag );
      //   WALBERLA_LOG_INFO( abs( diagValConst ) );
      //   WALBERLA_CHECK_LESS( abs( diagValConst ), toleranceDiag );

      auto elementwiseDoFDiagType =
          hyteg::applyComputeInverseDiagonalOperatorValuesScaled< p2_k_mass_visc_centroid_blending_q6_ElementwiseOperator >(
              *CentroidMass, real_c( 2 ) );
      auto elementwiseDiagType =
          hyteg::applyComputeInverseDiagonalOperatorValuesScaled< P2ElementwiseMassOperator >( *P2ElementwiseMass, real_c( 2 ) );
      //   auto constantDiagType =
      //       hyteg::applyComputeInverseDiagonalOperatorValuesScaled< P2ConstantMassOperator >( *P2ConstantMass, real_c( 2 ) );

      CentroidMass->getInverseDiagonalValues()->assign(
          { real_c( 1 ), real_c( -0.25 ), real_c( -0.25 ) },
          { *CentroidMass->getInverseDiagonalValues(), elementwiseDoFDiagStore, elementwiseDoFDiagStore },
          maxLevel,
          All );

      P2ElementwiseMass->getInverseDiagonalValues()->assign(
          { real_c( 1 ), real_c( -0.25 ), real_c( -0.25 ) },
          { *P2ElementwiseMass->getInverseDiagonalValues(), elementwiseDiagStore, elementwiseDiagStore },
          maxLevel,
          All );

      //   P2ConstantMassScaled->getInverseDiagonalValues()->assign(
      //       { real_c( 1 ), real_c( -0.25 ), real_c( -0.25 ) },
      //       { *P2ConstantMassScaled->getInverseDiagonalValues(), constDiagStore, constDiagStore },
      //       maxLevel,
      //       All );

      diagValDoFElem =
          CentroidMass->getInverseDiagonalValues()->dotGlobal( *CentroidMass->getInverseDiagonalValues(), maxLevel, All );
      diagValElem = P2ElementwiseMass->getInverseDiagonalValues()->dotGlobal(
          *P2ElementwiseMass->getInverseDiagonalValues(), maxLevel, All );
      //   diagValConst = P2ConstantMassScaled->getInverseDiagonalValues()->dotGlobal(
      //       *P2ConstantMassScaled->getInverseDiagonalValues(), maxLevel, All );

      WALBERLA_LOG_INFO( abs( diagValDoFElem ) );
      WALBERLA_CHECK_LESS( abs( diagValDoFElem ), toleranceDiag );
      WALBERLA_LOG_INFO( abs( diagValElem ) );
      WALBERLA_CHECK_LESS( abs( diagValElem ), toleranceDiag );
      //   WALBERLA_LOG_INFO( abs( diagValConst ) );
      //   WALBERLA_CHECK_LESS( abs( diagValConst ), toleranceDiag );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled, smooth_jac_scaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseDoFDiagType );
      WALBERLA_CHECK_EQUAL( elementwiseDoFDiagType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( elementwiseDiagType );
      WALBERLA_CHECK_EQUAL( elementwiseDiagType, GEMVType::SCALED );
      //   WALBERLA_LOG_INFO( constantDiagType );
      //   WALBERLA_CHECK_EQUAL( constantDiagType, GEMVType::MANUAL );

      // create dummy matrix proxy
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyDoFElem;
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyElem;
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyConst;

      dummyMatrixProxyDoFElem = std::make_shared< SumMatrixProxy >( real_c( 0 ) );
      dummyMatrixProxyElem    = std::make_shared< SumMatrixProxy >( real_c( 0 ) );
      dummyMatrixProxyConst   = std::make_shared< SumMatrixProxy >( real_c( 0 ) );

      P2Function< idx_t > enumerator( "enumerator", storage, minLevel, maxLevel );
      enumerator.enumerate( maxLevel );

      CentroidMass->toMatrix( dummyMatrixProxyDoFElem, enumerator, enumerator, maxLevel, All );
      P2ElementwiseMass->toMatrix( dummyMatrixProxyElem, enumerator, enumerator, maxLevel, All );
      P2ConstantMass->toMatrix( dummyMatrixProxyConst, enumerator, enumerator, maxLevel, All );

      real_t matrixSumValDoFElem = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal();
      real_t matrixSumValElem    = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal();
      real_t matrixSumValConst   = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal();

      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->setVal( real_c( 0 ) );

      CentroidMassScaled->toMatrix( dummyMatrixProxyDoFElem, enumerator, enumerator, maxLevel, All );
      P2ElementwiseMassScaled->toMatrix( dummyMatrixProxyElem, enumerator, enumerator, maxLevel, All );
      P2ConstantMassScaled->toMatrix( dummyMatrixProxyConst, enumerator, enumerator, maxLevel, All );

      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                              matrixSumValDoFElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                                matrixSumValDoFElem * real_c( 2 ) ),
                           toleranceMatrix );
      WALBERLA_LOG_INFO(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ),
          toleranceMatrix );
      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() -
                              matrixSumValConst * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() - matrixSumValConst * real_c( 2 ) ),
          toleranceMatrix );

      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->setVal( real_c( 0 ) );

      auto elementwiseDoFMatrixType = hyteg::applyToMatrixScaled(
          *CentroidMass, real_c( 2 ), dummyMatrixProxyDoFElem, enumerator, enumerator, maxLevel, All );
      auto elementwiseMatrixType = hyteg::applyToMatrixScaled(
          *P2ElementwiseMass, real_c( 2 ), dummyMatrixProxyElem, enumerator, enumerator, maxLevel, All );
      auto constMatrixType = hyteg::applyToMatrixScaled(
          *P2ConstantMass, real_c( 2 ), dummyMatrixProxyConst, enumerator, enumerator, maxLevel, All );

      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                              matrixSumValDoFElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                                matrixSumValDoFElem * real_c( 2 ) ),
                           toleranceMatrix );
      WALBERLA_LOG_INFO(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ),
          toleranceMatrix );
      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() -
                              matrixSumValConst * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() - matrixSumValConst * real_c( 2 ) ),
          toleranceMatrix );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled, smooth_jac_scaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseDoFMatrixType );
      WALBERLA_CHECK_EQUAL( elementwiseDoFMatrixType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( elementwiseMatrixType );
      WALBERLA_CHECK_EQUAL( elementwiseMatrixType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( constMatrixType );
      WALBERLA_CHECK_EQUAL( constMatrixType, GEMVType::MANUAL );
   }

   {
      // P2 3D Scaling Test

      const uint_t minLevel        = 2;
      const uint_t maxLevel        = 2;
      const real_t tolerance       = real_c( 5e-9 );
      const real_t toleranceDiag   = real_c( 1e-15 );
      const real_t toleranceMatrix = real_c( 1e-15 );

      // Init setup storage
      MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid(
          Point3D( real_c( 0 ), real_c( 0 ), real_c( 0 ) ), Point3D( real_c( 1 ), real_c( 1 ), real_c( 1 ) ), 3, 3, 3 );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

      // Create storage
      std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

      // Create functions
      P2Function< real_t > T( "T", storage, minLevel, maxLevel );
      P2Function< real_t > u( "u", storage, minLevel, maxLevel );
      P2Function< real_t > resConst( "resConst", storage, minLevel, maxLevel );
      P2Function< real_t > resElementwise( "resElementwise", storage, minLevel, maxLevel );
      P2Function< real_t > resElementwiseDoFValue( "resElementwiseDoFValue", storage, minLevel, maxLevel );
      P2Function< real_t > O( "O", storage, minLevel, maxLevel );

      std::function< real_t( const hyteg::Point3D& ) > uFct = []( const hyteg::Point3D& x ) {
         WALBERLA_UNUSED( x );
         return std::sin( x[0] ) + x[1] * x[1] + cos( x[2] );
      };

      u.interpolate( uFct, maxLevel, hyteg::All );
      O.interpolate( real_c( 1 ), maxLevel, hyteg::All );

      // operators
      std::function< real_t( const Point3D&, real_t, real_t, real_t, real_t, real_t ) > massScaling =
          [=]( const hyteg::Point3D& x, real_t temp, real_t tempCentroid, real_t centX, real_t centY, real_t centZ ) {
             return real_c( 1 );
          };

      auto CentroidMass = std::make_shared< p2_k_mass_visc_centroid_blending_q6_ElementwiseOperator >(
          storage, minLevel, maxLevel, T, massScaling );
      auto P2ElementwiseMass = std::make_shared< P2ElementwiseMassOperator >( storage, minLevel, maxLevel );
      auto P2ConstantMass    = std::make_shared< P2ConstantMassOperator >( storage, minLevel, maxLevel );

      CentroidMass->apply( u, resElementwiseDoFValue, maxLevel, All, Replace );
      P2ElementwiseMass->apply( u, resElementwise, maxLevel, All, Replace );
      P2ConstantMass->apply( u, resConst, maxLevel, All, Replace );

      real_t integValDoFValue = O.dotGlobal( resElementwiseDoFValue, maxLevel, All );
      real_t integValElem     = O.dotGlobal( resElementwise, maxLevel, All );
      real_t integValConst    = O.dotGlobal( resConst, maxLevel, All );
      real_t integValCtrl     = real_c( 4 ) / real_c( 3 ) + std::sin( 1 ) - std::cos( 1 );

      WALBERLA_LOG_INFO( abs( integValDoFValue - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValDoFValue - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      auto CentroidMassScaled = std::make_shared<
          ScaledOperator< p2_k_mass_visc_centroid_blending_q6_ElementwiseOperator, P2Function< real_t >, P2Function< real_t > > >(
          storage, minLevel, maxLevel, CentroidMass, real_c( 2 ) );
      auto P2ElementwiseMassScaled =
          std::make_shared< ScaledOperator< P2ElementwiseMassOperator, P2Function< real_t >, P2Function< real_t > > >(
              storage, minLevel, maxLevel, P2ElementwiseMass, real_c( 2 ) );
      auto P2ConstantMassScaled =
          std::make_shared< ScaledOperator< P2ConstantMassOperator, P2Function< real_t >, P2Function< real_t > > >(
              storage, minLevel, maxLevel, P2ConstantMass, real_c( 2 ) );

      CentroidMassScaled->apply( u, resElementwiseDoFValue, maxLevel, All, Replace );
      P2ElementwiseMassScaled->apply( u, resElementwise, maxLevel, All, Replace );
      P2ConstantMassScaled->apply( u, resConst, maxLevel, All, Replace );

      integValDoFValue = O.dotGlobal( resElementwiseDoFValue, maxLevel, All );
      integValElem     = O.dotGlobal( resElementwise, maxLevel, All );
      integValConst    = O.dotGlobal( resConst, maxLevel, All );
      integValCtrl     = real_c( 8 ) / real_c( 3 ) + real_c( 2 ) * std::sin( 1 ) - real_c( 2 ) * std::cos( 1 );

      WALBERLA_LOG_INFO( abs( integValDoFValue - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValDoFValue - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      auto elementwiseDoFType = hyteg::applyGEMV< p2_k_mass_visc_centroid_blending_q6_ElementwiseOperator >(
          *CentroidMass, real_c( 2 ), u, real_c( 1 ), resElementwiseDoFValue, maxLevel, All );
      auto elementwiseType = hyteg::applyGEMV< P2ElementwiseMassOperator >(
          *P2ElementwiseMass, real_c( 2 ), u, real_c( 1 ), resElementwise, maxLevel, All );
      auto constantType =
          hyteg::applyGEMV< P2ConstantMassOperator >( *P2ConstantMass, real_c( 2 ), u, real_c( 1 ), resConst, maxLevel, All );

      integValDoFValue = O.dotGlobal( resElementwiseDoFValue, maxLevel, All );
      integValElem     = O.dotGlobal( resElementwise, maxLevel, All );
      integValConst    = O.dotGlobal( resConst, maxLevel, All );
      integValCtrl     = real_c( 16 ) / real_c( 3 ) + real_c( 4 ) * std::sin( 1 ) - real_c( 4 ) * std::cos( 1 );

      WALBERLA_LOG_INFO( abs( integValDoFValue - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValDoFValue - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled, smooth_jac_scaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseDoFType );
      WALBERLA_CHECK_EQUAL( elementwiseDoFType, GEMVType::GEMV );
      WALBERLA_LOG_INFO( elementwiseType );
      WALBERLA_CHECK_EQUAL( elementwiseType, GEMVType::GEMV );
      WALBERLA_LOG_INFO( constantType );
      WALBERLA_CHECK_EQUAL( constantType, GEMVType::MANUAL );

      // Now check the inverse diagonal, P2 constant operators for some reason do not have a computeInverseDiagonalOperatorValues method.
      // Leaving the respective calls as a comment in case they can be reactivated later on.
      P2Function< real_t > elementwiseDoFDiagStore( "elementwiseDoFDiagStore", storage, minLevel, maxLevel );
      P2Function< real_t > elementwiseDiagStore( "elementwiseDiagStore", storage, minLevel, maxLevel );
      //   P2Function< real_t > constDiagStore( "constDiagStore", storage, minLevel, maxLevel );

      CentroidMass->computeInverseDiagonalOperatorValues();
      P2ElementwiseMass->computeInverseDiagonalOperatorValues();
      //   P2ConstantMass->computeInverseDiagonalOperatorValues();

      elementwiseDoFDiagStore.assign( { real_c( 1 ) }, { *CentroidMass->getInverseDiagonalValues() }, maxLevel, All );
      elementwiseDiagStore.assign( { real_c( 1 ) }, { *P2ElementwiseMass->getInverseDiagonalValues() }, maxLevel, All );
      //   constDiagStore.assign( { real_c( 1 ) }, { *P2ConstantMass->getInverseDiagonalValues() }, maxLevel, All );

      CentroidMassScaled->computeInverseDiagonalOperatorValues();
      P2ElementwiseMassScaled->computeInverseDiagonalOperatorValues();
      //   P2ConstantMassScaled->computeInverseDiagonalOperatorValues();

      CentroidMass->getInverseDiagonalValues()->assign(
          { real_c( 1 ), real_c( -0.25 ), real_c( -0.25 ) },
          { *CentroidMass->getInverseDiagonalValues(), elementwiseDoFDiagStore, elementwiseDoFDiagStore },
          maxLevel,
          All );

      P2ElementwiseMass->getInverseDiagonalValues()->assign(
          { real_c( 1 ), real_c( -0.25 ), real_c( -0.25 ) },
          { *P2ElementwiseMass->getInverseDiagonalValues(), elementwiseDiagStore, elementwiseDiagStore },
          maxLevel,
          All );

      //   P2ConstantMassScaled->getInverseDiagonalValues()->assign(
      //       { real_c( 1 ), real_c( -0.25 ), real_c( -0.25 ) },
      //       { *P2ConstantMassScaled->getInverseDiagonalValues(), constDiagStore, constDiagStore },
      //       maxLevel,
      //       All );

      real_t diagValDoFElem =
          CentroidMass->getInverseDiagonalValues()->dotGlobal( *CentroidMass->getInverseDiagonalValues(), maxLevel, All );
      real_t diagValElem = P2ElementwiseMass->getInverseDiagonalValues()->dotGlobal(
          *P2ElementwiseMass->getInverseDiagonalValues(), maxLevel, All );
      //   real_t diagValConst = P2ConstantMassScaled->getInverseDiagonalValues()->dotGlobal(
      //       *P2ConstantMassScaled->getInverseDiagonalValues(), maxLevel, All );

      WALBERLA_LOG_INFO( abs( diagValDoFElem ) );
      WALBERLA_CHECK_LESS( abs( diagValDoFElem ), toleranceDiag );
      WALBERLA_LOG_INFO( abs( diagValElem ) );
      WALBERLA_CHECK_LESS( abs( diagValElem ), toleranceDiag );
      //   WALBERLA_LOG_INFO( abs( diagValConst ) );
      //   WALBERLA_CHECK_LESS( abs( diagValConst ), toleranceDiag );

      auto elementwiseDoFDiagType =
          hyteg::applyComputeInverseDiagonalOperatorValuesScaled< p2_k_mass_visc_centroid_blending_q6_ElementwiseOperator >(
              *CentroidMass, real_c( 2 ) );
      auto elementwiseDiagType =
          hyteg::applyComputeInverseDiagonalOperatorValuesScaled< P2ElementwiseMassOperator >( *P2ElementwiseMass, real_c( 2 ) );
      //   auto constantDiagType =
      //       hyteg::applyComputeInverseDiagonalOperatorValuesScaled< P2ConstantMassOperator >( *P2ConstantMass, real_c( 2 ) );

      CentroidMass->getInverseDiagonalValues()->assign(
          { real_c( 1 ), real_c( -0.25 ), real_c( -0.25 ) },
          { *CentroidMass->getInverseDiagonalValues(), elementwiseDoFDiagStore, elementwiseDoFDiagStore },
          maxLevel,
          All );

      P2ElementwiseMass->getInverseDiagonalValues()->assign(
          { real_c( 1 ), real_c( -0.25 ), real_c( -0.25 ) },
          { *P2ElementwiseMass->getInverseDiagonalValues(), elementwiseDiagStore, elementwiseDiagStore },
          maxLevel,
          All );

      //   P2ConstantMassScaled->getInverseDiagonalValues()->assign(
      //       { real_c( 1 ), real_c( -0.25 ), real_c( -0.25 ) },
      //       { *P2ConstantMassScaled->getInverseDiagonalValues(), constDiagStore, constDiagStore },
      //       maxLevel,
      //       All );

      diagValDoFElem =
          CentroidMass->getInverseDiagonalValues()->dotGlobal( *CentroidMass->getInverseDiagonalValues(), maxLevel, All );
      diagValElem = P2ElementwiseMass->getInverseDiagonalValues()->dotGlobal(
          *P2ElementwiseMass->getInverseDiagonalValues(), maxLevel, All );
      //   diagValConst = P2ConstantMassScaled->getInverseDiagonalValues()->dotGlobal(
      //       *P2ConstantMassScaled->getInverseDiagonalValues(), maxLevel, All );

      WALBERLA_LOG_INFO( abs( diagValDoFElem ) );
      WALBERLA_CHECK_LESS( abs( diagValDoFElem ), toleranceDiag );
      WALBERLA_LOG_INFO( abs( diagValElem ) );
      WALBERLA_CHECK_LESS( abs( diagValElem ), toleranceDiag );
      //   WALBERLA_LOG_INFO( abs( diagValConst ) );
      //   WALBERLA_CHECK_LESS( abs( diagValConst ), toleranceDiag );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled, smooth_jac_scaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseDoFDiagType );
      WALBERLA_CHECK_EQUAL( elementwiseDoFDiagType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( elementwiseDiagType );
      WALBERLA_CHECK_EQUAL( elementwiseDiagType, GEMVType::SCALED );
      //   WALBERLA_LOG_INFO( constantDiagType );
      //   WALBERLA_CHECK_EQUAL( constantDiagType, GEMVType::MANUAL );

      // create dummy matrix proxy
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyDoFElem;
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyElem;
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyConst;

      dummyMatrixProxyDoFElem = std::make_shared< SumMatrixProxy >( real_c( 0 ) );
      dummyMatrixProxyElem    = std::make_shared< SumMatrixProxy >( real_c( 0 ) );
      dummyMatrixProxyConst   = std::make_shared< SumMatrixProxy >( real_c( 0 ) );

      P2Function< idx_t > enumerator( "enumerator", storage, minLevel, maxLevel );
      enumerator.enumerate( maxLevel );

      CentroidMass->toMatrix( dummyMatrixProxyDoFElem, enumerator, enumerator, maxLevel, All );
      P2ElementwiseMass->toMatrix( dummyMatrixProxyElem, enumerator, enumerator, maxLevel, All );
      P2ConstantMass->toMatrix( dummyMatrixProxyConst, enumerator, enumerator, maxLevel, All );

      real_t matrixSumValDoFElem = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal();
      real_t matrixSumValElem    = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal();
      real_t matrixSumValConst   = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal();

      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->setVal( real_c( 0 ) );

      CentroidMassScaled->toMatrix( dummyMatrixProxyDoFElem, enumerator, enumerator, maxLevel, All );
      P2ElementwiseMassScaled->toMatrix( dummyMatrixProxyElem, enumerator, enumerator, maxLevel, All );
      P2ConstantMassScaled->toMatrix( dummyMatrixProxyConst, enumerator, enumerator, maxLevel, All );

      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                              matrixSumValDoFElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                                matrixSumValDoFElem * real_c( 2 ) ),
                           toleranceMatrix );
      WALBERLA_LOG_INFO(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ),
          toleranceMatrix );
      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() -
                              matrixSumValConst * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() - matrixSumValConst * real_c( 2 ) ),
          toleranceMatrix );

      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->setVal( real_c( 0 ) );

      auto elementwiseDoFMatrixType = hyteg::applyToMatrixScaled(
          *CentroidMass, real_c( 2 ), dummyMatrixProxyDoFElem, enumerator, enumerator, maxLevel, All );
      auto elementwiseMatrixType = hyteg::applyToMatrixScaled(
          *P2ElementwiseMass, real_c( 2 ), dummyMatrixProxyElem, enumerator, enumerator, maxLevel, All );
      auto constMatrixType = hyteg::applyToMatrixScaled(
          *P2ConstantMass, real_c( 2 ), dummyMatrixProxyConst, enumerator, enumerator, maxLevel, All );

      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                              matrixSumValDoFElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                                matrixSumValDoFElem * real_c( 2 ) ),
                           toleranceMatrix );
      WALBERLA_LOG_INFO(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ),
          toleranceMatrix );
      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() -
                              matrixSumValConst * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() - matrixSumValConst * real_c( 2 ) ),
          toleranceMatrix );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled, smooth_jac_scaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseDoFMatrixType );
      WALBERLA_CHECK_EQUAL( elementwiseDoFMatrixType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( elementwiseMatrixType );
      WALBERLA_CHECK_EQUAL( elementwiseMatrixType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( constMatrixType );
      WALBERLA_CHECK_EQUAL( constMatrixType, GEMVType::MANUAL );
   }

   {
      // P1 To P2 2D Scaling Test

      const uint_t minLevel        = 3;
      const uint_t maxLevel        = 3;
      const real_t tolerance       = real_c( 2e-12 );
      const real_t toleranceMatrix = real_c( 1e-15 );

      // Init setup storage
      MeshInfo meshInfo = MeshInfo::meshRectangle(
          Point2D( real_c( 0 ), real_c( 0 ) ), Point2D( real_c( 1 ), real_c( 1 ) ), MeshInfo::CROSS, 3, 3 );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

      // Create storage
      std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

      // Create functions
      P2Function< real_t > T( "T", storage, minLevel, maxLevel );
      P1Function< real_t > u( "u", storage, minLevel, maxLevel );
      P2Function< real_t > resConst( "resConst", storage, minLevel, maxLevel );
      P2Function< real_t > resElementwise( "resElementwise", storage, minLevel, maxLevel );
      P2Function< real_t > resElementwiseDoFValue( "resElementwiseDoFValue", storage, minLevel, maxLevel );
      P2Function< real_t > O( "O", storage, minLevel, maxLevel );

      std::function< real_t( const hyteg::Point3D& ) > uFct = []( const hyteg::Point3D& x ) {
         WALBERLA_UNUSED( x );
         return std::sin( x[0] ) + x[1] * x[1];
      };

      u.interpolate( uFct, maxLevel, hyteg::All );
      O.interpolate( real_c( 1 ), maxLevel, hyteg::All );

      // operators
      std::function< real_t( const Point3D&, real_t, real_t, real_t, real_t, real_t ) > massScaling =
          [=]( const hyteg::Point3D& x, real_t temp, real_t tempCentroid, real_t centX, real_t centY, real_t centZ ) {
             return real_c( 1 );
          };

      auto divXElemDoF = std::make_shared< p1_to_p2_div_0_blending_q6_ElementwiseOperator >( storage, minLevel, maxLevel, T );
      auto divYElemDoF = std::make_shared< p1_to_p2_div_1_blending_q6_ElementwiseOperator >( storage, minLevel, maxLevel, T );
      auto divXElem    = std::make_shared< P1ToP2ElementwiseDivxOperator >( storage, minLevel, maxLevel );
      auto divYElem    = std::make_shared< P1ToP2ElementwiseDivyOperator >( storage, minLevel, maxLevel );
      auto divXConst   = std::make_shared< P1ToP2ConstantDivxOperator >( storage, minLevel, maxLevel );
      auto divYConst   = std::make_shared< P1ToP2ConstantDivyOperator >( storage, minLevel, maxLevel );

      divXElemDoF->apply( u, resElementwiseDoFValue, maxLevel, All, Replace );
      divYElemDoF->apply( u, resElementwiseDoFValue, maxLevel, All, Add );
      divXElem->apply( u, resElementwise, maxLevel, All, Replace );
      divYElem->apply( u, resElementwise, maxLevel, All, Add );
      divXConst->apply( u, resConst, maxLevel, All, Replace );
      divYConst->apply( u, resConst, maxLevel, All, Add );

      real_t integValDoFValue = O.dotGlobal( resElementwiseDoFValue, maxLevel, All );
      real_t integValElem     = O.dotGlobal( resElementwise, maxLevel, All );
      real_t integValConst    = O.dotGlobal( resConst, maxLevel, All );
      real_t integValCtrl     = real_c( -1.841470984808 );

      WALBERLA_LOG_INFO( abs( integValDoFValue - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValDoFValue - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      auto divXElemDoFScaled = std::make_shared<
          ScaledOperator< p1_to_p2_div_0_blending_q6_ElementwiseOperator, P1Function< real_t >, P2Function< real_t > > >(
          storage, minLevel, maxLevel, divXElemDoF, real_c( 2 ) );
      auto divYElemDoFScaled = std::make_shared<
          ScaledOperator< p1_to_p2_div_1_blending_q6_ElementwiseOperator, P1Function< real_t >, P2Function< real_t > > >(
          storage, minLevel, maxLevel, divYElemDoF, real_c( 2 ) );

      auto divXElemScaled =
          std::make_shared< ScaledOperator< P1ToP2ElementwiseDivxOperator, P1Function< real_t >, P2Function< real_t > > >(
              storage, minLevel, maxLevel, divXElem, real_c( 2 ) );
      auto divYElemScaled =
          std::make_shared< ScaledOperator< P1ToP2ElementwiseDivyOperator, P1Function< real_t >, P2Function< real_t > > >(
              storage, minLevel, maxLevel, divYElem, real_c( 2 ) );

      auto divXConstScaled =
          std::make_shared< ScaledOperator< P1ToP2ConstantDivxOperator, P1Function< real_t >, P2Function< real_t > > >(
              storage, minLevel, maxLevel, divXConst, real_c( 2 ) );
      auto divYConstScaled =
          std::make_shared< ScaledOperator< P1ToP2ConstantDivyOperator, P1Function< real_t >, P2Function< real_t > > >(
              storage, minLevel, maxLevel, divYConst, real_c( 2 ) );

      divXElemDoFScaled->apply( u, resElementwiseDoFValue, maxLevel, All, Replace );
      divYElemDoFScaled->apply( u, resElementwiseDoFValue, maxLevel, All, Add );
      divXElemScaled->apply( u, resElementwise, maxLevel, All, Replace );
      divYElemScaled->apply( u, resElementwise, maxLevel, All, Add );
      divXConstScaled->apply( u, resConst, maxLevel, All, Replace );
      divYConstScaled->apply( u, resConst, maxLevel, All, Add );

      integValDoFValue = O.dotGlobal( resElementwiseDoFValue, maxLevel, All );
      integValElem     = O.dotGlobal( resElementwise, maxLevel, All );
      integValConst    = O.dotGlobal( resConst, maxLevel, All );
      integValCtrl     = real_c( -3.682941969616 );

      WALBERLA_LOG_INFO( abs( integValDoFValue - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValDoFValue - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      auto elementwiseDoFTypeX =
          hyteg::applyGEMV( *divXElemDoF, real_c( 2 ), u, real_c( 1 ), resElementwiseDoFValue, maxLevel, All );
      auto elementwiseDoFTypeY =
          hyteg::applyGEMV( *divYElemDoF, real_c( 2 ), u, real_c( 1 ), resElementwiseDoFValue, maxLevel, All );

      auto elementwiseTypeX = hyteg::applyGEMV( *divXElem, real_c( 2 ), u, real_c( 1 ), resElementwise, maxLevel, All );
      auto elementwiseTypeY = hyteg::applyGEMV( *divYElem, real_c( 2 ), u, real_c( 1 ), resElementwise, maxLevel, All );

      auto constantTypeX = hyteg::applyGEMV( *divXConst, real_c( 2 ), u, real_c( 1 ), resConst, maxLevel, All );
      auto constantTypeY = hyteg::applyGEMV( *divYConst, real_c( 2 ), u, real_c( 1 ), resConst, maxLevel, All );

      integValDoFValue = O.dotGlobal( resElementwiseDoFValue, maxLevel, All );
      integValElem     = O.dotGlobal( resElementwise, maxLevel, All );
      integValConst    = O.dotGlobal( resConst, maxLevel, All );
      integValCtrl     = real_c( -7.36588393923 );

      WALBERLA_LOG_INFO( abs( integValDoFValue - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValDoFValue - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled, smooth_jac_scaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseDoFTypeX );
      WALBERLA_CHECK_EQUAL( elementwiseDoFTypeX, GEMVType::GEMV );
      WALBERLA_LOG_INFO( elementwiseDoFTypeY );
      WALBERLA_CHECK_EQUAL( elementwiseDoFTypeY, GEMVType::GEMV );
      WALBERLA_LOG_INFO( elementwiseTypeX );
      WALBERLA_CHECK_EQUAL( elementwiseTypeX, GEMVType::GEMV );
      WALBERLA_LOG_INFO( elementwiseTypeY );
      WALBERLA_CHECK_EQUAL( elementwiseTypeY, GEMVType::GEMV );
      WALBERLA_LOG_INFO( constantTypeX );
      WALBERLA_CHECK_EQUAL( constantTypeX, GEMVType::MANUAL );
      WALBERLA_LOG_INFO( constantTypeY );
      WALBERLA_CHECK_EQUAL( constantTypeY, GEMVType::MANUAL );

      // create dummy matrix proxy, only checking divX here
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyDoFElem;
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyElem;
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyConst;

      dummyMatrixProxyDoFElem = std::make_shared< SumMatrixProxy >( real_c( 0 ) );
      dummyMatrixProxyElem    = std::make_shared< SumMatrixProxy >( real_c( 0 ) );
      dummyMatrixProxyConst   = std::make_shared< SumMatrixProxy >( real_c( 0 ) );

      P1Function< idx_t > enumeratorSrc( "enumeratorSrc", storage, minLevel, maxLevel );
      P2Function< idx_t > enumeratorDst( "enumeratorDst", storage, minLevel, maxLevel );
      enumeratorSrc.enumerate( maxLevel );
      enumeratorDst.enumerate( maxLevel );

      divXElemDoF->toMatrix( dummyMatrixProxyDoFElem, enumeratorSrc, enumeratorDst, maxLevel, All );
      divXElem->toMatrix( dummyMatrixProxyElem, enumeratorSrc, enumeratorDst, maxLevel, All );
      divXConst->toMatrix( dummyMatrixProxyConst, enumeratorSrc, enumeratorDst, maxLevel, All );

      real_t matrixSumValDoFElem = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal();
      real_t matrixSumValElem    = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal();
      real_t matrixSumValConst   = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal();

      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->setVal( real_c( 0 ) );

      divXElemDoFScaled->toMatrix( dummyMatrixProxyDoFElem, enumeratorSrc, enumeratorDst, maxLevel, All );
      divXElemScaled->toMatrix( dummyMatrixProxyElem, enumeratorSrc, enumeratorDst, maxLevel, All );
      divXConstScaled->toMatrix( dummyMatrixProxyConst, enumeratorSrc, enumeratorDst, maxLevel, All );

      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                              matrixSumValDoFElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                                matrixSumValDoFElem * real_c( 2 ) ),
                           toleranceMatrix );
      WALBERLA_LOG_INFO(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ),
          toleranceMatrix );
      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() -
                              matrixSumValConst * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() - matrixSumValConst * real_c( 2 ) ),
          toleranceMatrix );

      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->setVal( real_c( 0 ) );

      auto elementwiseDoFMatrixType = hyteg::applyToMatrixScaled(
          *divXElemDoF, real_c( 2 ), dummyMatrixProxyDoFElem, enumeratorSrc, enumeratorDst, maxLevel, All );
      auto elementwiseMatrixType =
          hyteg::applyToMatrixScaled( *divXElem, real_c( 2 ), dummyMatrixProxyElem, enumeratorSrc, enumeratorDst, maxLevel, All );
      auto constMatrixType = hyteg::applyToMatrixScaled(
          *divXConst, real_c( 2 ), dummyMatrixProxyConst, enumeratorSrc, enumeratorDst, maxLevel, All );

      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                              matrixSumValDoFElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                                matrixSumValDoFElem * real_c( 2 ) ),
                           toleranceMatrix );
      WALBERLA_LOG_INFO(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ),
          toleranceMatrix );
      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() -
                              matrixSumValConst * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() - matrixSumValConst * real_c( 2 ) ),
          toleranceMatrix );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled, smooth_jac_scaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseDoFMatrixType );
      WALBERLA_CHECK_EQUAL( elementwiseDoFMatrixType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( elementwiseMatrixType );
      WALBERLA_CHECK_EQUAL( elementwiseMatrixType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( constMatrixType );
      WALBERLA_CHECK_EQUAL( constMatrixType, GEMVType::MANUAL );
   }

   {
      // P1 To P2 3D Scaling Test

      const uint_t minLevel        = 2;
      const uint_t maxLevel        = 2;
      const real_t tolerance       = real_c( 5e-12 );
      const real_t toleranceMatrix = real_c( 1e-15 );

      // Init setup storage
      MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid(
          Point3D( real_c( 0 ), real_c( 0 ), real_c( 0 ) ), Point3D( real_c( 1 ), real_c( 1 ), real_c( 1 ) ), 3, 3, 3 );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

      // Create storage
      std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

      // Create functions
      P2Function< real_t > T( "T", storage, minLevel, maxLevel );
      P1Function< real_t > u( "u", storage, minLevel, maxLevel );
      P2Function< real_t > resConst( "resConst", storage, minLevel, maxLevel );
      P2Function< real_t > resElementwise( "resElementwise", storage, minLevel, maxLevel );
      P2Function< real_t > resElementwiseDoFValue( "resElementwiseDoFValue", storage, minLevel, maxLevel );
      P2Function< real_t > O( "O", storage, minLevel, maxLevel );

      std::function< real_t( const hyteg::Point3D& ) > uFct = []( const hyteg::Point3D& x ) {
         WALBERLA_UNUSED( x );
         return std::sin( x[0] ) + x[1] * x[1] + cos( x[2] );
      };

      u.interpolate( uFct, maxLevel, hyteg::All );
      O.interpolate( real_c( 1 ), maxLevel, hyteg::All );

      // operators
      std::function< real_t( const Point3D&, real_t, real_t, real_t, real_t, real_t ) > massScaling =
          [=]( const hyteg::Point3D& x, real_t temp, real_t tempCentroid, real_t centX, real_t centY, real_t centZ ) {
             return real_c( 1 );
          };

      auto divXElemDoF = std::make_shared< p1_to_p2_div_0_blending_q6_ElementwiseOperator >( storage, minLevel, maxLevel, T );
      auto divYElemDoF = std::make_shared< p1_to_p2_div_1_blending_q6_ElementwiseOperator >( storage, minLevel, maxLevel, T );
      auto divZElemDoF = std::make_shared< p1_to_p2_div_2_blending_q6_ElementwiseOperator >( storage, minLevel, maxLevel, T );
      auto divXElem    = std::make_shared< P1ToP2ElementwiseDivxOperator >( storage, minLevel, maxLevel );
      auto divYElem    = std::make_shared< P1ToP2ElementwiseDivyOperator >( storage, minLevel, maxLevel );
      auto divZElem    = std::make_shared< P1ToP2ElementwiseDivzOperator >( storage, minLevel, maxLevel );
      auto divXConst   = std::make_shared< P1ToP2ConstantDivxOperator >( storage, minLevel, maxLevel );
      auto divYConst   = std::make_shared< P1ToP2ConstantDivyOperator >( storage, minLevel, maxLevel );
      auto divZConst   = std::make_shared< P1ToP2ConstantDivzOperator >( storage, minLevel, maxLevel );

      divXElemDoF->apply( u, resElementwiseDoFValue, maxLevel, All, Replace );
      divYElemDoF->apply( u, resElementwiseDoFValue, maxLevel, All, Add );
      divZElemDoF->apply( u, resElementwiseDoFValue, maxLevel, All, Add );
      divXElem->apply( u, resElementwise, maxLevel, All, Replace );
      divYElem->apply( u, resElementwise, maxLevel, All, Add );
      divZElem->apply( u, resElementwise, maxLevel, All, Add );
      divXConst->apply( u, resConst, maxLevel, All, Replace );
      divYConst->apply( u, resConst, maxLevel, All, Add );
      divZConst->apply( u, resConst, maxLevel, All, Add );

      real_t integValDoFValue = O.dotGlobal( resElementwiseDoFValue, maxLevel, All );
      real_t integValElem     = O.dotGlobal( resElementwise, maxLevel, All );
      real_t integValConst    = O.dotGlobal( resConst, maxLevel, All );
      real_t integValCtrl     = real_c( -1.381773290675 );

      WALBERLA_LOG_INFO( abs( integValDoFValue - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValDoFValue - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      auto divXElemDoFScaled = std::make_shared<
          ScaledOperator< p1_to_p2_div_0_blending_q6_ElementwiseOperator, P1Function< real_t >, P2Function< real_t > > >(
          storage, minLevel, maxLevel, divXElemDoF, real_c( 2 ) );
      auto divYElemDoFScaled = std::make_shared<
          ScaledOperator< p1_to_p2_div_1_blending_q6_ElementwiseOperator, P1Function< real_t >, P2Function< real_t > > >(
          storage, minLevel, maxLevel, divYElemDoF, real_c( 2 ) );
      auto divZElemDoFScaled = std::make_shared<
          ScaledOperator< p1_to_p2_div_2_blending_q6_ElementwiseOperator, P1Function< real_t >, P2Function< real_t > > >(
          storage, minLevel, maxLevel, divZElemDoF, real_c( 2 ) );

      auto divXElemScaled =
          std::make_shared< ScaledOperator< P1ToP2ElementwiseDivxOperator, P1Function< real_t >, P2Function< real_t > > >(
              storage, minLevel, maxLevel, divXElem, real_c( 2 ) );
      auto divYElemScaled =
          std::make_shared< ScaledOperator< P1ToP2ElementwiseDivyOperator, P1Function< real_t >, P2Function< real_t > > >(
              storage, minLevel, maxLevel, divYElem, real_c( 2 ) );
      auto divZElemScaled =
          std::make_shared< ScaledOperator< P1ToP2ElementwiseDivzOperator, P1Function< real_t >, P2Function< real_t > > >(
              storage, minLevel, maxLevel, divZElem, real_c( 2 ) );

      auto divXConstScaled =
          std::make_shared< ScaledOperator< P1ToP2ConstantDivxOperator, P1Function< real_t >, P2Function< real_t > > >(
              storage, minLevel, maxLevel, divXConst, real_c( 2 ) );
      auto divYConstScaled =
          std::make_shared< ScaledOperator< P1ToP2ConstantDivyOperator, P1Function< real_t >, P2Function< real_t > > >(
              storage, minLevel, maxLevel, divYConst, real_c( 2 ) );
      auto divZConstScaled =
          std::make_shared< ScaledOperator< P1ToP2ConstantDivzOperator, P1Function< real_t >, P2Function< real_t > > >(
              storage, minLevel, maxLevel, divZConst, real_c( 2 ) );

      divXElemDoFScaled->apply( u, resElementwiseDoFValue, maxLevel, All, Replace );
      divYElemDoFScaled->apply( u, resElementwiseDoFValue, maxLevel, All, Add );
      divZElemDoFScaled->apply( u, resElementwiseDoFValue, maxLevel, All, Add );
      divXElemScaled->apply( u, resElementwise, maxLevel, All, Replace );
      divYElemScaled->apply( u, resElementwise, maxLevel, All, Add );
      divZElemScaled->apply( u, resElementwise, maxLevel, All, Add );
      divXConstScaled->apply( u, resConst, maxLevel, All, Replace );
      divYConstScaled->apply( u, resConst, maxLevel, All, Add );
      divZConstScaled->apply( u, resConst, maxLevel, All, Add );

      integValDoFValue = O.dotGlobal( resElementwiseDoFValue, maxLevel, All );
      integValElem     = O.dotGlobal( resElementwise, maxLevel, All );
      integValConst    = O.dotGlobal( resConst, maxLevel, All );
      integValCtrl     = real_c( -2.763546581350 );

      WALBERLA_LOG_INFO( abs( integValDoFValue - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValDoFValue - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      auto elementwiseDoFTypeX =
          hyteg::applyGEMV( *divXElemDoF, real_c( 2 ), u, real_c( 1 ), resElementwiseDoFValue, maxLevel, All );
      auto elementwiseDoFTypeY =
          hyteg::applyGEMV( *divYElemDoF, real_c( 2 ), u, real_c( 1 ), resElementwiseDoFValue, maxLevel, All );
      auto elementwiseDoFTypeZ =
          hyteg::applyGEMV( *divZElemDoF, real_c( 2 ), u, real_c( 1 ), resElementwiseDoFValue, maxLevel, All );

      auto elementwiseTypeX = hyteg::applyGEMV( *divXElem, real_c( 2 ), u, real_c( 1 ), resElementwise, maxLevel, All );
      auto elementwiseTypeY = hyteg::applyGEMV( *divYElem, real_c( 2 ), u, real_c( 1 ), resElementwise, maxLevel, All );
      auto elementwiseTypeZ = hyteg::applyGEMV( *divZElem, real_c( 2 ), u, real_c( 1 ), resElementwise, maxLevel, All );

      auto constantTypeX = hyteg::applyGEMV( *divXConst, real_c( 2 ), u, real_c( 1 ), resConst, maxLevel, All );
      auto constantTypeY = hyteg::applyGEMV( *divYConst, real_c( 2 ), u, real_c( 1 ), resConst, maxLevel, All );
      auto constantTypeZ = hyteg::applyGEMV( *divZConst, real_c( 2 ), u, real_c( 1 ), resConst, maxLevel, All );

      integValDoFValue = O.dotGlobal( resElementwiseDoFValue, maxLevel, All );
      integValElem     = O.dotGlobal( resElementwise, maxLevel, All );
      integValConst    = O.dotGlobal( resConst, maxLevel, All );
      integValCtrl     = real_c( -5.5270931627 );

      WALBERLA_LOG_INFO( abs( integValDoFValue - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValDoFValue - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled, smooth_jac_scaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseDoFTypeX );
      WALBERLA_CHECK_EQUAL( elementwiseDoFTypeX, GEMVType::GEMV );
      WALBERLA_LOG_INFO( elementwiseDoFTypeY );
      WALBERLA_CHECK_EQUAL( elementwiseDoFTypeY, GEMVType::GEMV );
      WALBERLA_LOG_INFO( elementwiseDoFTypeZ );
      WALBERLA_CHECK_EQUAL( elementwiseDoFTypeZ, GEMVType::GEMV );
      WALBERLA_LOG_INFO( elementwiseTypeX );
      WALBERLA_CHECK_EQUAL( elementwiseTypeX, GEMVType::GEMV );
      WALBERLA_LOG_INFO( elementwiseTypeY );
      WALBERLA_CHECK_EQUAL( elementwiseTypeY, GEMVType::GEMV );
      WALBERLA_LOG_INFO( elementwiseTypeZ );
      WALBERLA_CHECK_EQUAL( elementwiseTypeZ, GEMVType::GEMV );
      WALBERLA_LOG_INFO( constantTypeX );
      WALBERLA_CHECK_EQUAL( constantTypeX, GEMVType::MANUAL );
      WALBERLA_LOG_INFO( constantTypeY );
      WALBERLA_CHECK_EQUAL( constantTypeY, GEMVType::MANUAL );
      WALBERLA_LOG_INFO( constantTypeZ );
      WALBERLA_CHECK_EQUAL( constantTypeZ, GEMVType::MANUAL );

      // create dummy matrix proxy, only checking divX here
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyDoFElem;
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyElem;
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyConst;

      dummyMatrixProxyDoFElem = std::make_shared< SumMatrixProxy >( real_c( 0 ) );
      dummyMatrixProxyElem    = std::make_shared< SumMatrixProxy >( real_c( 0 ) );
      dummyMatrixProxyConst   = std::make_shared< SumMatrixProxy >( real_c( 0 ) );

      P1Function< idx_t > enumeratorSrc( "enumeratorSrc", storage, minLevel, maxLevel );
      P2Function< idx_t > enumeratorDst( "enumeratorDst", storage, minLevel, maxLevel );
      enumeratorSrc.enumerate( maxLevel );
      enumeratorDst.enumerate( maxLevel );

      divXElemDoF->toMatrix( dummyMatrixProxyDoFElem, enumeratorSrc, enumeratorDst, maxLevel, All );
      divXElem->toMatrix( dummyMatrixProxyElem, enumeratorSrc, enumeratorDst, maxLevel, All );
      divXConst->toMatrix( dummyMatrixProxyConst, enumeratorSrc, enumeratorDst, maxLevel, All );

      real_t matrixSumValDoFElem = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal();
      real_t matrixSumValElem    = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal();
      real_t matrixSumValConst   = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal();

      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->setVal( real_c( 0 ) );

      divXElemDoFScaled->toMatrix( dummyMatrixProxyDoFElem, enumeratorSrc, enumeratorDst, maxLevel, All );
      divXElemScaled->toMatrix( dummyMatrixProxyElem, enumeratorSrc, enumeratorDst, maxLevel, All );
      divXConstScaled->toMatrix( dummyMatrixProxyConst, enumeratorSrc, enumeratorDst, maxLevel, All );

      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                              matrixSumValDoFElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                                matrixSumValDoFElem * real_c( 2 ) ),
                           toleranceMatrix );
      WALBERLA_LOG_INFO(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ),
          toleranceMatrix );
      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() -
                              matrixSumValConst * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() - matrixSumValConst * real_c( 2 ) ),
          toleranceMatrix );

      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->setVal( real_c( 0 ) );

      auto elementwiseDoFMatrixType = hyteg::applyToMatrixScaled(
          *divXElemDoF, real_c( 2 ), dummyMatrixProxyDoFElem, enumeratorSrc, enumeratorDst, maxLevel, All );
      auto elementwiseMatrixType =
          hyteg::applyToMatrixScaled( *divXElem, real_c( 2 ), dummyMatrixProxyElem, enumeratorSrc, enumeratorDst, maxLevel, All );
      auto constMatrixType = hyteg::applyToMatrixScaled(
          *divXConst, real_c( 2 ), dummyMatrixProxyConst, enumeratorSrc, enumeratorDst, maxLevel, All );

      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                              matrixSumValDoFElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                                matrixSumValDoFElem * real_c( 2 ) ),
                           toleranceMatrix );
      WALBERLA_LOG_INFO(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ),
          toleranceMatrix );
      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() -
                              matrixSumValConst * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() - matrixSumValConst * real_c( 2 ) ),
          toleranceMatrix );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled, smooth_jac_scaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseDoFMatrixType );
      WALBERLA_CHECK_EQUAL( elementwiseDoFMatrixType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( elementwiseMatrixType );
      WALBERLA_CHECK_EQUAL( elementwiseMatrixType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( constMatrixType );
      WALBERLA_CHECK_EQUAL( constMatrixType, GEMVType::MANUAL );
   }

   {
      // P2 To P1 2D Scaling Test

      const uint_t minLevel        = 3;
      const uint_t maxLevel        = 3;
      const real_t tolerance       = real_c( 2e-12 );
      const real_t toleranceMatrix = real_c( 1e-15 );

      // Init setup storage
      MeshInfo meshInfo = MeshInfo::meshRectangle(
          Point2D( real_c( 0 ), real_c( 0 ) ), Point2D( real_c( 1 ), real_c( 1 ) ), MeshInfo::CROSS, 3, 3 );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

      // Create storage
      std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

      // Create functions
      P2Function< real_t > T( "T", storage, minLevel, maxLevel );
      P2Function< real_t > u( "u", storage, minLevel, maxLevel );
      P1Function< real_t > resConst( "resConst", storage, minLevel, maxLevel );
      P1Function< real_t > resElementwise( "resElementwise", storage, minLevel, maxLevel );
      P1Function< real_t > resElementwiseDoFValue( "resElementwiseDoFValue", storage, minLevel, maxLevel );
      P1Function< real_t > O( "O", storage, minLevel, maxLevel );

      std::function< real_t( const hyteg::Point3D& ) > uFct = []( const hyteg::Point3D& x ) {
         WALBERLA_UNUSED( x );
         return std::sin( x[0] ) + x[1] * x[1];
      };

      u.interpolate( uFct, maxLevel, hyteg::All );
      O.interpolate( real_c( 1 ), maxLevel, hyteg::All );

      // operators
      std::function< real_t( const Point3D&, real_t, real_t, real_t, real_t, real_t ) > massScaling =
          [=]( const hyteg::Point3D& x, real_t temp, real_t tempCentroid, real_t centX, real_t centY, real_t centZ ) {
             return real_c( 1 );
          };

      auto divXElemDoF = std::make_shared< p2_to_p1_div_0_blending_q6_ElementwiseOperator >( storage, minLevel, maxLevel, T );
      auto divYElemDoF = std::make_shared< p2_to_p1_div_1_blending_q6_ElementwiseOperator >( storage, minLevel, maxLevel, T );
      auto divXElem    = std::make_shared< P2ToP1ElementwiseDivxOperator >( storage, minLevel, maxLevel );
      auto divYElem    = std::make_shared< P2ToP1ElementwiseDivyOperator >( storage, minLevel, maxLevel );
      auto divXConst   = std::make_shared< P2ToP1ConstantDivxOperator >( storage, minLevel, maxLevel );
      auto divYConst   = std::make_shared< P2ToP1ConstantDivyOperator >( storage, minLevel, maxLevel );

      divXElemDoF->apply( u, resElementwiseDoFValue, maxLevel, All, Replace );
      divYElemDoF->apply( u, resElementwiseDoFValue, maxLevel, All, Add );
      divXElem->apply( u, resElementwise, maxLevel, All, Replace );
      divYElem->apply( u, resElementwise, maxLevel, All, Add );
      divXConst->apply( u, resConst, maxLevel, All, Replace );
      divYConst->apply( u, resConst, maxLevel, All, Add );

      real_t integValDoFValue = O.dotGlobal( resElementwiseDoFValue, maxLevel, All );
      real_t integValElem     = O.dotGlobal( resElementwise, maxLevel, All );
      real_t integValConst    = O.dotGlobal( resConst, maxLevel, All );
      real_t integValCtrl     = real_c( -1.841470984808 );

      WALBERLA_LOG_INFO( abs( integValDoFValue - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValDoFValue - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      auto divXElemDoFScaled = std::make_shared<
          ScaledOperator< p2_to_p1_div_0_blending_q6_ElementwiseOperator, P2Function< real_t >, P1Function< real_t > > >(
          storage, minLevel, maxLevel, divXElemDoF, real_c( 2 ) );
      auto divYElemDoFScaled = std::make_shared<
          ScaledOperator< p2_to_p1_div_1_blending_q6_ElementwiseOperator, P2Function< real_t >, P1Function< real_t > > >(
          storage, minLevel, maxLevel, divYElemDoF, real_c( 2 ) );

      auto divXElemScaled =
          std::make_shared< ScaledOperator< P2ToP1ElementwiseDivxOperator, P2Function< real_t >, P1Function< real_t > > >(
              storage, minLevel, maxLevel, divXElem, real_c( 2 ) );
      auto divYElemScaled =
          std::make_shared< ScaledOperator< P2ToP1ElementwiseDivyOperator, P2Function< real_t >, P1Function< real_t > > >(
              storage, minLevel, maxLevel, divYElem, real_c( 2 ) );

      auto divXConstScaled =
          std::make_shared< ScaledOperator< P2ToP1ConstantDivxOperator, P2Function< real_t >, P1Function< real_t > > >(
              storage, minLevel, maxLevel, divXConst, real_c( 2 ) );
      auto divYConstScaled =
          std::make_shared< ScaledOperator< P2ToP1ConstantDivyOperator, P2Function< real_t >, P1Function< real_t > > >(
              storage, minLevel, maxLevel, divYConst, real_c( 2 ) );

      divXElemDoFScaled->apply( u, resElementwiseDoFValue, maxLevel, All, Replace );
      divYElemDoFScaled->apply( u, resElementwiseDoFValue, maxLevel, All, Add );
      divXElemScaled->apply( u, resElementwise, maxLevel, All, Replace );
      divYElemScaled->apply( u, resElementwise, maxLevel, All, Add );
      divXConstScaled->apply( u, resConst, maxLevel, All, Replace );
      divYConstScaled->apply( u, resConst, maxLevel, All, Add );

      integValDoFValue = O.dotGlobal( resElementwiseDoFValue, maxLevel, All );
      integValElem     = O.dotGlobal( resElementwise, maxLevel, All );
      integValConst    = O.dotGlobal( resConst, maxLevel, All );
      integValCtrl     = real_c( -3.682941969616 );

      WALBERLA_LOG_INFO( abs( integValDoFValue - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValDoFValue - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      auto elementwiseDoFTypeX =
          hyteg::applyGEMV( *divXElemDoF, real_c( 2 ), u, real_c( 1 ), resElementwiseDoFValue, maxLevel, All );
      auto elementwiseDoFTypeY =
          hyteg::applyGEMV( *divYElemDoF, real_c( 2 ), u, real_c( 1 ), resElementwiseDoFValue, maxLevel, All );

      auto elementwiseTypeX = hyteg::applyGEMV( *divXElem, real_c( 2 ), u, real_c( 1 ), resElementwise, maxLevel, All );
      auto elementwiseTypeY = hyteg::applyGEMV( *divYElem, real_c( 2 ), u, real_c( 1 ), resElementwise, maxLevel, All );

      auto constantTypeX = hyteg::applyGEMV( *divXConst, real_c( 2 ), u, real_c( 1 ), resConst, maxLevel, All );
      auto constantTypeY = hyteg::applyGEMV( *divYConst, real_c( 2 ), u, real_c( 1 ), resConst, maxLevel, All );

      integValDoFValue = O.dotGlobal( resElementwiseDoFValue, maxLevel, All );
      integValElem     = O.dotGlobal( resElementwise, maxLevel, All );
      integValConst    = O.dotGlobal( resConst, maxLevel, All );
      integValCtrl     = real_c( -7.36588393923 );

      WALBERLA_LOG_INFO( abs( integValDoFValue - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValDoFValue - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled, smooth_jac_scaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseDoFTypeX );
      WALBERLA_CHECK_EQUAL( elementwiseDoFTypeX, GEMVType::GEMV );
      WALBERLA_LOG_INFO( elementwiseDoFTypeY );
      WALBERLA_CHECK_EQUAL( elementwiseDoFTypeY, GEMVType::GEMV );
      WALBERLA_LOG_INFO( elementwiseTypeX );
      WALBERLA_CHECK_EQUAL( elementwiseTypeX, GEMVType::GEMV );
      WALBERLA_LOG_INFO( elementwiseTypeY );
      WALBERLA_CHECK_EQUAL( elementwiseTypeY, GEMVType::GEMV );
      WALBERLA_LOG_INFO( constantTypeX );
      WALBERLA_CHECK_EQUAL( constantTypeX, GEMVType::MANUAL );
      WALBERLA_LOG_INFO( constantTypeY );
      WALBERLA_CHECK_EQUAL( constantTypeY, GEMVType::MANUAL );

      // create dummy matrix proxy, only checking divX here
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyDoFElem;
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyElem;
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyConst;

      dummyMatrixProxyDoFElem = std::make_shared< SumMatrixProxy >( real_c( 0 ) );
      dummyMatrixProxyElem    = std::make_shared< SumMatrixProxy >( real_c( 0 ) );
      dummyMatrixProxyConst   = std::make_shared< SumMatrixProxy >( real_c( 0 ) );

      P2Function< idx_t > enumeratorSrc( "enumeratorSrc", storage, minLevel, maxLevel );
      P1Function< idx_t > enumeratorDst( "enumeratorDst", storage, minLevel, maxLevel );
      enumeratorSrc.enumerate( maxLevel );
      enumeratorDst.enumerate( maxLevel );

      divXElemDoF->toMatrix( dummyMatrixProxyDoFElem, enumeratorSrc, enumeratorDst, maxLevel, All );
      divXElem->toMatrix( dummyMatrixProxyElem, enumeratorSrc, enumeratorDst, maxLevel, All );
      divXConst->toMatrix( dummyMatrixProxyConst, enumeratorSrc, enumeratorDst, maxLevel, All );

      real_t matrixSumValDoFElem = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal();
      real_t matrixSumValElem    = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal();
      real_t matrixSumValConst   = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal();

      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->setVal( real_c( 0 ) );

      divXElemDoFScaled->toMatrix( dummyMatrixProxyDoFElem, enumeratorSrc, enumeratorDst, maxLevel, All );
      divXElemScaled->toMatrix( dummyMatrixProxyElem, enumeratorSrc, enumeratorDst, maxLevel, All );
      divXConstScaled->toMatrix( dummyMatrixProxyConst, enumeratorSrc, enumeratorDst, maxLevel, All );

      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                              matrixSumValDoFElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                                matrixSumValDoFElem * real_c( 2 ) ),
                           toleranceMatrix );
      WALBERLA_LOG_INFO(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ),
          toleranceMatrix );
      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() -
                              matrixSumValConst * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() - matrixSumValConst * real_c( 2 ) ),
          toleranceMatrix );

      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->setVal( real_c( 0 ) );

      auto elementwiseDoFMatrixType = hyteg::applyToMatrixScaled(
          *divXElemDoF, real_c( 2 ), dummyMatrixProxyDoFElem, enumeratorSrc, enumeratorDst, maxLevel, All );
      auto elementwiseMatrixType =
          hyteg::applyToMatrixScaled( *divXElem, real_c( 2 ), dummyMatrixProxyElem, enumeratorSrc, enumeratorDst, maxLevel, All );
      auto constMatrixType = hyteg::applyToMatrixScaled(
          *divXConst, real_c( 2 ), dummyMatrixProxyConst, enumeratorSrc, enumeratorDst, maxLevel, All );

      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                              matrixSumValDoFElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                                matrixSumValDoFElem * real_c( 2 ) ),
                           toleranceMatrix );
      WALBERLA_LOG_INFO(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ),
          toleranceMatrix );
      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() -
                              matrixSumValConst * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() - matrixSumValConst * real_c( 2 ) ),
          toleranceMatrix );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled, smooth_jac_scaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseDoFMatrixType );
      WALBERLA_CHECK_EQUAL( elementwiseDoFMatrixType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( elementwiseMatrixType );
      WALBERLA_CHECK_EQUAL( elementwiseMatrixType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( constMatrixType );
      WALBERLA_CHECK_EQUAL( constMatrixType, GEMVType::MANUAL );
   }

   {
      // P2 To P1 3D Scaling Test

      const uint_t minLevel        = 2;
      const uint_t maxLevel        = 2;
      const real_t tolerance       = real_c( 5e-12 );
      const real_t toleranceMatrix = real_c( 1e-15 );

      // Init setup storage
      MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid(
          Point3D( real_c( 0 ), real_c( 0 ), real_c( 0 ) ), Point3D( real_c( 1 ), real_c( 1 ), real_c( 1 ) ), 3, 3, 3 );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

      // Create storage
      std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

      // Create functions
      P2Function< real_t > T( "T", storage, minLevel, maxLevel );
      P2Function< real_t > u( "u", storage, minLevel, maxLevel );
      P1Function< real_t > resConst( "resConst", storage, minLevel, maxLevel );
      P1Function< real_t > resElementwise( "resElementwise", storage, minLevel, maxLevel );
      P1Function< real_t > resElementwiseDoFValue( "resElementwiseDoFValue", storage, minLevel, maxLevel );
      P1Function< real_t > O( "O", storage, minLevel, maxLevel );

      std::function< real_t( const hyteg::Point3D& ) > uFct = []( const hyteg::Point3D& x ) {
         WALBERLA_UNUSED( x );
         return std::sin( x[0] ) + x[1] * x[1] + cos( x[2] );
      };

      u.interpolate( uFct, maxLevel, hyteg::All );
      O.interpolate( real_c( 1 ), maxLevel, hyteg::All );

      // operators
      std::function< real_t( const Point3D&, real_t, real_t, real_t, real_t, real_t ) > massScaling =
          [=]( const hyteg::Point3D& x, real_t temp, real_t tempCentroid, real_t centX, real_t centY, real_t centZ ) {
             return real_c( 1 );
          };

      auto divXElemDoF = std::make_shared< p2_to_p1_div_0_blending_q6_ElementwiseOperator >( storage, minLevel, maxLevel, T );
      auto divYElemDoF = std::make_shared< p2_to_p1_div_1_blending_q6_ElementwiseOperator >( storage, minLevel, maxLevel, T );
      auto divZElemDoF = std::make_shared< p2_to_p1_div_2_blending_q6_ElementwiseOperator >( storage, minLevel, maxLevel, T );
      auto divXElem    = std::make_shared< P2ToP1ElementwiseDivxOperator >( storage, minLevel, maxLevel );
      auto divYElem    = std::make_shared< P2ToP1ElementwiseDivyOperator >( storage, minLevel, maxLevel );
      auto divZElem    = std::make_shared< P2ToP1ElementwiseDivzOperator >( storage, minLevel, maxLevel );
      auto divXConst   = std::make_shared< P2ToP1ConstantDivxOperator >( storage, minLevel, maxLevel );
      auto divYConst   = std::make_shared< P2ToP1ConstantDivyOperator >( storage, minLevel, maxLevel );
      auto divZConst   = std::make_shared< P2ToP1ConstantDivzOperator >( storage, minLevel, maxLevel );

      divXElemDoF->apply( u, resElementwiseDoFValue, maxLevel, All, Replace );
      divYElemDoF->apply( u, resElementwiseDoFValue, maxLevel, All, Add );
      divZElemDoF->apply( u, resElementwiseDoFValue, maxLevel, All, Add );
      divXElem->apply( u, resElementwise, maxLevel, All, Replace );
      divYElem->apply( u, resElementwise, maxLevel, All, Add );
      divZElem->apply( u, resElementwise, maxLevel, All, Add );
      divXConst->apply( u, resConst, maxLevel, All, Replace );
      divYConst->apply( u, resConst, maxLevel, All, Add );
      divZConst->apply( u, resConst, maxLevel, All, Add );

      real_t integValDoFValue = O.dotGlobal( resElementwiseDoFValue, maxLevel, All );
      real_t integValElem     = O.dotGlobal( resElementwise, maxLevel, All );
      real_t integValConst    = O.dotGlobal( resConst, maxLevel, All );
      real_t integValCtrl     = real_c( -1.381773290675 );

      WALBERLA_LOG_INFO( abs( integValDoFValue - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValDoFValue - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      auto divXElemDoFScaled = std::make_shared<
          ScaledOperator< p2_to_p1_div_0_blending_q6_ElementwiseOperator, P2Function< real_t >, P1Function< real_t > > >(
          storage, minLevel, maxLevel, divXElemDoF, real_c( 2 ) );
      auto divYElemDoFScaled = std::make_shared<
          ScaledOperator< p2_to_p1_div_1_blending_q6_ElementwiseOperator, P2Function< real_t >, P1Function< real_t > > >(
          storage, minLevel, maxLevel, divYElemDoF, real_c( 2 ) );
      auto divZElemDoFScaled = std::make_shared<
          ScaledOperator< p2_to_p1_div_2_blending_q6_ElementwiseOperator, P2Function< real_t >, P1Function< real_t > > >(
          storage, minLevel, maxLevel, divZElemDoF, real_c( 2 ) );

      auto divXElemScaled =
          std::make_shared< ScaledOperator< P2ToP1ElementwiseDivxOperator, P2Function< real_t >, P1Function< real_t > > >(
              storage, minLevel, maxLevel, divXElem, real_c( 2 ) );
      auto divYElemScaled =
          std::make_shared< ScaledOperator< P2ToP1ElementwiseDivyOperator, P2Function< real_t >, P1Function< real_t > > >(
              storage, minLevel, maxLevel, divYElem, real_c( 2 ) );
      auto divZElemScaled =
          std::make_shared< ScaledOperator< P2ToP1ElementwiseDivzOperator, P2Function< real_t >, P1Function< real_t > > >(
              storage, minLevel, maxLevel, divZElem, real_c( 2 ) );

      auto divXConstScaled =
          std::make_shared< ScaledOperator< P2ToP1ConstantDivxOperator, P2Function< real_t >, P1Function< real_t > > >(
              storage, minLevel, maxLevel, divXConst, real_c( 2 ) );
      auto divYConstScaled =
          std::make_shared< ScaledOperator< P2ToP1ConstantDivyOperator, P2Function< real_t >, P1Function< real_t > > >(
              storage, minLevel, maxLevel, divYConst, real_c( 2 ) );
      auto divZConstScaled =
          std::make_shared< ScaledOperator< P2ToP1ConstantDivzOperator, P2Function< real_t >, P1Function< real_t > > >(
              storage, minLevel, maxLevel, divZConst, real_c( 2 ) );

      divXElemDoFScaled->apply( u, resElementwiseDoFValue, maxLevel, All, Replace );
      divYElemDoFScaled->apply( u, resElementwiseDoFValue, maxLevel, All, Add );
      divZElemDoFScaled->apply( u, resElementwiseDoFValue, maxLevel, All, Add );
      divXElemScaled->apply( u, resElementwise, maxLevel, All, Replace );
      divYElemScaled->apply( u, resElementwise, maxLevel, All, Add );
      divZElemScaled->apply( u, resElementwise, maxLevel, All, Add );
      divXConstScaled->apply( u, resConst, maxLevel, All, Replace );
      divYConstScaled->apply( u, resConst, maxLevel, All, Add );
      divZConstScaled->apply( u, resConst, maxLevel, All, Add );

      integValDoFValue = O.dotGlobal( resElementwiseDoFValue, maxLevel, All );
      integValElem     = O.dotGlobal( resElementwise, maxLevel, All );
      integValConst    = O.dotGlobal( resConst, maxLevel, All );
      integValCtrl     = real_c( -2.763546581350 );

      WALBERLA_LOG_INFO( abs( integValDoFValue - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValDoFValue - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      auto elementwiseDoFTypeX =
          hyteg::applyGEMV( *divXElemDoF, real_c( 2 ), u, real_c( 1 ), resElementwiseDoFValue, maxLevel, All );
      auto elementwiseDoFTypeY =
          hyteg::applyGEMV( *divYElemDoF, real_c( 2 ), u, real_c( 1 ), resElementwiseDoFValue, maxLevel, All );
      auto elementwiseDoFTypeZ =
          hyteg::applyGEMV( *divZElemDoF, real_c( 2 ), u, real_c( 1 ), resElementwiseDoFValue, maxLevel, All );

      auto elementwiseTypeX = hyteg::applyGEMV( *divXElem, real_c( 2 ), u, real_c( 1 ), resElementwise, maxLevel, All );
      auto elementwiseTypeY = hyteg::applyGEMV( *divYElem, real_c( 2 ), u, real_c( 1 ), resElementwise, maxLevel, All );
      auto elementwiseTypeZ = hyteg::applyGEMV( *divZElem, real_c( 2 ), u, real_c( 1 ), resElementwise, maxLevel, All );

      auto constantTypeX = hyteg::applyGEMV( *divXConst, real_c( 2 ), u, real_c( 1 ), resConst, maxLevel, All );
      auto constantTypeY = hyteg::applyGEMV( *divYConst, real_c( 2 ), u, real_c( 1 ), resConst, maxLevel, All );
      auto constantTypeZ = hyteg::applyGEMV( *divZConst, real_c( 2 ), u, real_c( 1 ), resConst, maxLevel, All );

      integValDoFValue = O.dotGlobal( resElementwiseDoFValue, maxLevel, All );
      integValElem     = O.dotGlobal( resElementwise, maxLevel, All );
      integValConst    = O.dotGlobal( resConst, maxLevel, All );
      integValCtrl     = real_c( -5.5270931627 );

      WALBERLA_LOG_INFO( abs( integValDoFValue - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValDoFValue - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled, smooth_jac_scaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseDoFTypeX );
      WALBERLA_CHECK_EQUAL( elementwiseDoFTypeX, GEMVType::GEMV );
      WALBERLA_LOG_INFO( elementwiseDoFTypeY );
      WALBERLA_CHECK_EQUAL( elementwiseDoFTypeY, GEMVType::GEMV );
      WALBERLA_LOG_INFO( elementwiseDoFTypeZ );
      WALBERLA_CHECK_EQUAL( elementwiseDoFTypeZ, GEMVType::GEMV );
      WALBERLA_LOG_INFO( elementwiseTypeX );
      WALBERLA_CHECK_EQUAL( elementwiseTypeX, GEMVType::GEMV );
      WALBERLA_LOG_INFO( elementwiseTypeY );
      WALBERLA_CHECK_EQUAL( elementwiseTypeY, GEMVType::GEMV );
      WALBERLA_LOG_INFO( elementwiseTypeZ );
      WALBERLA_CHECK_EQUAL( elementwiseTypeZ, GEMVType::GEMV );
      WALBERLA_LOG_INFO( constantTypeX );
      WALBERLA_CHECK_EQUAL( constantTypeX, GEMVType::MANUAL );
      WALBERLA_LOG_INFO( constantTypeY );
      WALBERLA_CHECK_EQUAL( constantTypeY, GEMVType::MANUAL );
      WALBERLA_LOG_INFO( constantTypeZ );
      WALBERLA_CHECK_EQUAL( constantTypeZ, GEMVType::MANUAL );

      // create dummy matrix proxy, only checking divX here
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyDoFElem;
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyElem;
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyConst;

      dummyMatrixProxyDoFElem = std::make_shared< SumMatrixProxy >( real_c( 0 ) );
      dummyMatrixProxyElem    = std::make_shared< SumMatrixProxy >( real_c( 0 ) );
      dummyMatrixProxyConst   = std::make_shared< SumMatrixProxy >( real_c( 0 ) );

      P2Function< idx_t > enumeratorSrc( "enumeratorSrc", storage, minLevel, maxLevel );
      P1Function< idx_t > enumeratorDst( "enumeratorDst", storage, minLevel, maxLevel );
      enumeratorSrc.enumerate( maxLevel );
      enumeratorDst.enumerate( maxLevel );

      divXElemDoF->toMatrix( dummyMatrixProxyDoFElem, enumeratorSrc, enumeratorDst, maxLevel, All );
      divXElem->toMatrix( dummyMatrixProxyElem, enumeratorSrc, enumeratorDst, maxLevel, All );
      divXConst->toMatrix( dummyMatrixProxyConst, enumeratorSrc, enumeratorDst, maxLevel, All );

      real_t matrixSumValDoFElem = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal();
      real_t matrixSumValElem    = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal();
      real_t matrixSumValConst   = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal();

      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->setVal( real_c( 0 ) );

      divXElemDoFScaled->toMatrix( dummyMatrixProxyDoFElem, enumeratorSrc, enumeratorDst, maxLevel, All );
      divXElemScaled->toMatrix( dummyMatrixProxyElem, enumeratorSrc, enumeratorDst, maxLevel, All );
      divXConstScaled->toMatrix( dummyMatrixProxyConst, enumeratorSrc, enumeratorDst, maxLevel, All );

      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                              matrixSumValDoFElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                                matrixSumValDoFElem * real_c( 2 ) ),
                           toleranceMatrix );
      WALBERLA_LOG_INFO(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ),
          toleranceMatrix );
      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() -
                              matrixSumValConst * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() - matrixSumValConst * real_c( 2 ) ),
          toleranceMatrix );

      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->setVal( real_c( 0 ) );
      std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->setVal( real_c( 0 ) );

      auto elementwiseDoFMatrixType = hyteg::applyToMatrixScaled(
          *divXElemDoF, real_c( 2 ), dummyMatrixProxyDoFElem, enumeratorSrc, enumeratorDst, maxLevel, All );
      auto elementwiseMatrixType =
          hyteg::applyToMatrixScaled( *divXElem, real_c( 2 ), dummyMatrixProxyElem, enumeratorSrc, enumeratorDst, maxLevel, All );
      auto constMatrixType = hyteg::applyToMatrixScaled(
          *divXConst, real_c( 2 ), dummyMatrixProxyConst, enumeratorSrc, enumeratorDst, maxLevel, All );

      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                              matrixSumValDoFElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyDoFElem )->getVal() -
                                matrixSumValDoFElem * real_c( 2 ) ),
                           toleranceMatrix );
      WALBERLA_LOG_INFO(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal() - matrixSumValElem * real_c( 2 ) ),
          toleranceMatrix );
      WALBERLA_LOG_INFO( abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() -
                              matrixSumValConst * real_c( 2 ) ) );
      WALBERLA_CHECK_LESS(
          abs( std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal() - matrixSumValConst * real_c( 2 ) ),
          toleranceMatrix );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled, smooth_jac_scaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseDoFMatrixType );
      WALBERLA_CHECK_EQUAL( elementwiseDoFMatrixType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( elementwiseMatrixType );
      WALBERLA_CHECK_EQUAL( elementwiseMatrixType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( constMatrixType );
      WALBERLA_CHECK_EQUAL( constMatrixType, GEMVType::MANUAL );
   }

   {
      // P1 2D smooth_jac_scaled Test

      const uint_t minLevel  = 3;
      const uint_t maxLevel  = 3;
      const real_t tolerance = real_c( 1e-15 );
      const real_t omega     = 0.5;

      // Init setup storage
      MeshInfo meshInfo = MeshInfo::meshRectangle(
          Point2D( real_c( 0 ), real_c( 0 ) ), Point2D( real_c( 1 ), real_c( 1 ) ), MeshInfo::CROSS, 3, 3 );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

      // Create storage
      std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

      // Create functions
      P2Function< real_t > T( "T", storage, minLevel, maxLevel );
      P1Function< real_t > u( "u", storage, minLevel, maxLevel );

      P1Function< real_t > rhsConst( "rhsConst", storage, minLevel, maxLevel );
      P1Function< real_t > rhsElementwise( "rhsElementwise", storage, minLevel, maxLevel );
      P1Function< real_t > rhsElementwiseDoFValue( "rhsElementwiseDoFValue", storage, minLevel, maxLevel );

      P1Function< real_t > manualInvDiagConst( "manualInvDiagConst", storage, minLevel, maxLevel );
      P1Function< real_t > manualInvDiagElementwise( "manualInvDiagElementwise", storage, minLevel, maxLevel );
      P1Function< real_t > manualInvDiagElementwiseDoFValue( "manualInvDiagElementwiseDoFValue", storage, minLevel, maxLevel );

      P1Function< real_t > resConst( "resConst", storage, minLevel, maxLevel );
      P1Function< real_t > resElementwise( "resElementwise", storage, minLevel, maxLevel );
      P1Function< real_t > resElementwiseDoFValue( "resElementwiseDoFValue", storage, minLevel, maxLevel );

      P1Function< real_t > resConstGEMV( "resConstGEMV", storage, minLevel, maxLevel );
      P1Function< real_t > resElementwiseGEMV( "resElementwiseGEMV", storage, minLevel, maxLevel );
      P1Function< real_t > resElementwiseDoFValueGEMV( "resElementwiseDoFValueGEMV", storage, minLevel, maxLevel );

      std::function< real_t( const hyteg::Point3D& ) > uFct = []( const hyteg::Point3D& x ) {
         WALBERLA_UNUSED( x );
         return std::sin( x[0] ) + x[1] * x[1];
      };

      u.interpolate( uFct, maxLevel, hyteg::All );

      // operators
      std::function< real_t( const Point3D&, real_t, real_t, real_t, real_t, real_t ) > massScaling =
          [=]( const hyteg::Point3D& x, real_t temp, real_t tempCentroid, real_t centX, real_t centY, real_t centZ ) {
             return real_c( 1 );
          };

      auto CentroidMass = std::make_shared< p1_k_mass_visc_centroid_blending_q4_ElementwiseOperator >(
          storage, minLevel, maxLevel, T, massScaling );
      auto ElementwiseMass = std::make_shared< P1ElementwiseMassOperator >( storage, minLevel, maxLevel );
      auto ConstantMass    = std::make_shared< P1ConstantMassOperator >( storage, minLevel, maxLevel );

      // scaled operators
      auto CentroidMassScaled = std::make_shared<
          ScaledOperator< p1_k_mass_visc_centroid_blending_q4_ElementwiseOperator, P1Function< real_t >, P1Function< real_t > > >(
          storage, minLevel, maxLevel, CentroidMass, real_c( 2 ) );
      auto ElementwiseMassScaled =
          std::make_shared< ScaledOperator< P1ElementwiseMassOperator, P1Function< real_t >, P1Function< real_t > > >(
              storage, minLevel, maxLevel, ElementwiseMass, real_c( 2 ) );
      auto ConstantMassScaled =
          std::make_shared< ScaledOperator< P1ConstantMassOperator, P1Function< real_t >, P1Function< real_t > > >(
              storage, minLevel, maxLevel, ConstantMass, real_c( 2 ) );

      // calculate some right hand side, in this case with scaling one
      CentroidMass->apply( u, rhsElementwiseDoFValue, maxLevel, All, Replace );
      ElementwiseMass->apply( u, resElementwise, maxLevel, All, Replace );
      ConstantMass->apply( u, resConst, maxLevel, All, Replace );

      // calculate manual inverse diagonals
      CentroidMass->computeInverseDiagonalOperatorValues();
      ElementwiseMass->computeInverseDiagonalOperatorValues();
      ConstantMass->computeInverseDiagonalOperatorValues();

      // save and scale inv diags
      manualInvDiagElementwiseDoFValue.assign( { real_c( 0.5 ) }, { *CentroidMass->getInverseDiagonalValues() }, maxLevel, All );
      manualInvDiagElementwise.assign( { real_c( 0.5 ) }, { *ElementwiseMass->getInverseDiagonalValues() }, maxLevel, All );
      manualInvDiagConst.assign( { real_c( 0.5 ) }, { *ConstantMass->getInverseDiagonalValues() }, maxLevel, All );

      // calculate scaled inverse diagonals
      CentroidMassScaled->computeInverseDiagonalOperatorValues();
      ElementwiseMassScaled->computeInverseDiagonalOperatorValues();
      ConstantMassScaled->computeInverseDiagonalOperatorValues();

      // manually calculate a scaled smooth jac

      // compute the current residual
      CentroidMass->apply( u, resElementwiseDoFValue, maxLevel, All, Replace );
      resElementwiseDoFValue.assign(
          { real_c( 1 ), real_c( -2 ) }, { rhsElementwiseDoFValue, resElementwiseDoFValue }, maxLevel, All );

      ElementwiseMass->apply( u, resElementwise, maxLevel, All, Replace );
      resElementwise.assign( { real_c( 1 ), real_c( -2 ) }, { rhsElementwise, resElementwise }, maxLevel, All );

      ConstantMass->apply( u, resConst, maxLevel, All, Replace );
      resConst.assign( { real_c( 1 ), real_c( -2 ) }, { rhsConst, resConst }, maxLevel, All );

      // perform Jacobi update step
      resElementwiseDoFValue.multElementwise( { manualInvDiagElementwiseDoFValue, resElementwiseDoFValue }, maxLevel, All );
      resElementwiseDoFValue.assign( { real_c( 1 ), omega }, { u, resElementwiseDoFValue }, maxLevel, All );

      resElementwise.multElementwise( { manualInvDiagElementwise, resElementwise }, maxLevel, All );
      resElementwise.assign( { real_c( 1 ), omega }, { u, resElementwise }, maxLevel, All );

      resConst.multElementwise( { manualInvDiagConst, resConst }, maxLevel, All );
      resConst.assign( { real_c( 1 ), omega }, { u, resConst }, maxLevel, All );

      // calculate smooth_jac_scaled via scaled operator wrapper
      CentroidMassScaled->smooth_jac( resElementwiseDoFValueGEMV, rhsElementwiseDoFValue, u, omega, maxLevel, All );
      ElementwiseMassScaled->smooth_jac( resElementwiseGEMV, rhsElementwise, u, omega, maxLevel, All );
      ConstantMassScaled->smooth_jac( resConstGEMV, rhsConst, u, omega, maxLevel, All );

      // subtract results
      resElementwiseDoFValueGEMV.assign(
          { real_c( 1 ), real_c( -1 ) }, { resElementwiseDoFValueGEMV, resElementwiseDoFValue }, maxLevel, All );
      resElementwiseGEMV.assign( { real_c( 1 ), real_c( -1 ) }, { resElementwiseGEMV, resElementwise }, maxLevel, All );
      resConstGEMV.assign( { real_c( 1 ), real_c( -1 ) }, { resConstGEMV, resConst }, maxLevel, All );

      real_t ValDoFValue = resElementwiseDoFValueGEMV.dotGlobal( resElementwiseDoFValueGEMV, maxLevel, All );
      real_t ValElem     = resElementwiseGEMV.dotGlobal( resElementwiseGEMV, maxLevel, All );
      real_t ValConst    = resConstGEMV.dotGlobal( resConstGEMV, maxLevel, All );

      WALBERLA_LOG_INFO( abs( ValDoFValue ) );
      WALBERLA_CHECK_LESS( abs( ValDoFValue ), tolerance );
      WALBERLA_LOG_INFO( abs( ValElem ) );
      WALBERLA_CHECK_LESS( abs( ValElem ), tolerance );
      WALBERLA_LOG_INFO( abs( ValConst ) );
      WALBERLA_CHECK_LESS( abs( ValConst ), tolerance );

      // reset
      resElementwiseDoFValueGEMV.setToZero( maxLevel );
      resElementwiseGEMV.setToZero( maxLevel );
      resConstGEMV.setToZero( maxLevel );

      // calculate smooth_jac_scaled via scaled operator wrapper
      auto dofvalueType = hyteg::applySmoothJacScaled(
          *CentroidMass, real_c( 2 ), resElementwiseDoFValueGEMV, rhsElementwiseDoFValue, u, omega, maxLevel, All );
      auto elemType = hyteg::applySmoothJacScaled(
          *ElementwiseMass, real_c( 2 ), resElementwiseGEMV, rhsElementwise, u, omega, maxLevel, All );
      auto constType = hyteg::applySmoothJacScaled( *ConstantMass, real_c( 2 ), resConstGEMV, rhsConst, u, omega, maxLevel, All );

      // subtract results
      resElementwiseDoFValueGEMV.assign(
          { real_c( 1 ), real_c( -1 ) }, { resElementwiseDoFValueGEMV, resElementwiseDoFValue }, maxLevel, All );
      resElementwiseGEMV.assign( { real_c( 1 ), real_c( -1 ) }, { resElementwiseGEMV, resElementwise }, maxLevel, All );
      resConstGEMV.assign( { real_c( 1 ), real_c( -1 ) }, { resConstGEMV, resConst }, maxLevel, All );

      ValDoFValue = resElementwiseDoFValueGEMV.dotGlobal( resElementwiseDoFValueGEMV, maxLevel, All );
      ValElem     = resElementwiseGEMV.dotGlobal( resElementwiseGEMV, maxLevel, All );
      ValConst    = resConstGEMV.dotGlobal( resConstGEMV, maxLevel, All );

      WALBERLA_LOG_INFO( abs( ValDoFValue ) );
      WALBERLA_CHECK_LESS( abs( ValDoFValue ), tolerance );
      WALBERLA_LOG_INFO( abs( ValElem ) );
      WALBERLA_CHECK_LESS( abs( ValElem ), tolerance );
      WALBERLA_LOG_INFO( abs( ValConst ) );
      WALBERLA_CHECK_LESS( abs( ValConst ), tolerance );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled, smooth_jac_scaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( dofvalueType );
      WALBERLA_CHECK_EQUAL( dofvalueType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( elemType );
      WALBERLA_CHECK_EQUAL( elemType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( constType );
      WALBERLA_CHECK_EQUAL( constType, GEMVType::MANUAL );
   }

   {
      // P1 3D smooth_jac_scaled Test

      const uint_t minLevel  = 2;
      const uint_t maxLevel  = 2;
      const real_t tolerance = real_c( 1e-15 );
      const real_t omega     = 0.5;

      // Init setup storage
      MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid(
          Point3D( real_c( 0 ), real_c( 0 ), real_c( 0 ) ), Point3D( real_c( 1 ), real_c( 1 ), real_c( 1 ) ), 3, 3, 3 );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

      // Create storage
      std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

      // Create functions
      P2Function< real_t > T( "T", storage, minLevel, maxLevel );
      P1Function< real_t > u( "u", storage, minLevel, maxLevel );

      P1Function< real_t > rhsConst( "rhsConst", storage, minLevel, maxLevel );
      P1Function< real_t > rhsElementwise( "rhsElementwise", storage, minLevel, maxLevel );
      P1Function< real_t > rhsElementwiseDoFValue( "rhsElementwiseDoFValue", storage, minLevel, maxLevel );

      P1Function< real_t > manualInvDiagConst( "manualInvDiagConst", storage, minLevel, maxLevel );
      P1Function< real_t > manualInvDiagElementwise( "manualInvDiagElementwise", storage, minLevel, maxLevel );
      P1Function< real_t > manualInvDiagElementwiseDoFValue( "manualInvDiagElementwiseDoFValue", storage, minLevel, maxLevel );

      P1Function< real_t > resConst( "resConst", storage, minLevel, maxLevel );
      P1Function< real_t > resElementwise( "resElementwise", storage, minLevel, maxLevel );
      P1Function< real_t > resElementwiseDoFValue( "resElementwiseDoFValue", storage, minLevel, maxLevel );

      P1Function< real_t > resConstGEMV( "resConstGEMV", storage, minLevel, maxLevel );
      P1Function< real_t > resElementwiseGEMV( "resElementwiseGEMV", storage, minLevel, maxLevel );
      P1Function< real_t > resElementwiseDoFValueGEMV( "resElementwiseDoFValueGEMV", storage, minLevel, maxLevel );

      std::function< real_t( const hyteg::Point3D& ) > uFct = []( const hyteg::Point3D& x ) {
         WALBERLA_UNUSED( x );
         return std::sin( x[0] ) + x[1] * x[1] + cos( x[2] );
      };

      u.interpolate( uFct, maxLevel, hyteg::All );

      // operators
      std::function< real_t( const Point3D&, real_t, real_t, real_t, real_t, real_t ) > massScaling =
          [=]( const hyteg::Point3D& x, real_t temp, real_t tempCentroid, real_t centX, real_t centY, real_t centZ ) {
             return real_c( 1 );
          };

      auto CentroidMass = std::make_shared< p1_k_mass_visc_centroid_blending_q4_ElementwiseOperator >(
          storage, minLevel, maxLevel, T, massScaling );
      auto ElementwiseMass = std::make_shared< P1ElementwiseMassOperator >( storage, minLevel, maxLevel );
      auto ConstantMass    = std::make_shared< P1ConstantMassOperator >( storage, minLevel, maxLevel );

      // scaled operators
      auto CentroidMassScaled = std::make_shared<
          ScaledOperator< p1_k_mass_visc_centroid_blending_q4_ElementwiseOperator, P1Function< real_t >, P1Function< real_t > > >(
          storage, minLevel, maxLevel, CentroidMass, real_c( 2 ) );
      auto ElementwiseMassScaled =
          std::make_shared< ScaledOperator< P1ElementwiseMassOperator, P1Function< real_t >, P1Function< real_t > > >(
              storage, minLevel, maxLevel, ElementwiseMass, real_c( 2 ) );
      auto ConstantMassScaled =
          std::make_shared< ScaledOperator< P1ConstantMassOperator, P1Function< real_t >, P1Function< real_t > > >(
              storage, minLevel, maxLevel, ConstantMass, real_c( 2 ) );

      // calculate some right hand side, in this case with scaling one
      CentroidMass->apply( u, rhsElementwiseDoFValue, maxLevel, All, Replace );
      ElementwiseMass->apply( u, resElementwise, maxLevel, All, Replace );
      ConstantMass->apply( u, resConst, maxLevel, All, Replace );

      // calculate manual inverse diagonals
      CentroidMass->computeInverseDiagonalOperatorValues();
      ElementwiseMass->computeInverseDiagonalOperatorValues();
      ConstantMass->computeInverseDiagonalOperatorValues();

      // save and scale inv diags
      manualInvDiagElementwiseDoFValue.assign( { real_c( 0.5 ) }, { *CentroidMass->getInverseDiagonalValues() }, maxLevel, All );
      manualInvDiagElementwise.assign( { real_c( 0.5 ) }, { *ElementwiseMass->getInverseDiagonalValues() }, maxLevel, All );
      manualInvDiagConst.assign( { real_c( 0.5 ) }, { *ConstantMass->getInverseDiagonalValues() }, maxLevel, All );

      // calculate scaled inverse diagonals
      CentroidMassScaled->computeInverseDiagonalOperatorValues();
      ElementwiseMassScaled->computeInverseDiagonalOperatorValues();
      ConstantMassScaled->computeInverseDiagonalOperatorValues();

      // manually calculate a scaled smooth jac

      // compute the current residual
      CentroidMass->apply( u, resElementwiseDoFValue, maxLevel, All, Replace );
      resElementwiseDoFValue.assign(
          { real_c( 1 ), real_c( -2 ) }, { rhsElementwiseDoFValue, resElementwiseDoFValue }, maxLevel, All );

      ElementwiseMass->apply( u, resElementwise, maxLevel, All, Replace );
      resElementwise.assign( { real_c( 1 ), real_c( -2 ) }, { rhsElementwise, resElementwise }, maxLevel, All );

      ConstantMass->apply( u, resConst, maxLevel, All, Replace );
      resConst.assign( { real_c( 1 ), real_c( -2 ) }, { rhsConst, resConst }, maxLevel, All );

      // perform Jacobi update step
      resElementwiseDoFValue.multElementwise( { manualInvDiagElementwiseDoFValue, resElementwiseDoFValue }, maxLevel, All );
      resElementwiseDoFValue.assign( { real_c( 1 ), omega }, { u, resElementwiseDoFValue }, maxLevel, All );

      resElementwise.multElementwise( { manualInvDiagElementwise, resElementwise }, maxLevel, All );
      resElementwise.assign( { real_c( 1 ), omega }, { u, resElementwise }, maxLevel, All );

      resConst.multElementwise( { manualInvDiagConst, resConst }, maxLevel, All );
      resConst.assign( { real_c( 1 ), omega }, { u, resConst }, maxLevel, All );

      // calculate smooth_jac_scaled via scaled operator wrapper
      CentroidMassScaled->smooth_jac( resElementwiseDoFValueGEMV, rhsElementwiseDoFValue, u, omega, maxLevel, All );
      ElementwiseMassScaled->smooth_jac( resElementwiseGEMV, rhsElementwise, u, omega, maxLevel, All );
      ConstantMassScaled->smooth_jac( resConstGEMV, rhsConst, u, omega, maxLevel, All );

      // subtract results
      resElementwiseDoFValueGEMV.assign(
          { real_c( 1 ), real_c( -1 ) }, { resElementwiseDoFValueGEMV, resElementwiseDoFValue }, maxLevel, All );
      resElementwiseGEMV.assign( { real_c( 1 ), real_c( -1 ) }, { resElementwiseGEMV, resElementwise }, maxLevel, All );
      resConstGEMV.assign( { real_c( 1 ), real_c( -1 ) }, { resConstGEMV, resConst }, maxLevel, All );

      real_t ValDoFValue = resElementwiseDoFValueGEMV.dotGlobal( resElementwiseDoFValueGEMV, maxLevel, All );
      real_t ValElem     = resElementwiseGEMV.dotGlobal( resElementwiseGEMV, maxLevel, All );
      real_t ValConst    = resConstGEMV.dotGlobal( resConstGEMV, maxLevel, All );

      WALBERLA_LOG_INFO( abs( ValDoFValue ) );
      WALBERLA_CHECK_LESS( abs( ValDoFValue ), tolerance );
      WALBERLA_LOG_INFO( abs( ValElem ) );
      WALBERLA_CHECK_LESS( abs( ValElem ), tolerance );
      WALBERLA_LOG_INFO( abs( ValConst ) );
      WALBERLA_CHECK_LESS( abs( ValConst ), tolerance );

      // reset
      resElementwiseDoFValueGEMV.setToZero( maxLevel );
      resElementwiseGEMV.setToZero( maxLevel );
      resConstGEMV.setToZero( maxLevel );

      // calculate smooth_jac_scaled via scaled operator wrapper
      auto dofvalueType = hyteg::applySmoothJacScaled(
          *CentroidMass, real_c( 2 ), resElementwiseDoFValueGEMV, rhsElementwiseDoFValue, u, omega, maxLevel, All );
      auto elemType = hyteg::applySmoothJacScaled(
          *ElementwiseMass, real_c( 2 ), resElementwiseGEMV, rhsElementwise, u, omega, maxLevel, All );
      auto constType = hyteg::applySmoothJacScaled( *ConstantMass, real_c( 2 ), resConstGEMV, rhsConst, u, omega, maxLevel, All );

      // subtract results
      resElementwiseDoFValueGEMV.assign(
          { real_c( 1 ), real_c( -1 ) }, { resElementwiseDoFValueGEMV, resElementwiseDoFValue }, maxLevel, All );
      resElementwiseGEMV.assign( { real_c( 1 ), real_c( -1 ) }, { resElementwiseGEMV, resElementwise }, maxLevel, All );
      resConstGEMV.assign( { real_c( 1 ), real_c( -1 ) }, { resConstGEMV, resConst }, maxLevel, All );

      ValDoFValue = resElementwiseDoFValueGEMV.dotGlobal( resElementwiseDoFValueGEMV, maxLevel, All );
      ValElem     = resElementwiseGEMV.dotGlobal( resElementwiseGEMV, maxLevel, All );
      ValConst    = resConstGEMV.dotGlobal( resConstGEMV, maxLevel, All );

      WALBERLA_LOG_INFO( abs( ValDoFValue ) );
      WALBERLA_CHECK_LESS( abs( ValDoFValue ), tolerance );
      WALBERLA_LOG_INFO( abs( ValElem ) );
      WALBERLA_CHECK_LESS( abs( ValElem ), tolerance );
      WALBERLA_LOG_INFO( abs( ValConst ) );
      WALBERLA_CHECK_LESS( abs( ValConst ), tolerance );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled, smooth_jac_scaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( dofvalueType );
      WALBERLA_CHECK_EQUAL( dofvalueType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( elemType );
      WALBERLA_CHECK_EQUAL( elemType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( constType );
      WALBERLA_CHECK_EQUAL( constType, GEMVType::MANUAL );
   }

   {
      // P2 2D smooth_jac_scaled Test

      const uint_t minLevel  = 3;
      const uint_t maxLevel  = 3;
      const real_t tolerance = real_c( 1e-15 );
      const real_t omega     = 0.5;

      // Init setup storage
      MeshInfo meshInfo = MeshInfo::meshRectangle(
          Point2D( real_c( 0 ), real_c( 0 ) ), Point2D( real_c( 1 ), real_c( 1 ) ), MeshInfo::CROSS, 3, 3 );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

      // Create storage
      std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

      // Create functions
      P2Function< real_t > T( "T", storage, minLevel, maxLevel );
      P2Function< real_t > u( "u", storage, minLevel, maxLevel );

      P2Function< real_t > rhsConst( "rhsConst", storage, minLevel, maxLevel );
      P2Function< real_t > rhsElementwise( "rhsElementwise", storage, minLevel, maxLevel );
      P2Function< real_t > rhsElementwiseDoFValue( "rhsElementwiseDoFValue", storage, minLevel, maxLevel );

      P2Function< real_t > manualInvDiagConst( "manualInvDiagConst", storage, minLevel, maxLevel );
      P2Function< real_t > manualInvDiagElementwise( "manualInvDiagElementwise", storage, minLevel, maxLevel );
      P2Function< real_t > manualInvDiagElementwiseDoFValue( "manualInvDiagElementwiseDoFValue", storage, minLevel, maxLevel );

      P2Function< real_t > resConst( "resConst", storage, minLevel, maxLevel );
      P2Function< real_t > resElementwise( "resElementwise", storage, minLevel, maxLevel );
      P2Function< real_t > resElementwiseDoFValue( "resElementwiseDoFValue", storage, minLevel, maxLevel );

      P2Function< real_t > resConstGEMV( "resConstGEMV", storage, minLevel, maxLevel );
      P2Function< real_t > resElementwiseGEMV( "resElementwiseGEMV", storage, minLevel, maxLevel );
      P2Function< real_t > resElementwiseDoFValueGEMV( "resElementwiseDoFValueGEMV", storage, minLevel, maxLevel );

      std::function< real_t( const hyteg::Point3D& ) > uFct = []( const hyteg::Point3D& x ) {
         WALBERLA_UNUSED( x );
         return std::sin( x[0] ) + x[1] * x[1];
      };

      u.interpolate( uFct, maxLevel, hyteg::All );

      // operators
      std::function< real_t( const Point3D&, real_t, real_t, real_t, real_t, real_t ) > massScaling =
          [=]( const hyteg::Point3D& x, real_t temp, real_t tempCentroid, real_t centX, real_t centY, real_t centZ ) {
             return real_c( 1 );
          };

      auto CentroidMass = std::make_shared< p2_k_mass_visc_centroid_blending_q6_ElementwiseOperator >(
          storage, minLevel, maxLevel, T, massScaling );
      auto ElementwiseMass = std::make_shared< P2ElementwiseMassOperator >( storage, minLevel, maxLevel );
      auto ConstantMass    = std::make_shared< P2ConstantMassOperator >( storage, minLevel, maxLevel );

      // scaled operators
      auto CentroidMassScaled = std::make_shared<
          ScaledOperator< p2_k_mass_visc_centroid_blending_q6_ElementwiseOperator, P2Function< real_t >, P2Function< real_t > > >(
          storage, minLevel, maxLevel, CentroidMass, real_c( 2 ) );
      auto ElementwiseMassScaled =
          std::make_shared< ScaledOperator< P2ElementwiseMassOperator, P2Function< real_t >, P2Function< real_t > > >(
              storage, minLevel, maxLevel, ElementwiseMass, real_c( 2 ) );
      auto ConstantMassScaled =
          std::make_shared< ScaledOperator< P2ConstantMassOperator, P2Function< real_t >, P2Function< real_t > > >(
              storage, minLevel, maxLevel, ConstantMass, real_c( 2 ) );

      // calculate some right hand side, in this case with scaling one
      CentroidMass->apply( u, rhsElementwiseDoFValue, maxLevel, All, Replace );
      ElementwiseMass->apply( u, resElementwise, maxLevel, All, Replace );
      ConstantMass->apply( u, resConst, maxLevel, All, Replace );

      // calculate manual inverse diagonals
      CentroidMass->computeInverseDiagonalOperatorValues();
      ElementwiseMass->computeInverseDiagonalOperatorValues();
      ConstantMass->computeInverseDiagonalOperatorValues();

      // save and scale inv diags
      manualInvDiagElementwiseDoFValue.assign( { real_c( 0.5 ) }, { *CentroidMass->getInverseDiagonalValues() }, maxLevel, All );
      manualInvDiagElementwise.assign( { real_c( 0.5 ) }, { *ElementwiseMass->getInverseDiagonalValues() }, maxLevel, All );
      manualInvDiagConst.assign( { real_c( 0.5 ) }, { *ConstantMass->getInverseDiagonalValues() }, maxLevel, All );

      // calculate scaled inverse diagonals
      CentroidMassScaled->computeInverseDiagonalOperatorValues();
      ElementwiseMassScaled->computeInverseDiagonalOperatorValues();
      ConstantMassScaled->computeInverseDiagonalOperatorValues();

      // manually calculate a scaled smooth jac

      // compute the current residual
      CentroidMass->apply( u, resElementwiseDoFValue, maxLevel, All, Replace );
      resElementwiseDoFValue.assign(
          { real_c( 1 ), real_c( -2 ) }, { rhsElementwiseDoFValue, resElementwiseDoFValue }, maxLevel, All );

      ElementwiseMass->apply( u, resElementwise, maxLevel, All, Replace );
      resElementwise.assign( { real_c( 1 ), real_c( -2 ) }, { rhsElementwise, resElementwise }, maxLevel, All );

      ConstantMass->apply( u, resConst, maxLevel, All, Replace );
      resConst.assign( { real_c( 1 ), real_c( -2 ) }, { rhsConst, resConst }, maxLevel, All );

      // perform Jacobi update step
      resElementwiseDoFValue.multElementwise( { manualInvDiagElementwiseDoFValue, resElementwiseDoFValue }, maxLevel, All );
      resElementwiseDoFValue.assign( { real_c( 1 ), omega }, { u, resElementwiseDoFValue }, maxLevel, All );

      resElementwise.multElementwise( { manualInvDiagElementwise, resElementwise }, maxLevel, All );
      resElementwise.assign( { real_c( 1 ), omega }, { u, resElementwise }, maxLevel, All );

      resConst.multElementwise( { manualInvDiagConst, resConst }, maxLevel, All );
      resConst.assign( { real_c( 1 ), omega }, { u, resConst }, maxLevel, All );

      // calculate smooth_jac_scaled via scaled operator wrapper
      CentroidMassScaled->smooth_jac( resElementwiseDoFValueGEMV, rhsElementwiseDoFValue, u, omega, maxLevel, All );
      ElementwiseMassScaled->smooth_jac( resElementwiseGEMV, rhsElementwise, u, omega, maxLevel, All );
      ConstantMassScaled->smooth_jac( resConstGEMV, rhsConst, u, omega, maxLevel, All );

      // subtract results
      resElementwiseDoFValueGEMV.assign(
          { real_c( 1 ), real_c( -1 ) }, { resElementwiseDoFValueGEMV, resElementwiseDoFValue }, maxLevel, All );
      resElementwiseGEMV.assign( { real_c( 1 ), real_c( -1 ) }, { resElementwiseGEMV, resElementwise }, maxLevel, All );
      resConstGEMV.assign( { real_c( 1 ), real_c( -1 ) }, { resConstGEMV, resConst }, maxLevel, All );

      real_t ValDoFValue = resElementwiseDoFValueGEMV.dotGlobal( resElementwiseDoFValueGEMV, maxLevel, All );
      real_t ValElem     = resElementwiseGEMV.dotGlobal( resElementwiseGEMV, maxLevel, All );
      real_t ValConst    = resConstGEMV.dotGlobal( resConstGEMV, maxLevel, All );

      WALBERLA_LOG_INFO( abs( ValDoFValue ) );
      WALBERLA_CHECK_LESS( abs( ValDoFValue ), tolerance );
      WALBERLA_LOG_INFO( abs( ValElem ) );
      WALBERLA_CHECK_LESS( abs( ValElem ), tolerance );
      WALBERLA_LOG_INFO( abs( ValConst ) );
      WALBERLA_CHECK_LESS( abs( ValConst ), tolerance );

      // reset
      resElementwiseDoFValueGEMV.setToZero( maxLevel );
      resElementwiseGEMV.setToZero( maxLevel );
      resConstGEMV.setToZero( maxLevel );

      // calculate smooth_jac_scaled via scaled operator wrapper
      auto dofvalueType = hyteg::applySmoothJacScaled(
          *CentroidMass, real_c( 2 ), resElementwiseDoFValueGEMV, rhsElementwiseDoFValue, u, omega, maxLevel, All );
      auto elemType = hyteg::applySmoothJacScaled(
          *ElementwiseMass, real_c( 2 ), resElementwiseGEMV, rhsElementwise, u, omega, maxLevel, All );
      auto constType = hyteg::applySmoothJacScaled( *ConstantMass, real_c( 2 ), resConstGEMV, rhsConst, u, omega, maxLevel, All );

      // subtract results
      resElementwiseDoFValueGEMV.assign(
          { real_c( 1 ), real_c( -1 ) }, { resElementwiseDoFValueGEMV, resElementwiseDoFValue }, maxLevel, All );
      resElementwiseGEMV.assign( { real_c( 1 ), real_c( -1 ) }, { resElementwiseGEMV, resElementwise }, maxLevel, All );
      resConstGEMV.assign( { real_c( 1 ), real_c( -1 ) }, { resConstGEMV, resConst }, maxLevel, All );

      ValDoFValue = resElementwiseDoFValueGEMV.dotGlobal( resElementwiseDoFValueGEMV, maxLevel, All );
      ValElem     = resElementwiseGEMV.dotGlobal( resElementwiseGEMV, maxLevel, All );
      ValConst    = resConstGEMV.dotGlobal( resConstGEMV, maxLevel, All );

      WALBERLA_LOG_INFO( abs( ValDoFValue ) );
      WALBERLA_CHECK_LESS( abs( ValDoFValue ), tolerance );
      WALBERLA_LOG_INFO( abs( ValElem ) );
      WALBERLA_CHECK_LESS( abs( ValElem ), tolerance );
      WALBERLA_LOG_INFO( abs( ValConst ) );
      WALBERLA_CHECK_LESS( abs( ValConst ), tolerance );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled, smooth_jac_scaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( dofvalueType );
      WALBERLA_CHECK_EQUAL( dofvalueType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( elemType );
      WALBERLA_CHECK_EQUAL( elemType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( constType );
      WALBERLA_CHECK_EQUAL( constType, GEMVType::MANUAL );
   }

   {
      // P2 3D smooth_jac_scaled Test

      const uint_t minLevel  = 2;
      const uint_t maxLevel  = 2;
      const real_t tolerance = real_c( 1e-15 );
      const real_t omega     = 0.5;

      // Init setup storage
      MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid(
          Point3D( real_c( 0 ), real_c( 0 ), real_c( 0 ) ), Point3D( real_c( 1 ), real_c( 1 ), real_c( 1 ) ), 3, 3, 3 );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

      // Create storage
      std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

      // Create functions
      P2Function< real_t > T( "T", storage, minLevel, maxLevel );
      P2Function< real_t > u( "u", storage, minLevel, maxLevel );

      P2Function< real_t > rhsConst( "rhsConst", storage, minLevel, maxLevel );
      P2Function< real_t > rhsElementwise( "rhsElementwise", storage, minLevel, maxLevel );
      P2Function< real_t > rhsElementwiseDoFValue( "rhsElementwiseDoFValue", storage, minLevel, maxLevel );

      P2Function< real_t > manualInvDiagConst( "manualInvDiagConst", storage, minLevel, maxLevel );
      P2Function< real_t > manualInvDiagElementwise( "manualInvDiagElementwise", storage, minLevel, maxLevel );
      P2Function< real_t > manualInvDiagElementwiseDoFValue( "manualInvDiagElementwiseDoFValue", storage, minLevel, maxLevel );

      P2Function< real_t > resConst( "resConst", storage, minLevel, maxLevel );
      P2Function< real_t > resElementwise( "resElementwise", storage, minLevel, maxLevel );
      P2Function< real_t > resElementwiseDoFValue( "resElementwiseDoFValue", storage, minLevel, maxLevel );

      P2Function< real_t > resConstGEMV( "resConstGEMV", storage, minLevel, maxLevel );
      P2Function< real_t > resElementwiseGEMV( "resElementwiseGEMV", storage, minLevel, maxLevel );
      P2Function< real_t > resElementwiseDoFValueGEMV( "resElementwiseDoFValueGEMV", storage, minLevel, maxLevel );

      std::function< real_t( const hyteg::Point3D& ) > uFct = []( const hyteg::Point3D& x ) {
         WALBERLA_UNUSED( x );
         return std::sin( x[0] ) + x[1] * x[1] + cos( x[2] );
      };

      u.interpolate( uFct, maxLevel, hyteg::All );

      // operators
      std::function< real_t( const Point3D&, real_t, real_t, real_t, real_t, real_t ) > massScaling =
          [=]( const hyteg::Point3D& x, real_t temp, real_t tempCentroid, real_t centX, real_t centY, real_t centZ ) {
             return real_c( 1 );
          };

      auto CentroidMass = std::make_shared< p2_k_mass_visc_centroid_blending_q6_ElementwiseOperator >(
          storage, minLevel, maxLevel, T, massScaling );
      auto ElementwiseMass = std::make_shared< P2ElementwiseMassOperator >( storage, minLevel, maxLevel );
      auto ConstantMass    = std::make_shared< P2ConstantMassOperator >( storage, minLevel, maxLevel );

      // scaled operators
      auto CentroidMassScaled = std::make_shared<
          ScaledOperator< p2_k_mass_visc_centroid_blending_q6_ElementwiseOperator, P2Function< real_t >, P2Function< real_t > > >(
          storage, minLevel, maxLevel, CentroidMass, real_c( 2 ) );
      auto ElementwiseMassScaled =
          std::make_shared< ScaledOperator< P2ElementwiseMassOperator, P2Function< real_t >, P2Function< real_t > > >(
              storage, minLevel, maxLevel, ElementwiseMass, real_c( 2 ) );
      auto ConstantMassScaled =
          std::make_shared< ScaledOperator< P2ConstantMassOperator, P2Function< real_t >, P2Function< real_t > > >(
              storage, minLevel, maxLevel, ConstantMass, real_c( 2 ) );

      // calculate some right hand side, in this case with scaling one
      CentroidMass->apply( u, rhsElementwiseDoFValue, maxLevel, All, Replace );
      ElementwiseMass->apply( u, resElementwise, maxLevel, All, Replace );
      ConstantMass->apply( u, resConst, maxLevel, All, Replace );

      // calculate manual inverse diagonals, not implemented for P2ConstantOperator in 3D yet, commented out for now
      CentroidMass->computeInverseDiagonalOperatorValues();
      ElementwiseMass->computeInverseDiagonalOperatorValues();
      //   ConstantMass->computeInverseDiagonalOperatorValues();

      // save and scale inv diags
      manualInvDiagElementwiseDoFValue.assign( { real_c( 0.5 ) }, { *CentroidMass->getInverseDiagonalValues() }, maxLevel, All );
      manualInvDiagElementwise.assign( { real_c( 0.5 ) }, { *ElementwiseMass->getInverseDiagonalValues() }, maxLevel, All );
      //   manualInvDiagConst.assign( { real_c( 0.5 ) }, { *ConstantMass->getInverseDiagonalValues() }, maxLevel, All );

      // calculate scaled inverse diagonals
      CentroidMassScaled->computeInverseDiagonalOperatorValues();
      ElementwiseMassScaled->computeInverseDiagonalOperatorValues();
      //   ConstantMassScaled->computeInverseDiagonalOperatorValues();

      // manually calculate a scaled smooth jac

      // compute the current residual
      CentroidMass->apply( u, resElementwiseDoFValue, maxLevel, All, Replace );
      resElementwiseDoFValue.assign(
          { real_c( 1 ), real_c( -2 ) }, { rhsElementwiseDoFValue, resElementwiseDoFValue }, maxLevel, All );

      ElementwiseMass->apply( u, resElementwise, maxLevel, All, Replace );
      resElementwise.assign( { real_c( 1 ), real_c( -2 ) }, { rhsElementwise, resElementwise }, maxLevel, All );

      //   ConstantMass->apply( u, resConst, maxLevel, All, Replace );
      //   resConst.assign( { real_c( 1 ), real_c( -2 ) }, { rhsConst, resConst }, maxLevel, All );

      // perform Jacobi update step
      resElementwiseDoFValue.multElementwise( { manualInvDiagElementwiseDoFValue, resElementwiseDoFValue }, maxLevel, All );
      resElementwiseDoFValue.assign( { real_c( 1 ), omega }, { u, resElementwiseDoFValue }, maxLevel, All );

      resElementwise.multElementwise( { manualInvDiagElementwise, resElementwise }, maxLevel, All );
      resElementwise.assign( { real_c( 1 ), omega }, { u, resElementwise }, maxLevel, All );

      //   resConst.multElementwise( { manualInvDiagConst, resConst }, maxLevel, All );
      //   resConst.assign( { real_c( 1 ), omega }, { u, resConst }, maxLevel, All );

      // calculate smooth_jac_scaled via scaled operator wrapper
      CentroidMassScaled->smooth_jac( resElementwiseDoFValueGEMV, rhsElementwiseDoFValue, u, omega, maxLevel, All );
      ElementwiseMassScaled->smooth_jac( resElementwiseGEMV, rhsElementwise, u, omega, maxLevel, All );
      //   ConstantMassScaled->smooth_jac( resConstGEMV, rhsConst, u, omega, maxLevel, All );

      // subtract results
      resElementwiseDoFValueGEMV.assign(
          { real_c( 1 ), real_c( -1 ) }, { resElementwiseDoFValueGEMV, resElementwiseDoFValue }, maxLevel, All );
      resElementwiseGEMV.assign( { real_c( 1 ), real_c( -1 ) }, { resElementwiseGEMV, resElementwise }, maxLevel, All );
      //   resConstGEMV.assign( { real_c( 1 ), real_c( -1 ) }, { resConstGEMV, resConst }, maxLevel, All );

      real_t ValDoFValue = resElementwiseDoFValueGEMV.dotGlobal( resElementwiseDoFValueGEMV, maxLevel, All );
      real_t ValElem     = resElementwiseGEMV.dotGlobal( resElementwiseGEMV, maxLevel, All );
      //   real_t ValConst    = resConstGEMV.dotGlobal( resConstGEMV, maxLevel, All );

      WALBERLA_LOG_INFO( abs( ValDoFValue ) );
      WALBERLA_CHECK_LESS( abs( ValDoFValue ), tolerance );
      WALBERLA_LOG_INFO( abs( ValElem ) );
      WALBERLA_CHECK_LESS( abs( ValElem ), tolerance );
      //   WALBERLA_LOG_INFO( abs( ValConst ) );
      //   WALBERLA_CHECK_LESS( abs( ValConst ), tolerance );

      // reset
      resElementwiseDoFValueGEMV.setToZero( maxLevel );
      resElementwiseGEMV.setToZero( maxLevel );
      resConstGEMV.setToZero( maxLevel );

      // calculate smooth_jac_scaled via scaled operator wrapper
      auto dofvalueType = hyteg::applySmoothJacScaled(
          *CentroidMass, real_c( 2 ), resElementwiseDoFValueGEMV, rhsElementwiseDoFValue, u, omega, maxLevel, All );
      auto elemType = hyteg::applySmoothJacScaled(
          *ElementwiseMass, real_c( 2 ), resElementwiseGEMV, rhsElementwise, u, omega, maxLevel, All );
      //   auto constType = hyteg::applySmoothJacScaled( *ConstantMass, real_c( 2 ), resConstGEMV, rhsConst, u, omega, maxLevel, All );

      // subtract results
      resElementwiseDoFValueGEMV.assign(
          { real_c( 1 ), real_c( -1 ) }, { resElementwiseDoFValueGEMV, resElementwiseDoFValue }, maxLevel, All );
      resElementwiseGEMV.assign( { real_c( 1 ), real_c( -1 ) }, { resElementwiseGEMV, resElementwise }, maxLevel, All );
      //   resConstGEMV.assign( { real_c( 1 ), real_c( -1 ) }, { resConstGEMV, resConst }, maxLevel, All );

      ValDoFValue = resElementwiseDoFValueGEMV.dotGlobal( resElementwiseDoFValueGEMV, maxLevel, All );
      ValElem     = resElementwiseGEMV.dotGlobal( resElementwiseGEMV, maxLevel, All );
      //   ValConst    = resConstGEMV.dotGlobal( resConstGEMV, maxLevel, All );

      WALBERLA_LOG_INFO( abs( ValDoFValue ) );
      WALBERLA_CHECK_LESS( abs( ValDoFValue ), tolerance );
      WALBERLA_LOG_INFO( abs( ValElem ) );
      WALBERLA_CHECK_LESS( abs( ValElem ), tolerance );
      //   WALBERLA_LOG_INFO( abs( ValConst ) );
      //   WALBERLA_CHECK_LESS( abs( ValConst ), tolerance );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled, smooth_jac_scaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( dofvalueType );
      WALBERLA_CHECK_EQUAL( dofvalueType, GEMVType::SCALED );
      WALBERLA_LOG_INFO( elemType );
      WALBERLA_CHECK_EQUAL( elemType, GEMVType::SCALED );
      //   WALBERLA_LOG_INFO( constType );
      //   WALBERLA_CHECK_EQUAL( constType, GEMVType::MANUAL );
   }

   WALBERLA_LOG_INFO( "Finish" );

   return EXIT_SUCCESS;
}