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

   void addValue( uint_t row, uint_t col, real_t value ) override { val_ += abs(value); }

   void addValues( const std::vector< uint_t >& rows,
                   const std::vector< uint_t >& cols,
                   const std::vector< real_t >& values ) override
   {
      WALBERLA_ASSERT_EQUAL( values.size(), rows.size() * cols.size() );

      for ( uint_t i = 0; i < rows.size(); i++ )
      {
         for ( uint_t j = 0; j < cols.size(); j++ )
         {
            val_ += abs(values[i * cols.size() + j]);
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

      const uint_t minLevel        = 1;
      const uint_t maxLevel        = 6;
      const real_t tolerance       = real_c( 2e-5 );
      const real_t toleranceDiag   = real_c( 1e-15 );
      const real_t toleranceMatrix = real_c( 1e-15 );

      // Init setup storage
      MeshInfo meshInfo = MeshInfo::meshRectangle(
          Point2D( real_c( 0 ), real_c( 0 ) ), Point2D( real_c( 1 ), real_c( 1 ) ), MeshInfo::CROSS, 3, 3 );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

      // Create storage
      std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

      // Create functions
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

      auto P1ElementwiseMass = std::make_shared< P1ElementwiseMassOperator >( storage, minLevel, maxLevel );
      auto P1ConstantMass    = std::make_shared< P1ConstantMassOperator >( storage, minLevel, maxLevel );

      P1ElementwiseMass->apply( u, resElementwise, maxLevel, All, Replace );
      P1ConstantMass->apply( u, resConst, maxLevel, All, Replace );

      real_t integValElem  = O.dotGlobal( resElementwise, maxLevel, All );
      real_t integValConst = O.dotGlobal( resConst, maxLevel, All );
      real_t integValCtrl  = real_c( 4 ) / real_c( 3 ) - std::cos( 1 );

      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      auto P1ElementwiseMassScaled =
          std::make_shared< ScaledOperator< P1ElementwiseMassOperator, P1Function< real_t >, P1Function< real_t > > >(
              storage, minLevel, maxLevel, P1ElementwiseMass, real_c( 2 ) );
      auto P1ConstantMassScaled =
          std::make_shared< ScaledOperator< P1ConstantMassOperator, P1Function< real_t >, P1Function< real_t > > >(
              storage, minLevel, maxLevel, P1ConstantMass, real_c( 2 ) );

      P1ElementwiseMassScaled->apply( u, resElementwise, maxLevel, All, Replace );
      P1ConstantMassScaled->apply( u, resConst, maxLevel, All, Replace );

      integValElem  = O.dotGlobal( resElementwise, maxLevel, All );
      integValConst = O.dotGlobal( resConst, maxLevel, All );
      integValCtrl  = real_c( 8 ) / real_c( 3 ) - real_c( 2 ) * std::cos( 1 );

      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      auto elementwiseType = hyteg::applyGEMV< P1ElementwiseMassOperator >(
          *P1ElementwiseMass, real_c( 2 ), u, real_c( 1 ), resElementwise, maxLevel, All );
      auto constantType =
          hyteg::applyGEMV< P1ConstantMassOperator >( *P1ConstantMass, real_c( 2 ), u, real_c( 1 ), resConst, maxLevel, All );

      integValElem  = O.dotGlobal( resElementwise, maxLevel, All );
      integValConst = O.dotGlobal( resConst, maxLevel, All );
      integValCtrl  = real_c( 16 ) / real_c( 3 ) - real_c( 4 ) * std::cos( 1 );

      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseType );
      WALBERLA_CHECK_EQUAL( elementwiseType, GEMVType::GEMV );
      WALBERLA_LOG_INFO( constantType );
      WALBERLA_CHECK_EQUAL( constantType, GEMVType::MANUAL );

      // Now check the inverse diagonal
      P1Function< real_t > elementwiseDiagStore( "elementwiseDiagStore", storage, minLevel, maxLevel );
      P1Function< real_t > constDiagStore( "constDiagStore", storage, minLevel, maxLevel );

      P1ElementwiseMass->computeInverseDiagonalOperatorValues();
      P1ConstantMass->computeInverseDiagonalOperatorValues();

      elementwiseDiagStore.assign( { real_c( 1 ) }, { *P1ElementwiseMass->getInverseDiagonalValues() }, maxLevel, All );
      constDiagStore.assign( { real_c( 1 ) }, { *P1ConstantMass->getInverseDiagonalValues() }, maxLevel, All );

      P1ElementwiseMassScaled->computeInverseDiagonalOperatorValues();
      P1ConstantMassScaled->computeInverseDiagonalOperatorValues();

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

      real_t diagValElem = P1ElementwiseMass->getInverseDiagonalValues()->dotGlobal(
          *P1ElementwiseMass->getInverseDiagonalValues(), maxLevel, All );
      real_t diagValConst = P1ConstantMassScaled->getInverseDiagonalValues()->dotGlobal(
          *P1ConstantMassScaled->getInverseDiagonalValues(), maxLevel, All );

      WALBERLA_LOG_INFO( abs( diagValElem ) );
      WALBERLA_CHECK_LESS( abs( diagValElem ), toleranceDiag );
      WALBERLA_LOG_INFO( abs( diagValConst ) );
      WALBERLA_CHECK_LESS( abs( diagValConst ), toleranceDiag );

      auto elementwiseDiagType =
          hyteg::applyComputeInverseDiagonalOperatorValuesScaled< P1ElementwiseMassOperator >( *P1ElementwiseMass, real_c( 2 ) );
      auto constantDiagType =
          hyteg::applyComputeInverseDiagonalOperatorValuesScaled< P1ConstantMassOperator >( *P1ConstantMass, real_c( 2 ) );

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

      diagValElem = P1ElementwiseMass->getInverseDiagonalValues()->dotGlobal(
          *P1ElementwiseMass->getInverseDiagonalValues(), maxLevel, All );
      diagValConst = P1ConstantMassScaled->getInverseDiagonalValues()->dotGlobal(
          *P1ConstantMassScaled->getInverseDiagonalValues(), maxLevel, All );

      WALBERLA_LOG_INFO( abs( diagValElem ) );
      WALBERLA_CHECK_LESS( abs( diagValElem ), toleranceDiag );
      WALBERLA_LOG_INFO( abs( diagValConst ) );
      WALBERLA_CHECK_LESS( abs( diagValConst ), toleranceDiag );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseDiagType );
      WALBERLA_CHECK_EQUAL( elementwiseDiagType, GEMVType::MANUAL );
      WALBERLA_LOG_INFO( constantDiagType );
      WALBERLA_CHECK_EQUAL( constantDiagType, GEMVType::MANUAL );

      // create dummy matrix proxy
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyElem;
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyConst;

      dummyMatrixProxyElem  = std::make_shared< SumMatrixProxy >( real_c( 0 ) );
      dummyMatrixProxyConst = std::make_shared< SumMatrixProxy >( real_c( 0 ) );

      P1Function< idx_t > enumerator( "enumerator", storage, minLevel, maxLevel );
      enumerator.enumerate( maxLevel );

      P1ElementwiseMass->toMatrix( dummyMatrixProxyElem, enumerator, enumerator, maxLevel, All );
      P1ConstantMass->toMatrix( dummyMatrixProxyConst, enumerator, enumerator, maxLevel, All );

      real_t matrixSumValElem  = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal();
      real_t matrixSumValConst = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal();

      P1ElementwiseMassScaled->toMatrix( dummyMatrixProxyElem, enumerator, enumerator, maxLevel, All );
      P1ConstantMassScaled->toMatrix( dummyMatrixProxyConst, enumerator, enumerator, maxLevel, All );

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

      auto elementwiseMatrixType = hyteg::applyToMatrixScaled(
          *P1ElementwiseMass, real_c( 2 ), dummyMatrixProxyElem, enumerator, enumerator, maxLevel, All );
      auto constMatrixType = hyteg::applyToMatrixScaled(
          *P1ConstantMass, real_c( 2 ), dummyMatrixProxyConst, enumerator, enumerator, maxLevel, All );

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
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseMatrixType );
      WALBERLA_CHECK_EQUAL( elementwiseMatrixType, GEMVType::MANUAL );
      WALBERLA_LOG_INFO( constMatrixType );
      WALBERLA_CHECK_EQUAL( constMatrixType, GEMVType::MANUAL );
   }

   {
      // P1 3D Scaling Test

      const uint_t minLevel        = 1;
      const uint_t maxLevel        = 4;
      const real_t tolerance       = real_c( 6e-5 );
      const real_t toleranceDiag   = real_c( 1e-15 );
      const real_t toleranceMatrix = real_c( 1e-15 );

      // Init setup storage
      MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid(
          Point3D( real_c( 0 ), real_c( 0 ), real_c( 0 ) ), Point3D( real_c( 1 ), real_c( 1 ), real_c( 1 ) ), 3, 3, 3 );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

      // Create storage
      std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

      // Create functions
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

      auto P1ElementwiseMass = std::make_shared< P1ElementwiseMassOperator >( storage, minLevel, maxLevel );
      auto P1ConstantMass    = std::make_shared< P1ConstantMassOperator >( storage, minLevel, maxLevel );

      P1ElementwiseMass->apply( u, resElementwise, maxLevel, All, Replace );
      P1ConstantMass->apply( u, resConst, maxLevel, All, Replace );

      real_t integValElem  = O.dotGlobal( resElementwise, maxLevel, All );
      real_t integValConst = O.dotGlobal( resConst, maxLevel, All );
      real_t integValCtrl  = real_c( 4 ) / real_c( 3 ) + std::sin( 1 ) - std::cos( 1 );

      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      auto P1ElementwiseMassScaled =
          std::make_shared< ScaledOperator< P1ElementwiseMassOperator, P1Function< real_t >, P1Function< real_t > > >(
              storage, minLevel, maxLevel, P1ElementwiseMass, real_c( 2 ) );
      auto P1ConstantMassScaled =
          std::make_shared< ScaledOperator< P1ConstantMassOperator, P1Function< real_t >, P1Function< real_t > > >(
              storage, minLevel, maxLevel, P1ConstantMass, real_c( 2 ) );

      P1ElementwiseMassScaled->apply( u, resElementwise, maxLevel, All, Replace );
      P1ConstantMassScaled->apply( u, resConst, maxLevel, All, Replace );

      integValElem  = O.dotGlobal( resElementwise, maxLevel, All );
      integValConst = O.dotGlobal( resConst, maxLevel, All );
      integValCtrl  = real_c( 8 ) / real_c( 3 ) + real_c( 2 ) * std::sin( 1 ) - real_c( 2 ) * std::cos( 1 );

      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      auto elementwiseType = hyteg::applyGEMV< P1ElementwiseMassOperator >(
          *P1ElementwiseMass, real_c( 2 ), u, real_c( 1 ), resElementwise, maxLevel, All );
      auto constantType =
          hyteg::applyGEMV< P1ConstantMassOperator >( *P1ConstantMass, real_c( 2 ), u, real_c( 1 ), resConst, maxLevel, All );

      integValElem  = O.dotGlobal( resElementwise, maxLevel, All );
      integValConst = O.dotGlobal( resConst, maxLevel, All );
      integValCtrl  = real_c( 16 ) / real_c( 3 ) + real_c( 4 ) * std::sin( 1 ) - real_c( 4 ) * std::cos( 1 );

      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseType );
      WALBERLA_CHECK_EQUAL( elementwiseType, GEMVType::GEMV );
      WALBERLA_LOG_INFO( constantType );
      WALBERLA_CHECK_EQUAL( constantType, GEMVType::MANUAL );

      // Now check the inverse diagonal
      P1Function< real_t > elementwiseDiagStore( "elementwiseDiagStore", storage, minLevel, maxLevel );
      P1Function< real_t > constDiagStore( "constDiagStore", storage, minLevel, maxLevel );

      P1ElementwiseMass->computeInverseDiagonalOperatorValues();
      P1ConstantMass->computeInverseDiagonalOperatorValues();

      elementwiseDiagStore.assign( { real_c( 1 ) }, { *P1ElementwiseMass->getInverseDiagonalValues() }, maxLevel, All );
      constDiagStore.assign( { real_c( 1 ) }, { *P1ConstantMass->getInverseDiagonalValues() }, maxLevel, All );

      P1ElementwiseMassScaled->computeInverseDiagonalOperatorValues();
      P1ConstantMassScaled->computeInverseDiagonalOperatorValues();

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

      real_t diagValElem = P1ElementwiseMass->getInverseDiagonalValues()->dotGlobal(
          *P1ElementwiseMass->getInverseDiagonalValues(), maxLevel, All );
      real_t diagValConst = P1ConstantMassScaled->getInverseDiagonalValues()->dotGlobal(
          *P1ConstantMassScaled->getInverseDiagonalValues(), maxLevel, All );

      WALBERLA_LOG_INFO( abs( diagValElem ) );
      WALBERLA_CHECK_LESS( abs( diagValElem ), toleranceDiag );
      WALBERLA_LOG_INFO( abs( diagValConst ) );
      WALBERLA_CHECK_LESS( abs( diagValConst ), toleranceDiag );

      auto elementwiseDiagType =
          hyteg::applyComputeInverseDiagonalOperatorValuesScaled< P1ElementwiseMassOperator >( *P1ElementwiseMass, real_c( 2 ) );
      auto constantDiagType =
          hyteg::applyComputeInverseDiagonalOperatorValuesScaled< P1ConstantMassOperator >( *P1ConstantMass, real_c( 2 ) );

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

      diagValElem = P1ElementwiseMass->getInverseDiagonalValues()->dotGlobal(
          *P1ElementwiseMass->getInverseDiagonalValues(), maxLevel, All );
      diagValConst = P1ConstantMassScaled->getInverseDiagonalValues()->dotGlobal(
          *P1ConstantMassScaled->getInverseDiagonalValues(), maxLevel, All );

      WALBERLA_LOG_INFO( abs( diagValElem ) );
      WALBERLA_CHECK_LESS( abs( diagValElem ), toleranceDiag );
      WALBERLA_LOG_INFO( abs( diagValConst ) );
      WALBERLA_CHECK_LESS( abs( diagValConst ), toleranceDiag );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseDiagType );
      WALBERLA_CHECK_EQUAL( elementwiseDiagType, GEMVType::MANUAL );
      WALBERLA_LOG_INFO( constantDiagType );
      WALBERLA_CHECK_EQUAL( constantDiagType, GEMVType::MANUAL );

      // create dummy matrix proxy
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyElem;
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyConst;

      dummyMatrixProxyElem  = std::make_shared< SumMatrixProxy >( real_c( 0 ) );
      dummyMatrixProxyConst = std::make_shared< SumMatrixProxy >( real_c( 0 ) );

      P1Function< idx_t > enumerator( "enumerator", storage, minLevel, maxLevel );
      enumerator.enumerate( maxLevel );

      P1ElementwiseMass->toMatrix( dummyMatrixProxyElem, enumerator, enumerator, maxLevel, All );
      P1ConstantMass->toMatrix( dummyMatrixProxyConst, enumerator, enumerator, maxLevel, All );

      real_t matrixSumValElem  = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal();
      real_t matrixSumValConst = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal();

      P1ElementwiseMassScaled->toMatrix( dummyMatrixProxyElem, enumerator, enumerator, maxLevel, All );
      P1ConstantMassScaled->toMatrix( dummyMatrixProxyConst, enumerator, enumerator, maxLevel, All );

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

      auto elementwiseMatrixType = hyteg::applyToMatrixScaled(
          *P1ElementwiseMass, real_c( 2 ), dummyMatrixProxyElem, enumerator, enumerator, maxLevel, All );
      auto constMatrixType = hyteg::applyToMatrixScaled(
          *P1ConstantMass, real_c( 2 ), dummyMatrixProxyConst, enumerator, enumerator, maxLevel, All );

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
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseMatrixType );
      WALBERLA_CHECK_EQUAL( elementwiseMatrixType, GEMVType::MANUAL );
      WALBERLA_LOG_INFO( constMatrixType );
      WALBERLA_CHECK_EQUAL( constMatrixType, GEMVType::MANUAL );
   }

   {
      // P2 2D Scaling Test

      const uint_t minLevel        = 1;
      const uint_t maxLevel        = 6;
      const real_t tolerance       = real_c( 1e-12 );
      const real_t toleranceDiag   = real_c( 1e-15 );
      const real_t toleranceMatrix = real_c( 1e-15 );

      // Init setup storage
      MeshInfo meshInfo = MeshInfo::meshRectangle(
          Point2D( real_c( 0 ), real_c( 0 ) ), Point2D( real_c( 1 ), real_c( 1 ) ), MeshInfo::CROSS, 3, 3 );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

      // Create storage
      std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

      // Create functions
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

      auto P2ElementwiseMass = std::make_shared< P2ElementwiseMassOperator >( storage, minLevel, maxLevel );
      auto P2ConstantMass    = std::make_shared< P2ConstantMassOperator >( storage, minLevel, maxLevel );

      P2ElementwiseMass->apply( u, resElementwise, maxLevel, All, Replace );
      P2ConstantMass->apply( u, resConst, maxLevel, All, Replace );

      real_t integValElem  = O.dotGlobal( resElementwise, maxLevel, All );
      real_t integValConst = O.dotGlobal( resConst, maxLevel, All );
      real_t integValCtrl  = real_c( 4 ) / real_c( 3 ) - std::cos( 1 );

      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      auto P2ElementwiseMassScaled =
          std::make_shared< ScaledOperator< P2ElementwiseMassOperator, P2Function< real_t >, P2Function< real_t > > >(
              storage, minLevel, maxLevel, P2ElementwiseMass, real_c( 2 ) );
      auto P2ConstantMassScaled =
          std::make_shared< ScaledOperator< P2ConstantMassOperator, P2Function< real_t >, P2Function< real_t > > >(
              storage, minLevel, maxLevel, P2ConstantMass, real_c( 2 ) );

      P2ElementwiseMassScaled->apply( u, resElementwise, maxLevel, All, Replace );
      P2ConstantMassScaled->apply( u, resConst, maxLevel, All, Replace );

      integValElem  = O.dotGlobal( resElementwise, maxLevel, All );
      integValConst = O.dotGlobal( resConst, maxLevel, All );
      integValCtrl  = real_c( 8 ) / real_c( 3 ) - real_c( 2 ) * std::cos( 1 );

      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      auto elementwiseType = hyteg::applyGEMV< P2ElementwiseMassOperator >(
          *P2ElementwiseMass, real_c( 2 ), u, real_c( 1 ), resElementwise, maxLevel, All );
      auto constantType =
          hyteg::applyGEMV< P2ConstantMassOperator >( *P2ConstantMass, real_c( 2 ), u, real_c( 1 ), resConst, maxLevel, All );

      integValElem  = O.dotGlobal( resElementwise, maxLevel, All );
      integValConst = O.dotGlobal( resConst, maxLevel, All );
      integValCtrl  = real_c( 16 ) / real_c( 3 ) - real_c( 4 ) * std::cos( 1 );

      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseType );
      WALBERLA_CHECK_EQUAL( elementwiseType, GEMVType::MANUAL );
      WALBERLA_LOG_INFO( constantType );
      WALBERLA_CHECK_EQUAL( constantType, GEMVType::MANUAL );

      // Now check the inverse diagonal, P2 constant operators for some reason do not have a computeInverseDiagonalOperatorValues method.
      // Leaving the respective calls as a comment in case they can be reactivated later on.
      P2Function< real_t > elementwiseDiagStore( "elementwiseDiagStore", storage, minLevel, maxLevel );
      //   P2Function< real_t > constDiagStore( "constDiagStore", storage, minLevel, maxLevel );

      P2ElementwiseMass->computeInverseDiagonalOperatorValues();
      //   P2ConstantMass->computeInverseDiagonalOperatorValues();

      elementwiseDiagStore.assign( { real_c( 1 ) }, { *P2ElementwiseMass->getInverseDiagonalValues() }, maxLevel, All );
      //   constDiagStore.assign( { real_c( 1 ) }, { *P2ConstantMass->getInverseDiagonalValues() }, maxLevel, All );

      P2ElementwiseMassScaled->computeInverseDiagonalOperatorValues();
      //   P2ConstantMassScaled->computeInverseDiagonalOperatorValues();

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

      real_t diagValElem = P2ElementwiseMass->getInverseDiagonalValues()->dotGlobal(
          *P2ElementwiseMass->getInverseDiagonalValues(), maxLevel, All );
      //   real_t diagValConst = P2ConstantMassScaled->getInverseDiagonalValues()->dotGlobal(
      //       *P2ConstantMassScaled->getInverseDiagonalValues(), maxLevel, All );

      WALBERLA_LOG_INFO( abs( diagValElem ) );
      WALBERLA_CHECK_LESS( abs( diagValElem ), toleranceDiag );
      //   WALBERLA_LOG_INFO( abs( diagValConst ) );
      //   WALBERLA_CHECK_LESS( abs( diagValConst ), toleranceDiag );

      auto elementwiseDiagType =
          hyteg::applyComputeInverseDiagonalOperatorValuesScaled< P2ElementwiseMassOperator >( *P2ElementwiseMass, real_c( 2 ) );
      //   auto constantDiagType =
      //       hyteg::applyComputeInverseDiagonalOperatorValuesScaled< P2ConstantMassOperator >( *P2ConstantMass, real_c( 2 ) );

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

      diagValElem = P2ElementwiseMass->getInverseDiagonalValues()->dotGlobal(
          *P2ElementwiseMass->getInverseDiagonalValues(), maxLevel, All );
      //   diagValConst = P2ConstantMassScaled->getInverseDiagonalValues()->dotGlobal(
      //       *P2ConstantMassScaled->getInverseDiagonalValues(), maxLevel, All );

      WALBERLA_LOG_INFO( abs( diagValElem ) );
      WALBERLA_CHECK_LESS( abs( diagValElem ), toleranceDiag );
      //   WALBERLA_LOG_INFO( abs( diagValConst ) );
      //   WALBERLA_CHECK_LESS( abs( diagValConst ), toleranceDiag );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseDiagType );
      WALBERLA_CHECK_EQUAL( elementwiseDiagType, GEMVType::MANUAL );
      //   WALBERLA_LOG_INFO( constantDiagType );
      //   WALBERLA_CHECK_EQUAL( constantDiagType, GEMVType::MANUAL );

      // create dummy matrix proxy
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyElem;
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyConst;

      dummyMatrixProxyElem  = std::make_shared< SumMatrixProxy >( real_c( 0 ) );
      dummyMatrixProxyConst = std::make_shared< SumMatrixProxy >( real_c( 0 ) );

      P2Function< idx_t > enumerator( "enumerator", storage, minLevel, maxLevel );
      enumerator.enumerate( maxLevel );

      P2ElementwiseMass->toMatrix( dummyMatrixProxyElem, enumerator, enumerator, maxLevel, All );
      P2ConstantMass->toMatrix( dummyMatrixProxyConst, enumerator, enumerator, maxLevel, All );

      real_t matrixSumValElem  = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal();
      real_t matrixSumValConst = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal();

      P2ElementwiseMassScaled->toMatrix( dummyMatrixProxyElem, enumerator, enumerator, maxLevel, All );
      P2ConstantMassScaled->toMatrix( dummyMatrixProxyConst, enumerator, enumerator, maxLevel, All );

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

      auto elementwiseMatrixType = hyteg::applyToMatrixScaled(
          *P2ElementwiseMass, real_c( 2 ), dummyMatrixProxyElem, enumerator, enumerator, maxLevel, All );
      auto constMatrixType = hyteg::applyToMatrixScaled(
          *P2ConstantMass, real_c( 2 ), dummyMatrixProxyConst, enumerator, enumerator, maxLevel, All );

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
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseMatrixType );
      WALBERLA_CHECK_EQUAL( elementwiseMatrixType, GEMVType::MANUAL );
      WALBERLA_LOG_INFO( constMatrixType );
      WALBERLA_CHECK_EQUAL( constMatrixType, GEMVType::MANUAL );
   }

   {
      // P2 3D Scaling Test

      const uint_t minLevel        = 1;
      const uint_t maxLevel        = 3;
      const real_t tolerance       = real_c( 3e-10 );
      const real_t toleranceDiag   = real_c( 1e-15 );
      const real_t toleranceMatrix = real_c( 1e-15 );

      // Init setup storage
      MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid(
          Point3D( real_c( 0 ), real_c( 0 ), real_c( 0 ) ), Point3D( real_c( 1 ), real_c( 1 ), real_c( 1 ) ), 3, 3, 3 );
      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

      // Create storage
      std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

      // Create functions
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

      auto P2ElementwiseMass = std::make_shared< P2ElementwiseMassOperator >( storage, minLevel, maxLevel );
      auto P2ConstantMass    = std::make_shared< P2ConstantMassOperator >( storage, minLevel, maxLevel );

      P2ElementwiseMass->apply( u, resElementwise, maxLevel, All, Replace );
      P2ConstantMass->apply( u, resConst, maxLevel, All, Replace );

      real_t integValElem  = O.dotGlobal( resElementwise, maxLevel, All );
      real_t integValConst = O.dotGlobal( resConst, maxLevel, All );
      real_t integValCtrl  = real_c( 4 ) / real_c( 3 ) + std::sin( 1 ) - std::cos( 1 );

      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      auto P2ElementwiseMassScaled =
          std::make_shared< ScaledOperator< P2ElementwiseMassOperator, P2Function< real_t >, P2Function< real_t > > >(
              storage, minLevel, maxLevel, P2ElementwiseMass, real_c( 2 ) );
      auto P2ConstantMassScaled =
          std::make_shared< ScaledOperator< P2ConstantMassOperator, P2Function< real_t >, P2Function< real_t > > >(
              storage, minLevel, maxLevel, P2ConstantMass, real_c( 2 ) );

      P2ElementwiseMassScaled->apply( u, resElementwise, maxLevel, All, Replace );
      P2ConstantMassScaled->apply( u, resConst, maxLevel, All, Replace );

      integValElem  = O.dotGlobal( resElementwise, maxLevel, All );
      integValConst = O.dotGlobal( resConst, maxLevel, All );
      integValCtrl  = real_c( 8 ) / real_c( 3 ) + real_c( 2 ) * std::sin( 1 ) - real_c( 2 ) * std::cos( 1 );

      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      auto elementwiseType = hyteg::applyGEMV< P2ElementwiseMassOperator >(
          *P2ElementwiseMass, real_c( 2 ), u, real_c( 1 ), resElementwise, maxLevel, All );
      auto constantType =
          hyteg::applyGEMV< P2ConstantMassOperator >( *P2ConstantMass, real_c( 2 ), u, real_c( 1 ), resConst, maxLevel, All );

      integValElem  = O.dotGlobal( resElementwise, maxLevel, All );
      integValConst = O.dotGlobal( resConst, maxLevel, All );
      integValCtrl  = real_c( 16 ) / real_c( 3 ) + real_c( 4 ) * std::sin( 1 ) - real_c( 4 ) * std::cos( 1 );

      WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
      WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
      WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseType );
      WALBERLA_CHECK_EQUAL( elementwiseType, GEMVType::MANUAL );
      WALBERLA_LOG_INFO( constantType );
      WALBERLA_CHECK_EQUAL( constantType, GEMVType::MANUAL );

      // Now check the inverse diagonal, P2 constant operators for some reason do not have a computeInverseDiagonalOperatorValues method.
      // Leaving the respective calls as a comment in case they can be reactivated later on.
      P2Function< real_t > elementwiseDiagStore( "elementwiseDiagStore", storage, minLevel, maxLevel );
      //   P2Function< real_t > constDiagStore( "constDiagStore", storage, minLevel, maxLevel );

      P2ElementwiseMass->computeInverseDiagonalOperatorValues();
      //   P2ConstantMass->computeInverseDiagonalOperatorValues();

      elementwiseDiagStore.assign( { real_c( 1 ) }, { *P2ElementwiseMass->getInverseDiagonalValues() }, maxLevel, All );
      //   constDiagStore.assign( { real_c( 1 ) }, { *P2ConstantMass->getInverseDiagonalValues() }, maxLevel, All );

      P2ElementwiseMassScaled->computeInverseDiagonalOperatorValues();
      //   P2ConstantMassScaled->computeInverseDiagonalOperatorValues();

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

      real_t diagValElem = P2ElementwiseMass->getInverseDiagonalValues()->dotGlobal(
          *P2ElementwiseMass->getInverseDiagonalValues(), maxLevel, All );
      //   real_t diagValConst = P2ConstantMassScaled->getInverseDiagonalValues()->dotGlobal(
      //       *P2ConstantMassScaled->getInverseDiagonalValues(), maxLevel, All );

      WALBERLA_LOG_INFO( abs( diagValElem ) );
      WALBERLA_CHECK_LESS( abs( diagValElem ), toleranceDiag );
      //   WALBERLA_LOG_INFO( abs( diagValConst ) );
      //   WALBERLA_CHECK_LESS( abs( diagValConst ), toleranceDiag );

      auto elementwiseDiagType =
          hyteg::applyComputeInverseDiagonalOperatorValuesScaled< P2ElementwiseMassOperator >( *P2ElementwiseMass, real_c( 2 ) );
      //   auto constantDiagType =
      //       hyteg::applyComputeInverseDiagonalOperatorValuesScaled< P2ConstantMassOperator >( *P2ConstantMass, real_c( 2 ) );

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

      diagValElem = P2ElementwiseMass->getInverseDiagonalValues()->dotGlobal(
          *P2ElementwiseMass->getInverseDiagonalValues(), maxLevel, All );
      //   diagValConst = P2ConstantMassScaled->getInverseDiagonalValues()->dotGlobal(
      //       *P2ConstantMassScaled->getInverseDiagonalValues(), maxLevel, All );

      WALBERLA_LOG_INFO( abs( diagValElem ) );
      WALBERLA_CHECK_LESS( abs( diagValElem ), toleranceDiag );
      //   WALBERLA_LOG_INFO( abs( diagValConst ) );
      //   WALBERLA_CHECK_LESS( abs( diagValConst ), toleranceDiag );

      // If this throws an error and you have just added the functionality that one of the operator types can use applyScaled
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseDiagType );
      WALBERLA_CHECK_EQUAL( elementwiseDiagType, GEMVType::MANUAL );
      //   WALBERLA_LOG_INFO( constantDiagType );
      //   WALBERLA_CHECK_EQUAL( constantDiagType, GEMVType::MANUAL );

      // create dummy matrix proxy
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyElem;
      std::shared_ptr< SparseMatrixProxy > dummyMatrixProxyConst;

      dummyMatrixProxyElem  = std::make_shared< SumMatrixProxy >( real_c( 0 ) );
      dummyMatrixProxyConst = std::make_shared< SumMatrixProxy >( real_c( 0 ) );

      P2Function< idx_t > enumerator( "enumerator", storage, minLevel, maxLevel );
      enumerator.enumerate( maxLevel );

      P2ElementwiseMass->toMatrix( dummyMatrixProxyElem, enumerator, enumerator, maxLevel, All );
      P2ConstantMass->toMatrix( dummyMatrixProxyConst, enumerator, enumerator, maxLevel, All );

      real_t matrixSumValElem  = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyElem )->getVal();
      real_t matrixSumValConst = std::dynamic_pointer_cast< SumMatrixProxy >( dummyMatrixProxyConst )->getVal();

      P2ElementwiseMassScaled->toMatrix( dummyMatrixProxyElem, enumerator, enumerator, maxLevel, All );
      P2ConstantMassScaled->toMatrix( dummyMatrixProxyConst, enumerator, enumerator, maxLevel, All );

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

      auto elementwiseMatrixType = hyteg::applyToMatrixScaled(
          *P2ElementwiseMass, real_c( 2 ), dummyMatrixProxyElem, enumerator, enumerator, maxLevel, All );
      auto constMatrixType = hyteg::applyToMatrixScaled(
          *P2ConstantMass, real_c( 2 ), dummyMatrixProxyConst, enumerator, enumerator, maxLevel, All );

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
      // computeInverseDiagonalOperatorValuesScaled, toMatrixScaled or GEMV, please update the check here.
      // It should check that the operator uses the best possible method. Thanks!
      WALBERLA_LOG_INFO( elementwiseMatrixType );
      WALBERLA_CHECK_EQUAL( elementwiseMatrixType, GEMVType::MANUAL );
      WALBERLA_LOG_INFO( constMatrixType );
      WALBERLA_CHECK_EQUAL( constMatrixType, GEMVType::MANUAL );
   }

   //    {
   //       // P1 To P2 2D Scaling Test, No Blending

   //       // Extract the required parameters
   //       const uint_t minLevel  = 1;
   //       const uint_t maxLevel  = 6;
   //       const real_t tolerance = real_c( 1e-12 );

   //       // Init setup storage
   //       MeshInfo meshInfo = MeshInfo::meshRectangle(
   //           Point2D( real_c( 0 ), real_c( 0 ) ), Point2D( real_c( 1 ), real_c( 1 ) ), MeshInfo::CROSS, 3, 3 );
   //       SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   //       // Create storage
   //       std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   //       // Create functions
   //       P1Function< real_t > u( "u", storage, minLevel, maxLevel );
   //       P2Function< real_t > resConst( "resConst", storage, minLevel, maxLevel );
   //       P2Function< real_t > resElementwise( "resElementwise", storage, minLevel, maxLevel );
   //       P2Function< real_t > O( "O", storage, minLevel, maxLevel );

   //       std::function< real_t( const hyteg::Point3D& ) > uFct = []( const hyteg::Point3D& x ) {
   //          WALBERLA_UNUSED( x );
   //          return std::sin( x[0] ) + x[1] * x[1];
   //       };

   //       u.interpolate( uFct, maxLevel, hyteg::All );
   //       O.interpolate( real_c( 1 ), maxLevel, hyteg::All );

   //       auto divXConst = std::make_shared< P1ToP2ConstantDivxOperator >( storage, minLevel, maxLevel );
   //       auto divXElem  = std::make_shared< P1ToP2ElementwiseDivxOperator >( storage, minLevel, maxLevel );
   //       auto divYConst = std::make_shared< P1ToP2ConstantDivyOperator >( storage, minLevel, maxLevel );
   //       auto divYElem  = std::make_shared< P1ToP2ElementwiseDivyOperator >( storage, minLevel, maxLevel );

   //       divXConst->apply( u, resConst, maxLevel, All, Replace );
   //       divXElem->apply( u, resElementwise, maxLevel, All, Replace );
   //       divYConst->apply( u, resConst, maxLevel, All, Add );
   //       divYElem->apply( u, resElementwise, maxLevel, All, Add );

   //       real_t integValElem  = O.dotGlobal( resElementwise, maxLevel, All );
   //       real_t integValConst = O.dotGlobal( resConst, maxLevel, All );
   //       real_t integValCtrl  = real_c( -1.841470984808 );

   //       WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
   //       WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
   //       WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
   //       WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

   //       divXConst->swapConstantScalingSign();
   //       divXElem->swapConstantScalingSign();
   //       divYConst->swapConstantScalingSign();
   //       divYElem->swapConstantScalingSign();
   //       divXConst->apply( u, resConst, maxLevel, All, Add );
   //       divXElem->apply( u, resElementwise, maxLevel, All, Add );
   //       divYConst->apply( u, resConst, maxLevel, All, Add );
   //       divYElem->apply( u, resElementwise, maxLevel, All, Add );

   //       real_t normElem  = resElementwise.dotGlobal( resElementwise, maxLevel, All );
   //       real_t normConst = resConst.dotGlobal( resConst, maxLevel, All );
   //       real_t normCTRL  = real_c( 0 );

   //       WALBERLA_LOG_INFO( abs( normElem - normCTRL ) );
   //       WALBERLA_CHECK_LESS( abs( normElem - normCTRL ), real_c( 1e-16 ) );
   //       WALBERLA_LOG_INFO( abs( normConst - normCTRL ) );
   //       WALBERLA_CHECK_LESS( abs( normConst - normCTRL ), real_c( 1e-16 ) );

   //       divXConst->setConstantScaling( real_c( 2 ) );
   //       divXElem->setConstantScaling( real_c( 2 ) );
   //       divYConst->setConstantScaling( real_c( 2 ) );
   //       divYElem->setConstantScaling( real_c( 2 ) );
   //       divXConst->apply( u, resConst, maxLevel, All, Replace );
   //       divXElem->apply( u, resElementwise, maxLevel, All, Replace );
   //       divYConst->apply( u, resConst, maxLevel, All, Add );
   //       divYElem->apply( u, resElementwise, maxLevel, All, Add );

   //       integValElem  = O.dotGlobal( resElementwise, maxLevel, All );
   //       integValConst = O.dotGlobal( resConst, maxLevel, All );
   //       integValCtrl  = real_c( -3.682941969616 );

   //       WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
   //       WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
   //       WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
   //       WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );
   //    }

   //    {
   //       // P1 To P2 3D Scaling Test, No Blending

   //       // Extract the required parameters
   //       const uint_t minLevel  = 1;
   //       const uint_t maxLevel  = 4;
   //       const real_t tolerance = real_c( 1e-11 );

   //       // Init setup storage
   //       MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid(
   //           Point3D( real_c( 0 ), real_c( 0 ), real_c( 0 ) ), Point3D( real_c( 1 ), real_c( 1 ), real_c( 1 ) ), 3, 3, 3 );
   //       SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   //       // Create storage
   //       std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   //       // Create functions
   //       P1Function< real_t > u( "u", storage, minLevel, maxLevel );
   //       P2Function< real_t > resConst( "resConst", storage, minLevel, maxLevel );
   //       P2Function< real_t > resElementwise( "resElementwise", storage, minLevel, maxLevel );
   //       P2Function< real_t > O( "O", storage, minLevel, maxLevel );

   //       std::function< real_t( const hyteg::Point3D& ) > uFct = []( const hyteg::Point3D& x ) {
   //          WALBERLA_UNUSED( x );
   //          return std::sin( x[0] ) + x[1] * x[1] + cos( x[2] );
   //       };

   //       u.interpolate( uFct, maxLevel, hyteg::All );
   //       O.interpolate( real_c( 1 ), maxLevel, hyteg::All );

   //       auto divXConst = std::make_shared< P1ToP2ConstantDivxOperator >( storage, minLevel, maxLevel );
   //       auto divXElem  = std::make_shared< P1ToP2ElementwiseDivxOperator >( storage, minLevel, maxLevel );
   //       auto divYConst = std::make_shared< P1ToP2ConstantDivyOperator >( storage, minLevel, maxLevel );
   //       auto divYElem  = std::make_shared< P1ToP2ElementwiseDivyOperator >( storage, minLevel, maxLevel );
   //       auto divZConst = std::make_shared< P1ToP2ConstantDivzOperator >( storage, minLevel, maxLevel );
   //       auto divZElem  = std::make_shared< P1ToP2ElementwiseDivzOperator >( storage, minLevel, maxLevel );

   //       divXConst->apply( u, resConst, maxLevel, All, Replace );
   //       divXElem->apply( u, resElementwise, maxLevel, All, Replace );
   //       divYConst->apply( u, resConst, maxLevel, All, Add );
   //       divYElem->apply( u, resElementwise, maxLevel, All, Add );
   //       divZConst->apply( u, resConst, maxLevel, All, Add );
   //       divZElem->apply( u, resElementwise, maxLevel, All, Add );

   //       real_t integValElem  = O.dotGlobal( resElementwise, maxLevel, All );
   //       real_t integValConst = O.dotGlobal( resConst, maxLevel, All );
   //       real_t integValCtrl  = real_c( -1.381773290675 );

   //       WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
   //       WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
   //       WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
   //       WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

   //       divXConst->swapConstantScalingSign();
   //       divXElem->swapConstantScalingSign();
   //       divYConst->swapConstantScalingSign();
   //       divYElem->swapConstantScalingSign();
   //       divZConst->swapConstantScalingSign();
   //       divZElem->swapConstantScalingSign();
   //       divXConst->apply( u, resConst, maxLevel, All, Add );
   //       divXElem->apply( u, resElementwise, maxLevel, All, Add );
   //       divYConst->apply( u, resConst, maxLevel, All, Add );
   //       divYElem->apply( u, resElementwise, maxLevel, All, Add );
   //       divZConst->apply( u, resConst, maxLevel, All, Add );
   //       divZElem->apply( u, resElementwise, maxLevel, All, Add );

   //       real_t normElem  = resElementwise.dotGlobal( resElementwise, maxLevel, All );
   //       real_t normConst = resConst.dotGlobal( resConst, maxLevel, All );
   //       real_t normCTRL  = real_c( 0 );

   //       WALBERLA_LOG_INFO( abs( normElem - normCTRL ) );
   //       WALBERLA_CHECK_LESS( abs( normElem - normCTRL ), real_c( 1e-16 ) );
   //       WALBERLA_LOG_INFO( abs( normConst - normCTRL ) );
   //       WALBERLA_CHECK_LESS( abs( normConst - normCTRL ), real_c( 1e-16 ) );

   //       divXConst->setConstantScaling( real_c( 2 ) );
   //       divXElem->setConstantScaling( real_c( 2 ) );
   //       divYConst->setConstantScaling( real_c( 2 ) );
   //       divYElem->setConstantScaling( real_c( 2 ) );
   //       divZConst->setConstantScaling( real_c( 2 ) );
   //       divZElem->setConstantScaling( real_c( 2 ) );
   //       divXConst->apply( u, resConst, maxLevel, All, Replace );
   //       divXElem->apply( u, resElementwise, maxLevel, All, Replace );
   //       divYConst->apply( u, resConst, maxLevel, All, Add );
   //       divYElem->apply( u, resElementwise, maxLevel, All, Add );
   //       divZConst->apply( u, resConst, maxLevel, All, Add );
   //       divZElem->apply( u, resElementwise, maxLevel, All, Add );

   //       integValElem  = O.dotGlobal( resElementwise, maxLevel, All );
   //       integValConst = O.dotGlobal( resConst, maxLevel, All );
   //       integValCtrl  = real_c( -2.763546581350 );

   //       WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
   //       WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
   //       WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
   //       WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );
   //    }

   //    {
   //       // P2 To P1 2D Scaling Test, No Blending

   //       // Extract the required parameters
   //       const uint_t minLevel  = 1;
   //       const uint_t maxLevel  = 6;
   //       const real_t tolerance = real_c( 1e-10 );

   //       // Init setup storage
   //       MeshInfo meshInfo = MeshInfo::meshRectangle(
   //           Point2D( real_c( 0 ), real_c( 0 ) ), Point2D( real_c( 1 ), real_c( 1 ) ), MeshInfo::CROSS, 3, 3 );
   //       SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   //       // Create storage
   //       std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   //       // Create functions
   //       P2Function< real_t > u( "u", storage, minLevel, maxLevel );
   //       P1Function< real_t > resConst( "resConst", storage, minLevel, maxLevel );
   //       P1Function< real_t > resElementwise( "resElementwise", storage, minLevel, maxLevel );
   //       P1Function< real_t > O( "O", storage, minLevel, maxLevel );

   //       std::function< real_t( const hyteg::Point3D& ) > uFct = []( const hyteg::Point3D& x ) {
   //          WALBERLA_UNUSED( x );
   //          return std::sin( x[0] ) + x[1] * x[1];
   //       };

   //       u.interpolate( uFct, maxLevel, hyteg::All );
   //       O.interpolate( real_c( 1 ), maxLevel, hyteg::All );

   //       auto divXConst = std::make_shared< P2ToP1ConstantDivxOperator >( storage, minLevel, maxLevel );
   //       auto divXElem  = std::make_shared< P2ToP1ElementwiseDivxOperator >( storage, minLevel, maxLevel );
   //       auto divYConst = std::make_shared< P2ToP1ConstantDivyOperator >( storage, minLevel, maxLevel );
   //       auto divYElem  = std::make_shared< P2ToP1ElementwiseDivyOperator >( storage, minLevel, maxLevel );

   //       divXConst->apply( u, resConst, maxLevel, All, Replace );
   //       divXElem->apply( u, resElementwise, maxLevel, All, Replace );
   //       divYConst->apply( u, resConst, maxLevel, All, Add );
   //       divYElem->apply( u, resElementwise, maxLevel, All, Add );

   //       real_t integValElem  = O.dotGlobal( resElementwise, maxLevel, All );
   //       real_t integValConst = O.dotGlobal( resConst, maxLevel, All );
   //       real_t integValCtrl  = real_c( -1.841470984785 );

   //       WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
   //       WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
   //       WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
   //       WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

   //       divXConst->swapConstantScalingSign();
   //       divXElem->swapConstantScalingSign();
   //       divYConst->swapConstantScalingSign();
   //       divYElem->swapConstantScalingSign();
   //       divXConst->apply( u, resConst, maxLevel, All, Add );
   //       divXElem->apply( u, resElementwise, maxLevel, All, Add );
   //       divYConst->apply( u, resConst, maxLevel, All, Add );
   //       divYElem->apply( u, resElementwise, maxLevel, All, Add );

   //       real_t normElem  = resElementwise.dotGlobal( resElementwise, maxLevel, All );
   //       real_t normConst = resConst.dotGlobal( resConst, maxLevel, All );
   //       real_t normCTRL  = real_c( 0 );

   //       WALBERLA_LOG_INFO( abs( normElem - normCTRL ) );
   //       WALBERLA_CHECK_LESS( abs( normElem - normCTRL ), real_c( 1e-16 ) );
   //       WALBERLA_LOG_INFO( abs( normConst - normCTRL ) );
   //       WALBERLA_CHECK_LESS( abs( normConst - normCTRL ), real_c( 1e-16 ) );

   //       divXConst->setConstantScaling( real_c( 2 ) );
   //       divXElem->setConstantScaling( real_c( 2 ) );
   //       divYConst->setConstantScaling( real_c( 2 ) );
   //       divYElem->setConstantScaling( real_c( 2 ) );
   //       divXConst->apply( u, resConst, maxLevel, All, Replace );
   //       divXElem->apply( u, resElementwise, maxLevel, All, Replace );
   //       divYConst->apply( u, resConst, maxLevel, All, Add );
   //       divYElem->apply( u, resElementwise, maxLevel, All, Add );

   //       integValElem  = O.dotGlobal( resElementwise, maxLevel, All );
   //       integValConst = O.dotGlobal( resConst, maxLevel, All );
   //       integValCtrl  = real_c( -3.682941969570 );

   //       WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
   //       WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
   //       WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
   //       WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );
   //    }

   //    {
   //       // P2 To P1 3D Scaling Test, No Blending

   //       // Extract the required parameters
   //       const uint_t minLevel  = 1;
   //       const uint_t maxLevel  = 4;
   //       const real_t tolerance = real_c( 1e-11 );

   //       // Init setup storage
   //       MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid(
   //           Point3D( real_c( 0 ), real_c( 0 ), real_c( 0 ) ), Point3D( real_c( 1 ), real_c( 1 ), real_c( 1 ) ), 3, 3, 3 );
   //       SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   //       // Create storage
   //       std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   //       // Create functions
   //       P2Function< real_t > u( "u", storage, minLevel, maxLevel );
   //       P1Function< real_t > resConst( "resConst", storage, minLevel, maxLevel );
   //       P1Function< real_t > resElementwise( "resElementwise", storage, minLevel, maxLevel );
   //       P1Function< real_t > O( "O", storage, minLevel, maxLevel );

   //       std::function< real_t( const hyteg::Point3D& ) > uFct = []( const hyteg::Point3D& x ) {
   //          WALBERLA_UNUSED( x );
   //          return std::sin( x[0] ) + x[1] * x[1] + cos( x[2] );
   //       };

   //       u.interpolate( uFct, maxLevel, hyteg::All );
   //       O.interpolate( real_c( 1 ), maxLevel, hyteg::All );

   //       auto divXConst = std::make_shared< P2ToP1ConstantDivxOperator >( storage, minLevel, maxLevel );
   //       auto divXElem  = std::make_shared< P2ToP1ElementwiseDivxOperator >( storage, minLevel, maxLevel );
   //       auto divYConst = std::make_shared< P2ToP1ConstantDivyOperator >( storage, minLevel, maxLevel );
   //       auto divYElem  = std::make_shared< P2ToP1ElementwiseDivyOperator >( storage, minLevel, maxLevel );
   //       auto divZConst = std::make_shared< P2ToP1ConstantDivzOperator >( storage, minLevel, maxLevel );
   //       auto divZElem  = std::make_shared< P2ToP1ElementwiseDivzOperator >( storage, minLevel, maxLevel );

   //       divXConst->apply( u, resConst, maxLevel, All, Replace );
   //       divXElem->apply( u, resElementwise, maxLevel, All, Replace );
   //       divYConst->apply( u, resConst, maxLevel, All, Add );
   //       divYElem->apply( u, resElementwise, maxLevel, All, Add );
   //       divZConst->apply( u, resConst, maxLevel, All, Add );
   //       divZElem->apply( u, resElementwise, maxLevel, All, Add );

   //       real_t integValElem  = O.dotGlobal( resElementwise, maxLevel, All );
   //       real_t integValConst = O.dotGlobal( resConst, maxLevel, All );
   //       real_t integValCtrl  = real_c( -1.381773290675 );

   //       WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
   //       WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
   //       WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
   //       WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );

   //       divXConst->swapConstantScalingSign();
   //       divXElem->swapConstantScalingSign();
   //       divYConst->swapConstantScalingSign();
   //       divYElem->swapConstantScalingSign();
   //       divZConst->swapConstantScalingSign();
   //       divZElem->swapConstantScalingSign();
   //       divXConst->apply( u, resConst, maxLevel, All, Add );
   //       divXElem->apply( u, resElementwise, maxLevel, All, Add );
   //       divYConst->apply( u, resConst, maxLevel, All, Add );
   //       divYElem->apply( u, resElementwise, maxLevel, All, Add );
   //       divZConst->apply( u, resConst, maxLevel, All, Add );
   //       divZElem->apply( u, resElementwise, maxLevel, All, Add );

   //       real_t normElem  = resElementwise.dotGlobal( resElementwise, maxLevel, All );
   //       real_t normConst = resConst.dotGlobal( resConst, maxLevel, All );
   //       real_t normCTRL  = real_c( 0 );

   //       WALBERLA_LOG_INFO( abs( normElem - normCTRL ) );
   //       WALBERLA_CHECK_LESS( abs( normElem - normCTRL ), real_c( 1e-16 ) );
   //       WALBERLA_LOG_INFO( abs( normConst - normCTRL ) );
   //       WALBERLA_CHECK_LESS( abs( normConst - normCTRL ), real_c( 1e-16 ) );

   //       divXConst->setConstantScaling( real_c( 2 ) );
   //       divXElem->setConstantScaling( real_c( 2 ) );
   //       divYConst->setConstantScaling( real_c( 2 ) );
   //       divYElem->setConstantScaling( real_c( 2 ) );
   //       divZConst->setConstantScaling( real_c( 2 ) );
   //       divZElem->setConstantScaling( real_c( 2 ) );
   //       divXConst->apply( u, resConst, maxLevel, All, Replace );
   //       divXElem->apply( u, resElementwise, maxLevel, All, Replace );
   //       divYConst->apply( u, resConst, maxLevel, All, Add );
   //       divYElem->apply( u, resElementwise, maxLevel, All, Add );
   //       divZConst->apply( u, resConst, maxLevel, All, Add );
   //       divZElem->apply( u, resElementwise, maxLevel, All, Add );

   //       integValElem  = O.dotGlobal( resElementwise, maxLevel, All );
   //       integValConst = O.dotGlobal( resConst, maxLevel, All );
   //       integValCtrl  = real_c( -2.763546581350 );

   //       WALBERLA_LOG_INFO( abs( integValElem - integValCtrl ) );
   //       WALBERLA_CHECK_LESS( abs( integValElem - integValCtrl ), tolerance );
   //       WALBERLA_LOG_INFO( abs( integValConst - integValCtrl ) );
   //       WALBERLA_CHECK_LESS( abs( integValConst - integValCtrl ), tolerance );
   //    }

   WALBERLA_LOG_INFO( "Finish" );

   return EXIT_SUCCESS;
}