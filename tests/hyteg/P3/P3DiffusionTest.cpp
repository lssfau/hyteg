/*
 * Copyright (c) 2026 Marcus Mohr.
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

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Constants.h"

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P3ElementwiseOperator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p3functionspace/P3Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg {

void printSeparator()
{
   WALBERLA_LOG_INFO_ON_ROOT( "-----------------------------------------------" );
}

void testEigenfunctionsOnUnitSquare( bool doVTKOutput = false )
{
   hyteg::printSeparator();
   WALBERLA_LOG_INFO_ON_ROOT( "Testing Eigenfunctions on Unit Square" );

   Point2D lowerLeft( real_c( 0.0 ), real_c( 0.0 ) );
   Point2D upperRight( real_c( 1.0 ), real_c( 1.0 ) );

   MeshInfo              mesh2D = MeshInfo::meshRectangle( lowerLeft, upperRight, MeshInfo::meshFlavour::CRISS, 1, 1 );
   SetupPrimitiveStorage setupStorage2D( mesh2D, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage2D = std::make_shared< PrimitiveStorage >( setupStorage2D );

   const uint_t level = 6;

   P3Function< real_t > eigenfunction( "Eigenfunction", storage2D, level, level );
   P3Function< real_t > rhsFEM( "weak RHS", storage2D, level, level );
   P3Function< real_t > lhsOpr( "LHS with operator", storage2D, level, level );
   P3Function< real_t > error( "Error", storage2D, level, level );

   real_t freqX = real_c( 0 );
   real_t freqY = real_c( 0 );

   using walberla::math::pi;

   std::function< real_t( const hyteg::Point3D& ) > eigenExpr = [&freqX, &freqY]( const Point3D& x ) -> real_t {
      real_t retVal = std::sin( freqX * pi * x[0] ) * std::sin( freqY * pi * x[1] );
      return retVal;
   };

   // we will need a mass operator for the rhs and a diffusion operator
   P3ElementwiseMassOperator      massOper( storage2D, level, level );
   P3ElementwiseDiffusionOperator diffusion( storage2D, level, level );

   for ( uint_t k = 1; k <= 3; ++k )
   {
      for ( uint_t m = 1; m <= 3; ++m )
      {
         freqX = real_c( k );
         freqY = real_c( m );
         WALBERLA_LOG_INFO_ON_ROOT( "* (freqX, freqY) = (" << freqX << ", " << freqY << ")" );

         real_t eigenvalue = ( freqX * freqX + freqY * freqY ) * pi * pi;
         WALBERLA_LOG_INFO_ON_ROOT( "  -> analytical eigenvalue = " << eigenvalue );

         eigenfunction.interpolate( eigenExpr, level, All );
         massOper.apply( eigenfunction, rhsFEM, level, Inner );
         diffusion.apply( eigenfunction, lhsOpr, level, Inner );

         error.assign( { real_c( -1 ), eigenvalue }, { lhsOpr, rhsFEM }, level, All );
         real_t errorMeasure = error.getMaxDoFMagnitude( level, All );
         WALBERLA_LOG_INFO_ON_ROOT( "  -> error (max dof magnitude) = " << errorMeasure );

         if ( doVTKOutput )
         {
            VTKOutput vtkOutput( ".", "P3DiffusionTest", storage2D );
            vtkOutput.add( eigenfunction );
            vtkOutput.add( rhsFEM );
            vtkOutput.add( lhsOpr );
            vtkOutput.add( error );
            vtkOutput.write( level );
         }

         real_t tolerance = std::is_same_v< real_t, double > ? real_c( 5e-6 ) : real_c( 5e-5 );
         WALBERLA_CHECK_LESS_EQUAL( errorMeasure, tolerance );
      }
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::INFO );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::testEigenfunctionsOnUnitSquare();

   hyteg::printSeparator();
   return EXIT_SUCCESS;
}
