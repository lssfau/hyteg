#include "tinyhhg_core/p1functionspace/P1PointwiseOperator.hpp"

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/timing/Timer.h"

#include "tinyhhg_core/FunctionIterator.hpp"
#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/generated/p1_tet_diffusion.h"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"
#include "tinyhhg_core/types/pointnd.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hhg;

static void testPointwiseOperator( const std::string& meshFile, const uint_t& minLevel, const uint_t& maxLevel )
{
   const bool   writeVTK   = true;
   const real_t errorLimit = 1e-13;

   const auto            meshInfo = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   const auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   WALBERLA_CHECK( !storage->hasGlobalCells() );

   if ( writeVTK )
      writeDomainPartitioningVTK( storage, "../../output", "P1PointwiseOperatorTest_Domain" );

   P1PointwiseOperator pointwiseOperator( storage, minLevel, maxLevel );

   // We fill all stencil entries of the operator at each point with values
   // that depend on the coordinates of the respective DoF.
   // Thus we can interpolate some function by multiplying the operator with the 1-vector

   auto someFunction = []( const Point3D& p ) -> real_t { return std::sin( p[0] ) + p[1]; };

   P1Function< real_t > u( "u", storage, minLevel, maxLevel );
   P1Function< real_t > Au( "Au", storage, minLevel, maxLevel );
   P1Function< real_t > AuExpected( "AuExpected", storage, minLevel, maxLevel );
   P1Function< real_t > error( "error", storage, minLevel, maxLevel );

   VTKOutput vtkOutput( "../../output", "P1PointwiseOperatorTest", storage );
   vtkOutput.add( u );
   vtkOutput.add( Au );
   vtkOutput.add( AuExpected );
   vtkOutput.add( error );

   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      // Fill operator:
      // 1. Looping over all DoFs of the domain
      // 2. Check which primitive the DoF is located on and obtain the pointer to the stencil memory
      // 3. Set each stencil weight to
      //      func( dof.coords ) / numStencilWeights

      for ( const auto& dof : FunctionIterator< P1Function< real_t > >( u, level ) )
      {
         const auto id = dof.primitiveID();

         if ( dof.isOnMacroVertex() )
         {
            auto       macroVertex       = storage->getVertex( id );
            const auto numStencilWeights = macroVertex->getNumNeighborEdges() + 1;
            auto       stencilPtr        = macroVertex->getData( pointwiseOperator.getVertexStencilID() )->getPointer( level );
            const auto stencilEntry      = someFunction( dof.coordinates() ) / real_c( numStencilWeights );
            for ( uint_t i = 0; i < numStencilWeights; i++ )
            {
               stencilPtr[i] = stencilEntry;
            }
         }
         else if ( dof.isOnMacroEdge() )
         {
            auto       macroEdge         = storage->getEdge( id );
            const auto numStencilWeights = 3 + 2 * macroEdge->getNumNeighborFaces();
            auto       stencilPtr        = macroEdge->getData( pointwiseOperator.getEdgeStencilID() )->getPointer( level );
            const auto stencilEntry      = someFunction( dof.coordinates() ) / real_c( numStencilWeights );
            for ( uint_t i = 0; i < numStencilWeights; i++ )
            {
               stencilPtr[numStencilWeights * vertexdof::macroedge::innerIndex( level, dof.index().x() ) + i] = stencilEntry;
            }
         }
         else if ( dof.isOnMacroFace() )
         {
            auto       macroFace         = storage->getFace( id );
            const auto numStencilWeights = 7;
            auto       stencilPtr        = macroFace->getData( pointwiseOperator.getFaceStencilID() )->getPointer( level );
            const auto stencilEntry      = someFunction( dof.coordinates() ) / real_c( numStencilWeights );
            for ( uint_t i = 0; i < numStencilWeights; i++ )
            {
               stencilPtr[numStencilWeights * vertexdof::macroface::innerIndex( level, dof.index().x(), dof.index().y() ) + i] =
                   stencilEntry;
            }
         }
      }

      u.interpolate( 1, level );
      pointwiseOperator.apply( u, Au, level, All );
      AuExpected.interpolate( someFunction, level );
      error.assign( {1.0, -1.0}, {Au, AuExpected}, level );

      if ( writeVTK )
         vtkOutput.write( level, 0 );

      const auto errorNorm = std::sqrt( error.dotGlobal( error, level ) );
      WALBERLA_CHECK_LESS( errorNorm, errorLimit );
   }
}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   testPointwiseOperator( "../../data/meshes/tri_1el.msh", 2, 4 );
   testPointwiseOperator( "../../data/meshes/quad_8el.msh", 2, 4 );
   testPointwiseOperator( "../../data/meshes/annulus_coarse.msh", 2, 4 );

   return 0;
}
