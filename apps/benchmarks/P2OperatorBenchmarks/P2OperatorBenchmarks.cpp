#include <iostream>
#include <vector>

#include "core/Environment.h"
#include "core/timing/Timer.h"

#include "tinyhhg_core/LikwidWrapper.hpp"
#include "tinyhhg_core/edgedofspace/generatedKernels/all.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/misc/dummy.hpp"
#include "tinyhhg_core/mixedoperators/EdgeDoFToVertexDoFOperator/generatedKernels/all.hpp"
#include "tinyhhg_core/mixedoperators/VertexDoFToEdgeDoFOperator/generatedKernels/all.hpp"
#include "tinyhhg_core/p1functionspace/generatedKernels/all.hpp"
#include "tinyhhg_core/p2functionspace/P2ConstantOperator.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/p2functionspace/generatedKernels/all.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using walberla::uint_t;

namespace hhg {

enum KernelType
{
   APPLY_V_TO_V_REPLACE,
   APPLY_V_TO_V_ADD,

   APPLY_V_TO_E_REPLACE,
   APPLY_V_TO_E_ADD,

   APPLY_E_TO_V_REPLACE,
   APPLY_E_TO_V_ADD,

   APPLY_E_TO_E_REPLACE,
   APPLY_E_TO_E_ADD,

   SOR_P1_V,
   SOR_P2_V,
   SOR_P2_E
};

const std::map< std::string, KernelType > strToKernelType = {
    {"APPLY_V_TO_V_REPLACE", APPLY_V_TO_V_REPLACE},
    {"APPLY_V_TO_V_ADD", APPLY_V_TO_V_ADD},

    {"APPLY_V_TO_E_REPLACE", APPLY_V_TO_E_REPLACE},
    {"APPLY_V_TO_E_ADD", APPLY_V_TO_E_ADD},

    {"APPLY_E_TO_V_REPLACE", APPLY_E_TO_V_REPLACE},
    {"APPLY_E_TO_V_ADD", APPLY_E_TO_V_ADD},

    {"APPLY_E_TO_E_REPLACE", APPLY_E_TO_E_REPLACE},
    {"APPLY_E_TO_E_ADD", APPLY_E_TO_E_ADD},

    {"SOR_P1_V", SOR_P1_V},
    {"SOR_P2_V", SOR_P2_V},
    {"SOR_P2_E", SOR_P2_E},
};

const std::map< KernelType, std::string > kernelTypeToString = {
    {APPLY_V_TO_V_REPLACE, "APPLY_V_TO_V_REPLACE"},
    {APPLY_V_TO_V_ADD, "APPLY_V_TO_V_ADD"},

    {APPLY_V_TO_E_REPLACE, "APPLY_V_TO_E_REPLACE"},
    {APPLY_V_TO_E_ADD, "APPLY_V_TO_E_ADD"},

    {APPLY_E_TO_V_REPLACE, "APPLY_E_TO_V_REPLACE"},
    {APPLY_E_TO_V_ADD, "APPLY_E_TO_V_ADD"},

    {APPLY_E_TO_E_REPLACE, "APPLY_E_TO_E_REPLACE"},
    {APPLY_E_TO_E_ADD, "APPLY_E_TO_E_ADD"},

    {SOR_P1_V, "SOR_P1_V"},
    {SOR_P2_V, "SOR_P2_V"},
    {SOR_P2_E, "SOR_P2_E"},
};

void printArithmeticIntensity( KernelType kernelType )
{
   // loadsSrcAll: all loads from src vector (excluding write allocate / load+store of dst vector) if there was no caching
   // loadsSrcMin: only the furthest entry must be loaded from each array (in regions where everything fits into the caches)
   // loadsSrcMax: all fronts
   uint_t adds( 0 ), mults( 0 ), loadsSrcAll( 0 ), loadsSrcMin( 0 ), loadsSrcMax( 0 ), loadsDst( 0 ), writesDst( 0 );
   switch ( kernelType )
   {
   case APPLY_V_TO_V_REPLACE:
      adds        = 14;
      mults       = 15;
      loadsSrcAll = 15;
      loadsSrcMin = 1;
      loadsDst    = 1;
      writesDst   = 1;
      break;
   case APPLY_V_TO_V_ADD:
      adds        = 15;
      mults       = 15;
      loadsSrcAll = 15;
      loadsSrcMin = 1;
      loadsDst    = 1;
      writesDst   = 1;
      break;

   case APPLY_E_TO_V_REPLACE:
      adds        = 49;
      mults       = 50;
      loadsSrcAll = 50;
      loadsSrcMin = 7;
      loadsDst    = 1;
      writesDst   = 1;
      break;
   case APPLY_E_TO_V_ADD:
      adds        = 49 + 1;
      mults       = 50;
      loadsSrcAll = 50;
      loadsSrcMax = 31;
      loadsSrcMin = 7;
      loadsDst    = 1;
      writesDst   = 1;
      break;

   case APPLY_V_TO_E_REPLACE:
      adds        = 7 + 5 + 7 + 7 + 5 + 7 + 5;
      mults       = 8 + 6 + 8 + 8 + 6 + 8 + 6;
      loadsSrcAll = mults;
      loadsSrcMin = 7;
      loadsDst    = 7;
      writesDst   = 7;
      break;
   case APPLY_V_TO_E_ADD:
      adds        = 7 + 5 + 7 + 7 + 5 + 7 + 5 + 7;
      mults       = 8 + 6 + 8 + 8 + 6 + 8 + 6;
      loadsSrcAll = mults;
      loadsSrcMin = 7;
      loadsDst    = 7;
      writesDst   = 7;
      break;

   case APPLY_E_TO_E_REPLACE:
      adds        = 18 + 12 + 18 + 18 + 12 + 18 + 12;
      mults       = 19 + 13 + 19 + 19 + 13 + 19 + 13;
      loadsSrcAll = mults;
      loadsSrcMin = 7;
      loadsDst    = 7;
      writesDst   = 7;
      break;

   case APPLY_E_TO_E_ADD:
      adds        = 18 + 12 + 18 + 18 + 12 + 18 + 12 + 7;
      mults       = 19 + 13 + 19 + 19 + 13 + 19 + 13;
      loadsSrcAll = mults;
      loadsSrcMin = 7;
      loadsDst    = 7;
      writesDst   = 7;
      break;

   case SOR_P1_V:
      adds        = 13 + 1;
      mults       = 14 + 1 + 1;
      loadsSrcAll = 16;
      loadsSrcMin = 2;
      loadsDst    = 0;
      writesDst   = 1;
      break;

   case SOR_P2_V:
      adds        = 63 + 1;
      mults       = 14 + 50 + 1 + 1;
      loadsSrcAll = 66;
      loadsSrcMin = 9;
      loadsDst    = 0;
      writesDst   = 1;
      break;

   case SOR_P2_E:
      adds        = 25 * 4 + 17 * 3 + 7;
      mults       = 26 * 4 + 18 * 3 + 7 + 7;
      loadsSrcAll = 27 * 4 + 19 * 3 + 7;
      loadsSrcMin = 8;
      loadsDst    = 0;
      writesDst   = 1;
      break;
   }

   real_t allIntensity = real_c( adds + mults ) / real_c( 8 * ( loadsDst + writesDst + loadsSrcAll ) );
   real_t minIntensity = real_c( adds + mults ) / real_c( 8 * ( loadsDst + writesDst + loadsSrcMax ) );
   real_t maxIntensity = real_c( adds + mults ) / real_c( 8 * ( loadsDst + writesDst + loadsSrcMin ) );

   WALBERLA_UNUSED( minIntensity );

   WALBERLA_LOG_INFO_ON_ROOT( "Arithmetic intensity:" )
   WALBERLA_LOG_INFO_ON_ROOT( "  - adds:            " << adds )
   WALBERLA_LOG_INFO_ON_ROOT( "  - mults:           " << mults )
   WALBERLA_LOG_INFO_ON_ROOT( "  - FLOPs:           " << adds + mults )
   WALBERLA_LOG_INFO_ON_ROOT( "  - loads all (src): " << loadsSrcAll )
   WALBERLA_LOG_INFO_ON_ROOT( "  - loads min (src): " << loadsSrcMin )
   WALBERLA_LOG_INFO_ON_ROOT( "  - loads (dst):     " << loadsDst )
   WALBERLA_LOG_INFO_ON_ROOT( "  - writes (dst):    " << writesDst )
   WALBERLA_LOG_INFO_ON_ROOT( "  -----------------------------------" )
   WALBERLA_LOG_INFO_ON_ROOT( "  - intensity no cache:                                            " << allIntensity )
   WALBERLA_LOG_INFO_ON_ROOT( "  - intensity large cache / perfect blocking / top of tetrahedron: " << maxIntensity )
}

void runBenchmark( uint_t level, real_t iterationMinTime, uint_t chunkSize, KernelType kernelType )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Running benchmark " << kernelTypeToString.at( kernelType ) )
   WALBERLA_LOG_INFO_ON_ROOT( "  - level: " << level )

   printArithmeticIntensity( kernelType );

   walberla::WcTimer timer;

   const MeshInfo              meshInfo = MeshInfo::meshSymmetricCuboid( Point3D( {0, 0, 0} ), Point3D( {1, 1, 1} ), 1, 1, 1 );
   const SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   const auto                  storage = std::make_shared< PrimitiveStorage >( setupStorage );

   WALBERLA_CHECK_GREATER_EQUAL( setupStorage.getNumberOfCells(),
                                 uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   P2Function< real_t > p2Src( "x", storage, level, level );
   P2Function< real_t > p2Dst( "x", storage, level, level );

   std::function< real_t( const hhg::Point3D& ) > someFunction = []( const hhg::Point3D& x ) {
      return cos( M_PI * x[0] ) - sin( 2.0 * M_PI * x[1] ) + cos( 2.0 * M_PI * x[2] );
   };

   WALBERLA_LOG_INFO_ON_ROOT( "alsdfkn" )
   p2Src.interpolate( someFunction, level );
   p2Dst.interpolate( someFunction, level );
   WALBERLA_LOG_INFO_ON_ROOT( "sssss" )

   P2ConstantLaplaceOperator p2Operator( storage, level, level );

   // get any tet from the storage
   const Cell& cell = *storage->getCell( storage->getCellIDs().front() );

   // gather all the pointers
   auto v2v_opr_data = cell.getData( p2Operator.getVertexToVertexOpr().getCellStencilID() )->getData( level );
   auto v2v_src_data = cell.getData( p2Src.getVertexDoFFunction().getCellDataID() )->getPointer( level );
   auto v2v_dst_data = cell.getData( p2Dst.getVertexDoFFunction().getCellDataID() )->getPointer( level );

   auto v2e_opr_data = cell.getData( p2Operator.getVertexToEdgeOpr().getCellStencilID() )->getData( level );
   auto v2e_src_data = cell.getData( p2Src.getVertexDoFFunction().getCellDataID() )->getPointer( level );
   auto v2e_dst_data = cell.getData( p2Dst.getEdgeDoFFunction().getCellDataID() )->getPointer( level );

   auto e2v_opr_data = cell.getData( p2Operator.getEdgeToVertexOpr().getCellStencilID() )->getData( level );
   auto e2v_src_data = cell.getData( p2Src.getEdgeDoFFunction().getCellDataID() )->getPointer( level );
   auto e2v_dst_data = cell.getData( p2Dst.getVertexDoFFunction().getCellDataID() )->getPointer( level );

   auto e2e_opr_data = cell.getData( p2Operator.getEdgeToEdgeOpr().getCellStencilID() )->getData( level );
   auto e2e_src_data = cell.getData( p2Src.getEdgeDoFFunction().getCellDataID() )->getPointer( level );
   auto e2e_dst_data = cell.getData( p2Dst.getEdgeDoFFunction().getCellDataID() )->getPointer( level );

   auto v_dst_data = cell.getData( p2Dst.getVertexDoFFunction().getCellDataID() )->getPointer( level );
   auto e_dst_data = cell.getData( p2Dst.getEdgeDoFFunction().getCellDataID() )->getPointer( level );

   auto v_rhs_data = cell.getData( p2Src.getVertexDoFFunction().getCellDataID() )->getPointer( level );
   auto e_rhs_data = cell.getData( p2Src.getEdgeDoFFunction().getCellDataID() )->getPointer( level );

   typedef edgedof::EdgeDoFOrientation eo;
   std::map< eo, uint_t >              firstEdgeIdx;
   for ( auto e : edgedof::allEdgeDoFOrientations )
      firstEdgeIdx[e] = edgedof::macrocell::index( level, 0, 0, 0, e );

   /////////////////////
   // Apply benchmark //
   /////////////////////

   uint_t chunk = 0;
   while ( timer.total() < iterationMinTime )
   {
      timer.start();
      LIKWID_MARKER_START( "OperatorBenchmark" );
      switch ( kernelType )
      {
      case APPLY_V_TO_V_REPLACE:

         for ( uint_t i = 0; i < chunkSize; i++ )
         {
            vertexdof::macrocell::generated::apply_3D_macrocell_vertexdof_to_vertexdof_replace(
                v2v_dst_data, v2v_src_data, static_cast< int32_t >( level ), v2v_opr_data );
            hhg::misc::dummy( v2v_src_data );
            hhg::misc::dummy( v2v_dst_data );
         }
         break;

      case APPLY_V_TO_V_ADD:
         for ( uint_t i = 0; i < chunkSize; i++ )
         {
            vertexdof::macrocell::generated::apply_3D_macrocell_vertexdof_to_vertexdof_add(
                v2v_dst_data, v2v_src_data, static_cast< int32_t >( level ), v2v_opr_data );
            hhg::misc::dummy( v2v_src_data );
            hhg::misc::dummy( v2v_dst_data );
         }
         break;

      case APPLY_E_TO_V_REPLACE:
         for ( uint_t i = 0; i < chunkSize; i++ )
         {
            EdgeDoFToVertexDoF::generated::apply_3D_macrocell_edgedof_to_vertexdof_replace( &e2v_src_data[firstEdgeIdx[eo::X]],
                                                                                            &e2v_src_data[firstEdgeIdx[eo::XY]],
                                                                                            &e2v_src_data[firstEdgeIdx[eo::XYZ]],
                                                                                            &e2v_src_data[firstEdgeIdx[eo::XZ]],
                                                                                            &e2v_src_data[firstEdgeIdx[eo::Y]],
                                                                                            &e2v_src_data[firstEdgeIdx[eo::YZ]],
                                                                                            &e2v_src_data[firstEdgeIdx[eo::Z]],
                                                                                            e2v_dst_data,
                                                                                            e2v_opr_data,
                                                                                            static_cast< int32_t >( level ) );
            hhg::misc::dummy( e2v_src_data );
            hhg::misc::dummy( e2v_dst_data );
         }
         break;

      case APPLY_E_TO_V_ADD:
         for ( uint_t i = 0; i < chunkSize; i++ )
         {
            EdgeDoFToVertexDoF::generated::apply_3D_macrocell_edgedof_to_vertexdof_add( &e2v_src_data[firstEdgeIdx[eo::X]],
                                                                                        &e2v_src_data[firstEdgeIdx[eo::XY]],
                                                                                        &e2v_src_data[firstEdgeIdx[eo::XYZ]],
                                                                                        &e2v_src_data[firstEdgeIdx[eo::XZ]],
                                                                                        &e2v_src_data[firstEdgeIdx[eo::Y]],
                                                                                        &e2v_src_data[firstEdgeIdx[eo::YZ]],
                                                                                        &e2v_src_data[firstEdgeIdx[eo::Z]],
                                                                                        e2v_dst_data,
                                                                                        e2v_opr_data,
                                                                                        static_cast< int32_t >( level ) );
            hhg::misc::dummy( e2v_src_data );
            hhg::misc::dummy( e2v_dst_data );
         }
         break;

      case APPLY_E_TO_E_REPLACE:
         for ( uint_t i = 0; i < chunkSize; i++ )
         {
            edgedof::macrocell::generated::apply_3D_macrocell_edgedof_to_edgedof_replace( &e2e_dst_data[firstEdgeIdx[eo::X]],
                                                                                          &e2e_dst_data[firstEdgeIdx[eo::XY]],
                                                                                          &e2e_dst_data[firstEdgeIdx[eo::XYZ]],
                                                                                          &e2e_dst_data[firstEdgeIdx[eo::XZ]],
                                                                                          &e2e_dst_data[firstEdgeIdx[eo::Y]],
                                                                                          &e2e_dst_data[firstEdgeIdx[eo::YZ]],
                                                                                          &e2e_dst_data[firstEdgeIdx[eo::Z]],
                                                                                          &e2e_src_data[firstEdgeIdx[eo::X]],
                                                                                          &e2e_src_data[firstEdgeIdx[eo::XY]],
                                                                                          &e2e_src_data[firstEdgeIdx[eo::XYZ]],
                                                                                          &e2e_src_data[firstEdgeIdx[eo::XZ]],
                                                                                          &e2e_src_data[firstEdgeIdx[eo::Y]],
                                                                                          &e2e_src_data[firstEdgeIdx[eo::YZ]],
                                                                                          &e2e_src_data[firstEdgeIdx[eo::Z]],
                                                                                          e2e_opr_data,
                                                                                          static_cast< int32_t >( level ) );
            hhg::misc::dummy( e2e_src_data );
            hhg::misc::dummy( e2e_dst_data );
         }
         break;
      case APPLY_E_TO_E_ADD:
         for ( uint_t i = 0; i < chunkSize; i++ )
         {
            edgedof::macrocell::generated::apply_3D_macrocell_edgedof_to_edgedof_add( &e2e_dst_data[firstEdgeIdx[eo::X]],
                                                                                      &e2e_dst_data[firstEdgeIdx[eo::XY]],
                                                                                      &e2e_dst_data[firstEdgeIdx[eo::XYZ]],
                                                                                      &e2e_dst_data[firstEdgeIdx[eo::XZ]],
                                                                                      &e2e_dst_data[firstEdgeIdx[eo::Y]],
                                                                                      &e2e_dst_data[firstEdgeIdx[eo::YZ]],
                                                                                      &e2e_dst_data[firstEdgeIdx[eo::Z]],
                                                                                      &e2e_src_data[firstEdgeIdx[eo::X]],
                                                                                      &e2e_src_data[firstEdgeIdx[eo::XY]],
                                                                                      &e2e_src_data[firstEdgeIdx[eo::XYZ]],
                                                                                      &e2e_src_data[firstEdgeIdx[eo::XZ]],
                                                                                      &e2e_src_data[firstEdgeIdx[eo::Y]],
                                                                                      &e2e_src_data[firstEdgeIdx[eo::YZ]],
                                                                                      &e2e_src_data[firstEdgeIdx[eo::Z]],
                                                                                      e2e_opr_data,
                                                                                      static_cast< int32_t >( level ) );
            hhg::misc::dummy( e2e_src_data );
            hhg::misc::dummy( e2e_dst_data );
         }
         break;
      case APPLY_V_TO_E_REPLACE:
         for ( uint_t i = 0; i < chunkSize; i++ )
         {
            VertexDoFToEdgeDoF::generated::apply_3D_macrocell_vertexdof_to_edgedof_replace( &v2e_dst_data[firstEdgeIdx[eo::X]],
                                                                                            &v2e_dst_data[firstEdgeIdx[eo::XY]],
                                                                                            &v2e_dst_data[firstEdgeIdx[eo::XYZ]],
                                                                                            &v2e_dst_data[firstEdgeIdx[eo::XZ]],
                                                                                            &v2e_dst_data[firstEdgeIdx[eo::Y]],
                                                                                            &v2e_dst_data[firstEdgeIdx[eo::YZ]],
                                                                                            &v2e_dst_data[firstEdgeIdx[eo::Z]],
                                                                                            v2e_src_data,
                                                                                            static_cast< int32_t >( level ),
                                                                                            v2e_opr_data );
            hhg::misc::dummy( v2e_src_data );
            hhg::misc::dummy( v2e_dst_data );
         }
         break;
      case APPLY_V_TO_E_ADD:
         for ( uint_t i = 0; i < chunkSize; i++ )
         {
            VertexDoFToEdgeDoF::generated::apply_3D_macrocell_vertexdof_to_edgedof_add( &v2e_dst_data[firstEdgeIdx[eo::X]],
                                                                                        &v2e_dst_data[firstEdgeIdx[eo::XY]],
                                                                                        &v2e_dst_data[firstEdgeIdx[eo::XYZ]],
                                                                                        &v2e_dst_data[firstEdgeIdx[eo::XZ]],
                                                                                        &v2e_dst_data[firstEdgeIdx[eo::Y]],
                                                                                        &v2e_dst_data[firstEdgeIdx[eo::YZ]],
                                                                                        &v2e_dst_data[firstEdgeIdx[eo::Z]],
                                                                                        v2e_src_data,
                                                                                        static_cast< int32_t >( level ),
                                                                                        v2e_opr_data );
            hhg::misc::dummy( v2e_src_data );
            hhg::misc::dummy( v2e_dst_data );
         }
         break;

      case SOR_P1_V:
         for ( uint_t i = 0; i < chunkSize; i++ )
         {
            vertexdof::macrocell::generated::sor_3D_macrocell_P1(
                v_dst_data, v_rhs_data, static_cast< int32_t >( level ), v2v_opr_data, 1.1 );
            hhg::misc::dummy( v_rhs_data );
            hhg::misc::dummy( v_dst_data );
         }
         break;

      case SOR_P2_V:
         for ( uint_t i = 0; i < chunkSize; i++ )
         {
            P2::macrocell::generated::sor_3D_macrocell_P2_update_vertexdofs( &e_dst_data[firstEdgeIdx[eo::X]],
                                                                             &e_dst_data[firstEdgeIdx[eo::XY]],
                                                                             &e_dst_data[firstEdgeIdx[eo::XYZ]],
                                                                             &e_dst_data[firstEdgeIdx[eo::XZ]],
                                                                             &e_dst_data[firstEdgeIdx[eo::Y]],
                                                                             &e_dst_data[firstEdgeIdx[eo::YZ]],
                                                                             &e_dst_data[firstEdgeIdx[eo::Z]],
                                                                             v_dst_data,
                                                                             v_rhs_data,
                                                                             e2v_opr_data,
                                                                             static_cast< int32_t >( level ),
                                                                             1.1,
                                                                             v2v_opr_data );
            hhg::misc::dummy( e_dst_data );
            hhg::misc::dummy( v_dst_data );
            hhg::misc::dummy( v_rhs_data );
         }
         break;

      case SOR_P2_E:
         for ( uint_t i = 0; i < chunkSize; i++ )
         {
           const real_t relax = 1.1;
           P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_X( &e_dst_data[firstEdgeIdx[eo::X]],
                                                                                    &e_dst_data[firstEdgeIdx[eo::XY]],
                                                                                    &e_dst_data[firstEdgeIdx[eo::XYZ]],
                                                                                    &e_dst_data[firstEdgeIdx[eo::XZ]],
                                                                                    &e_dst_data[firstEdgeIdx[eo::Y]],
                                                                                    &e_dst_data[firstEdgeIdx[eo::YZ]],
                                                                                    &e_dst_data[firstEdgeIdx[eo::Z]],
                                                                                    &e_rhs_data[firstEdgeIdx[eo::X]],
                                                                                    v_dst_data,
                                                                                    e2e_opr_data,
                                                                                    static_cast< int32_t >( level ),
                                                                                    relax,
                                                                                    v2e_opr_data );
           P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_Y( &e_dst_data[firstEdgeIdx[eo::X]],
                                                                                    &e_dst_data[firstEdgeIdx[eo::XY]],
                                                                                    &e_dst_data[firstEdgeIdx[eo::XYZ]],
                                                                                    &e_dst_data[firstEdgeIdx[eo::XZ]],
                                                                                    &e_dst_data[firstEdgeIdx[eo::Y]],
                                                                                    &e_dst_data[firstEdgeIdx[eo::YZ]],
                                                                                    &e_dst_data[firstEdgeIdx[eo::Z]],
                                                                                    &e_rhs_data[firstEdgeIdx[eo::Y]],
                                                                                    v_dst_data,
                                                                                    e2e_opr_data,
                                                                                    static_cast< int32_t >( level ),
                                                                                    relax,
                                                                                    v2e_opr_data );
           P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_Z( &e_dst_data[firstEdgeIdx[eo::X]],
                                                                                    &e_dst_data[firstEdgeIdx[eo::XY]],
                                                                                    &e_dst_data[firstEdgeIdx[eo::XYZ]],
                                                                                    &e_dst_data[firstEdgeIdx[eo::XZ]],
                                                                                    &e_dst_data[firstEdgeIdx[eo::Y]],
                                                                                    &e_dst_data[firstEdgeIdx[eo::YZ]],
                                                                                    &e_dst_data[firstEdgeIdx[eo::Z]],
                                                                                    &e_rhs_data[firstEdgeIdx[eo::Z]],
                                                                                    v_dst_data,
                                                                                    e2e_opr_data,
                                                                                    static_cast< int32_t >( level ),
                                                                                    relax,
                                                                                    v2e_opr_data );
           P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_XY( &e_dst_data[firstEdgeIdx[eo::X]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::XY]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::XYZ]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::XZ]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::Y]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::YZ]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::Z]],
                                                                                     &e_rhs_data[firstEdgeIdx[eo::XY]],
                                                                                     v_dst_data,
                                                                                     e2e_opr_data,
                                                                                     static_cast< int32_t >( level ),
                                                                                     relax,
                                                                                     v2e_opr_data );
           P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_XZ( &e_dst_data[firstEdgeIdx[eo::X]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::XY]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::XYZ]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::XZ]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::Y]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::YZ]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::Z]],
                                                                                     &e_rhs_data[firstEdgeIdx[eo::XZ]],
                                                                                     v_dst_data,
                                                                                     e2e_opr_data,
                                                                                     static_cast< int32_t >( level ),
                                                                                     relax,
                                                                                     v2e_opr_data );
           P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_YZ( &e_dst_data[firstEdgeIdx[eo::X]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::XY]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::XYZ]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::XZ]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::Y]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::YZ]],
                                                                                     &e_dst_data[firstEdgeIdx[eo::Z]],
                                                                                     &e_rhs_data[firstEdgeIdx[eo::YZ]],
                                                                                     v_dst_data,
                                                                                     e2e_opr_data,
                                                                                     static_cast< int32_t >( level ),
                                                                                     relax,
                                                                                     v2e_opr_data );
           P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_XYZ( &e_dst_data[firstEdgeIdx[eo::X]],
                                                                                      &e_dst_data[firstEdgeIdx[eo::XY]],
                                                                                      &e_dst_data[firstEdgeIdx[eo::XYZ]],
                                                                                      &e_dst_data[firstEdgeIdx[eo::XZ]],
                                                                                      &e_dst_data[firstEdgeIdx[eo::Y]],
                                                                                      &e_dst_data[firstEdgeIdx[eo::YZ]],
                                                                                      &e_dst_data[firstEdgeIdx[eo::Z]],
                                                                                      &e_rhs_data[firstEdgeIdx[eo::XYZ]],
                                                                                      v_dst_data,
                                                                                      e2e_opr_data,
                                                                                      static_cast< int32_t >( level ),
                                                                                      relax,
                                                                                      v2e_opr_data );
            hhg::misc::dummy( e_dst_data );
            hhg::misc::dummy( v_dst_data );
            hhg::misc::dummy( e_rhs_data );
         }
         break;

      default:
         WALBERLA_ABORT( "Kernel type not implemented." );
         break;
      }
      LIKWID_MARKER_STOP( "OperatorBenchmark" );
      timer.end();

      chunk++;
      // WALBERLA_LOG_INFO_ON_ROOT( "Iteration: " << chunk * chunkSize << ": " << timer.last() << "sec" );
   }

   LIKWID_MARKER_CLOSE;

   WALBERLA_LOG_INFO_ON_ROOT( "Total number of iterations: " << chunk * chunkSize << ", total time: " << timer.total() << "sec" );
}
} // namespace hhg
int main( int argc, char** argv )
{
   LIKWID_MARKER_INIT;

   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   LIKWID_MARKER_THREADINIT;

   LIKWID_MARKER_REGISTER( "OperatorBenchmark" );

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./P2OperatorBenchmarks.prm";
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf         = cfg->getBlock( "Parameters" );
   const uint_t                        level            = mainConf.getParameter< uint_t >( "level" );
   const real_t                        iterationMinTime = mainConf.getParameter< real_t >( "iterationMinTime" );
   const uint_t                        chunkSize        = mainConf.getParameter< uint_t >( "chunkSize" );
   const std::string                   kernelType       = mainConf.getParameter< std::string >( "kernelType" );

   hhg::runBenchmark( level, iterationMinTime, chunkSize, hhg::strToKernelType.at( kernelType ) );
}