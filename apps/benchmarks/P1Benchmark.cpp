#include "core/Environment.h"

#include "tinyhhg_core/LikwidWrapper.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_c;
using walberla::real_t;
using namespace hhg;

void kernel( double* fd_p1FaceDst, double* fd_p1FaceSrc, double* fd_p1FaceStencil, int64_t level )
{
   for( int ctr_2 = 1; ctr_2 < ( 1 << level ) - 1; ctr_2 += 1 )
      for( int ctr_1 = 1; ctr_1 < -ctr_2 + ( 1 << level ); ctr_1 += 1 )
      {
         fd_p1FaceDst[ctr_1 + ctr_2 * ( ( 1 << level ) + 2 ) - ( ctr_2 * ( ctr_2 + 1 ) / 2 )] =
             fd_p1FaceSrc[ctr_1 + ( ctr_2 + 1 ) * ( ( 1 << level ) + 2 ) - ( ( ctr_2 + 1 ) * ( ctr_2 + 2 ) / 2 ) - 1] *
                 fd_p1FaceStencil[5] +
             fd_p1FaceSrc[ctr_1 + ( ctr_2 + 1 ) * ( ( 1 << level ) + 2 ) - ( ( ctr_2 + 1 ) * ( ctr_2 + 2 ) / 2 )] *
                 fd_p1FaceStencil[6] +
             fd_p1FaceSrc[ctr_1 + ( ctr_2 - 1 ) * ( ( 1 << level ) + 2 ) - ( ctr_2 * ( ctr_2 - 1 ) / 2 ) + 1] *
                 fd_p1FaceStencil[1] +
             fd_p1FaceSrc[ctr_1 + ( ctr_2 - 1 ) * ( ( 1 << level ) + 2 ) - ( ctr_2 * ( ctr_2 - 1 ) / 2 )] * fd_p1FaceStencil[0] +
             fd_p1FaceSrc[ctr_1 + ctr_2 * ( ( 1 << level ) + 2 ) - ( ctr_2 * ( ctr_2 + 1 ) / 2 ) + 1] * fd_p1FaceStencil[4] +
             fd_p1FaceSrc[ctr_1 + ctr_2 * ( ( 1 << level ) + 2 ) - ( ctr_2 * ( ctr_2 + 1 ) / 2 ) - 1] * fd_p1FaceStencil[2] +
             fd_p1FaceSrc[ctr_1 + ctr_2 * ( ( 1 << level ) + 2 ) - ( ctr_2 * ( ctr_2 + 1 ) / 2 )] * fd_p1FaceStencil[3];
      }
}

int main( int argc, char** argv )
{
   LIKWID_MARKER_INIT;

   walberla::debug::enterTestMode();
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   LIKWID_MARKER_THREADINIT;

   MeshInfo                            meshInfo = MeshInfo::fromGmshFile( "../data/meshes/tri_1el.msh" );
   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const size_t level = 10;

   auto src = std::make_shared< hhg::P1Function< real_t > >( "src", storage, level, level );
   auto dst = std::make_shared< hhg::P1Function< real_t > >( "dst", storage, level, level );

   hhg::P1MassOperator M( storage, level, level );

   std::shared_ptr< Face > face = storage->getFaces().begin().operator*().second;

   std::function< real_t( const hhg::Point3D& ) > exactFunc = [&]( const hhg::Point3D& point ) {
      return sqrt( point[0] * point[0] + point[1] * point[1] );
   };

   //P1Function< real_t > x("x", storage, level, level);
   src->interpolate( exactFunc, level );

   walberla::WcTimer timer;

   LIKWID_MARKER_START( "apply" );
   timer.reset();
   vertexdof::macroface::apply< real_t >(
       level, *face, M.getFaceStencilID(), src->getFaceDataID(), dst->getFaceDataID(), Replace );
   timer.end();
   LIKWID_MARKER_STOP( "apply" );
   WALBERLA_LOG_INFO_ON_ROOT( "apply time: " << timer.last() );

   real_t check1 = vertexdof::macroface::dot< real_t >( level, *face, dst->getFaceDataID(), dst->getFaceDataID() );

   LIKWID_MARKER_START( "assign" );
   timer.reset();
   vertexdof::macroface::assign< real_t >(
         level, *face, {13}, {src->getFaceDataID()}, dst->getFaceDataID() );
   timer.end();
   LIKWID_MARKER_STOP( "assign" );
   WALBERLA_LOG_INFO_ON_ROOT( "assign timer: " << timer.last() );


   real_t* opr_data = face->getData(M.getFaceStencilID())->getPointer( level );
   real_t* srcPtr = face->getData(src->getFaceDataID())->getPointer( level );
   real_t* dstPtr = face->getData(dst->getFaceDataID())->getPointer( level );

   LIKWID_MARKER_START( "apply_gen" );
   timer.reset();
   kernel(dstPtr,srcPtr,opr_data,level);
   timer.end();
   LIKWID_MARKER_STOP( "apply_gen" );
   WALBERLA_LOG_INFO_ON_ROOT( "time with walberla timer: " << timer.last() );

   /// do something with the result to prevent the compiler from removing all the computations
   real_t check2 = vertexdof::macroface::dot< real_t >( level, *face, dst->getFaceDataID(), dst->getFaceDataID() );
   WALBERLA_LOG_INFO_ON_ROOT(check2);
   WALBERLA_CHECK_FLOAT_EQUAL( check1 , check2 );

   LIKWID_MARKER_CLOSE;
}
