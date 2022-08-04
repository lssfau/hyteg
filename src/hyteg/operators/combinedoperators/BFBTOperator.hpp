

#include <type_traits>

//#include "hyteg/forms/form_hyteg_generated/p2/p2_invk_mass_affine_q6.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.cpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_invk_mass_affine_q6.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_k_mass_affine_q4.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/p2functionspace/P2EpsilonOperator.hpp"
#include "hyteg/p2functionspace/P2FunctionApplyOperator.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/solvers/CGSolver.hpp"

namespace hyteg {

using walberla::real_t;
template < class SaddlePointOpType, class VelocityBlockOpType, class DivOpType, class DivOpTType >
class BFBTOperator : public Operator< typename DivOpType::dstType, typename DivOpType::dstType >
{
   // sub operator consistency
   static_assert( std::is_same< typename VelocityBlockOpType::srcType, typename VelocityBlockOpType::dstType >::value == true );
   static_assert( std::is_same< typename VelocityBlockOpType::srcType, typename DivOpType::srcType >::value == true );
   static_assert( std::is_same< typename DivOpType::dstType, typename DivOpTType::srcType >::value == true );
   using PressureFunctionType    = typename DivOpType::dstType;
   using VelocityFunctionType    = typename DivOpType::srcType;
   using VectorVelocityMassOp    = VectorToVectorOperator< walberla::real_t, P2VectorFunction, P2VectorFunction >;
   using SaddlePointFunctionType = typename SaddlePointOpType::srcType;

 public:
   BFBTOperator( const std::shared_ptr< PrimitiveStorage >&       storage,
                 size_t                                           minLevel,
                 size_t                                           maxLevel,
                 std::function< real_t( const hyteg::Point3D& ) > viscosity,
                 const std::vector< BoundaryCondition >&          VelocitySpaceBCs )
   : Operator< PressureFunctionType, PressureFunctionType >( storage, maxLevel, maxLevel )
   , A( storage, maxLevel, maxLevel, viscosity ) // A should be a discretization of the velocity block with non constant viscosity
   , Z( storage, maxLevel, viscosity, VelocitySpaceBCs )
   , ZSolver( std::make_shared< PETScLUSolver< BFBTSubOperator > >( storage, maxLevel ) )
   , storage_( storage )
   , pTmp( "pTmp", storage, maxLevel, maxLevel )
   , vTmp( "vTmp", storage, maxLevel, maxLevel )
   , vTmp2( "vTmp2", storage, maxLevel, maxLevel )
   {
      A.computeInverseDiagonalOperatorValues();
      Z.Minv.diagonal_ = A.getInverseDiagonalValues();

      for ( uint_t c = 0; c < VelocitySpaceBCs.size(); c++ )
      {
         vTmp.setBoundaryCondition( VelocitySpaceBCs[c], c );
         vTmp2.setBoundaryCondition( VelocitySpaceBCs[c], c );
      }
      // test for Minv
      /*
      // Tmp funcs
      VelocityFunctionType tmpf1( "tmp0", storage_, maxLevel, maxLevel );
      VelocityFunctionType tmpf2( "tmp1", storage_, maxLevel, maxLevel );
      
      // mult Minv with 1 
      tmpf1.interpolate( []( const hyteg::Point3D& ) { return real_c( 1 ); }, maxLevel, hyteg::All );
      Z.Minv.apply(tmpf1, tmpf2, maxLevel, hyteg::All);
      P2FunctionApplyOperator tmpOp( storage, maxLevel );
      tmpOp.diagonal_ = std::make_shared<VelocityFunctionType>(tmpf2);
      
      // print
      typename VelocityFunctionType::template FunctionType< idx_t > VelocityNum(
          "VelocityNum", storage , maxLevel, maxLevel );
      VelocityNum.enumerate( maxLevel );

      PETScSparseMatrix< P2FunctionApplyOperator > tmpMat;
      tmpMat.createMatrixFromOperator( tmpOp, maxLevel, VelocityNum );
      tmpMat.print( "tmpMat", false, PETSC_VIEWER_DEFAULT );
      printComponentMatrices(maxLevel, storage);
      WALBERLA_ABORT("bye");
    */

      //ZSolver->setPrintInfo( true );
      //printComponentMatrices(maxLevel, storage);
   }

   void apply( const PressureFunctionType& src,
               const PressureFunctionType& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace )
   {
      //S_{w-BFBT}^{-1} = Z^{-1}*(BM(w)^{-1}AM(w)^{-1}B^t)*Z^{-1}
      pTmp.copyBoundaryConditionFromFunction( src );
      pTmp.interpolate( []( const hyteg::Point3D& ) { return real_c( 0 ); }, level, hyteg::Inner );
      vTmp.interpolate( []( const hyteg::Point3D& ) { return real_c( 0 ); }, level, hyteg::Inner );
      vTmp2.interpolate( []( const hyteg::Point3D& ) { return real_c( 0 ); }, level, hyteg::Inner );

      VTKOutput vtkOutput_pTmp(
          "../../output", "BFBT_intres_BFBTApply_pTmp", Operator< PressureFunctionType, PressureFunctionType >::storage_ );
      VTKOutput vtkOutput_vTmp(
          "../../output", "BFBT_intres_BFBTApply_vTmp", Operator< PressureFunctionType, PressureFunctionType >::storage_ );
      VTKOutput vtkOutput_vTmp2(
          "../../output", "BFBT_intres_BFBTApply_vTmp2", Operator< PressureFunctionType, PressureFunctionType >::storage_ );
      vtkOutput_pTmp.add( pTmp );
      vtkOutput_pTmp.add( dst );

      vtkOutput_vTmp.add( vTmp );
      vtkOutput_vTmp2.add( vTmp2 );
      vtkOutput_pTmp.write( level, 0 );

      ZSolver->solve( Z, pTmp, src, level ); // pTmp = Z^{-1}*src

      //vtkOutput_pTmp.write( level, 1 );

      Z.BT.apply( pTmp, vTmp, level, hyteg::Inner ); // vTmp = B^t*pTmp

      //vtkOutput_vTmp.write( level, 0 );

      Z.Minv.apply( vTmp, vTmp2, level, hyteg::Inner ); // vTmp2 = M^{-1}*vTmp
      //Z.MSolver->solve( Z.M, vTmp2, vTmp, level ); // vTmp2= M(w)^{-1}*vTmp

      //vtkOutput_vTmp2.write( level, 0 );
      A.apply( vTmp2, vTmp, level, hyteg::Inner ); // vtmp= A*vTmp2
      //vtkOutput_vTmp.write( level, 1 );

      Z.Minv.apply( vTmp, vTmp2, level, hyteg::Inner ); // vTmp2 = M^{-1}*vTmp
      //Z.MSolver->solve( Z.M, vTmp2, vTmp, level ); // vTmp2= M(w)^{-1}*vTmp

      //vtkOutput_vTmp2.write( level, 1 );

      Z.B.apply( vTmp2, pTmp, level, hyteg::Inner ); // pTmp = B*vTmp2
      //vtkOutput_pTmp.write( level, 2 );

      ZSolver->solve( Z, dst, pTmp, level ); // dst = Z^{-1}*pTmp

      //vtkOutput_pTmp.write( level, 3 );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >&                          mat,
                  const typename PressureFunctionType::template FunctionType< idx_t >& src,
                  const typename PressureFunctionType::template FunctionType< idx_t >& dst,
                  size_t                                                               level,
                  DoFType                                                              flag ) const
   {
      WALBERLA_ABORT( "Requires exact inversion of A, should not be used!" );
   }

   void printComponentMatrices( size_t level, std::shared_ptr< PrimitiveStorage > storage )
   {
      typename VelocityFunctionType::template FunctionType< idx_t > VelocityNum(
          "VelocityNum", storage /*Operator< PressureFunctionType, PressureFunctionType >::storage_*/, level, level );

      typename PressureFunctionType::template FunctionType< idx_t > PressureNum(
          "PressureNum", Operator< PressureFunctionType, PressureFunctionType >::storage_, level, level );

      VelocityNum.enumerate( level );
      PressureNum.enumerate( level );

      PETScSparseMatrix< VectorVelocityMassOp >               MMat;
      PETScSparseMatrix< DivOpType >                          BMat;
      PETScSparseMatrix< DivOpTType >                         BTMat;
      PETScSparseMatrix< P2ElementwiseAffineEpsilonOperator > AMat;

      MMat.createMatrixFromOperator( Z.Minv, level, VelocityNum );
      AMat.createMatrixFromOperator( A, level, VelocityNum );
      BMat.createMatrixFromOperator( Z.B, level, VelocityNum, PressureNum );
      BTMat.createMatrixFromOperator( Z.BT, level, PressureNum, VelocityNum );
      MMat.applyDirichletBCSymmetrically( VelocityNum, level );
      AMat.applyDirichletBCSymmetrically( VelocityNum, level );
      BMat.applyDirichletBCSymmetrically( VelocityNum, level );
      BTMat.applyDirichletBCSymmetrically( PressureNum, level );

      MMat.print( "FOR_MATLAB_MinvMat.m", false, PETSC_VIEWER_ASCII_MATLAB );
      AMat.print( "FOR_MATLAB_AMat.m", false, PETSC_VIEWER_ASCII_MATLAB );
      BMat.print( "FOR_MATLAB_BMat.m", false, PETSC_VIEWER_ASCII_MATLAB );
      BTMat.print( "FOR_MATLAB_BTMat.m", false, PETSC_VIEWER_ASCII_MATLAB );
   }

   // this implements the suboperator Z = BM(w)^{-1}B^T
   class BFBTSubOperator : public Operator< PressureFunctionType, PressureFunctionType >
   {
    public:
      BFBTSubOperator( const std::shared_ptr< PrimitiveStorage >&       storage,
                       size_t                                           level,
                       std::function< real_t( const hyteg::Point3D& ) > viscosity,
                       const std::vector< BoundaryCondition >&          VelocitySpaceBCs
                       //std::shared_ptr<P2VectorFunction< real_t >> diagonal
                       )
      : Operator< PressureFunctionType, PressureFunctionType >( storage, level, level )
      , B( storage, level, level )
      , BT( storage, level, level )
      , Minv( storage, level )
      , K( storage, level, level, viscosity )
      //, M( storage, level, level )
      //, MSolver( std::make_shared< PETScLUSolver< VectorVelocityMassOp > >( storage, level ) )
      , storage_( storage )
      , tmp( "tmp", storage, level, level )
      , tmp2( "tmp2", storage, level, level )
      {
         for ( uint_t c = 0; c < VelocitySpaceBCs.size(); c++ )
         {
            tmp.setBoundaryCondition( VelocitySpaceBCs[c], c );
            tmp2.setBoundaryCondition( VelocitySpaceBCs[c], c );
         }

         //using massForm = hyteg::forms::p2_sqrtk_mass_affine_q4;
         //auto P2Mass =
         //   std::make_shared< P2ElementwiseOperator< massForm > >( storage, level, level, massForm( viscosity, viscosity ) );
         //M.setSubOperator( 0, 0, P2Mass );
         //M.setSubOperator( 1, 1, P2Mass );
      }

      void apply( const PressureFunctionType& src,
                  const PressureFunctionType& dst,
                  size_t                      level,
                  DoFType                     flag,
                  UpdateType                  updateType = Replace ) const
      {
         // tmp0.interpolate( []( const hyteg::Point3D& ) { return real_c( 0 ); }, level, hyteg::All );
         //  tmp1.interpolate( []( const hyteg::Point3D& ) { return real_c( 0 ); }, level, hyteg::All );
         //tmp0.interpolate( []( const hyteg::Point3D& ) { return real_c( 0 ); }, level, hyteg::DirichletBoundary );
         // tmp1.interpolate( []( const hyteg::Point3D& ) { return real_c( 0 ); }, level, hyteg::DirichletBoundary );

         //tmp0.setBoundaryCondition(BoundaryCondition(hyteg::Inner));
         //tmp1.setBoundaryCondition(BoundaryCondition(hyteg::Inner));

         //WALBERLA_LOG_INFO_ON_ROOT( "tmp0 bc type for 1: " << tmp0.getBoundaryCondition().getBoundaryType( 1 ) );

         //WALBERLA_LOG_INFO_ON_ROOT( "tmp0 bc type for 0: " << tmp0.getBoundaryCondition().getBoundaryType( 0 ) );
         BT.apply( src, tmp, level, hyteg::All ); // tmp = B^T*src

         Minv.apply( tmp, tmp2, level, hyteg::All ); // tmp2 = M^{-1}*tmp
         //MSolver->solve( M, tmp2, tmp, level );
         B.apply( tmp2, dst, level, hyteg::All ); // dst = B*tmp2
      }
      void catchPetscError( PetscErrorCode PetscErr ){};

      void toMatrix( const std::shared_ptr< SparseMatrixProxy >&                          mat,
                     const typename PressureFunctionType::template FunctionType< idx_t >& src,
                     const typename PressureFunctionType::template FunctionType< idx_t >& dst,
                     size_t                                                               level,
                     DoFType                                                              flag ) const
      {
         // Z = BM(w)^{-1}B^T to matrix

       //  if ( !assembledOnce_ )
       //  {
            // Sparse matrix objects for sub operators
            PETScSparseMatrix< P2FunctionApplyOperator > MinvMat( "MinvMat" );

            // numerators
            typename VelocityFunctionType::template FunctionType< idx_t >    VelocityNum( "VelocityNum", storage_, level, level );
            typename PressureFunctionType::template FunctionType< idx_t >    PressureNum( "PressureNum", storage_, level, level );
            typename SaddlePointFunctionType::template FunctionType< idx_t > SPNum( "SPNum", storage_, level, level );
            SPNum.enumerate( level );
            VelocityNum.enumerate( level );
            PressureNum.enumerate( level );

            // assemble Minv
            // MinvMat.createMatrixFromOperator( Minv, level, VelocityNum );
            // MinvMat.applyDirichletBCSymmetrically( VelocityNum, level );

            // assemble K, apply boundary conditions
            PETScSparseMatrix< SaddlePointOpType > KMat;
            KMat.createMatrixFromOperator( K, level, SPNum, All );
            KMat.applyDirichletBCSymmetrically( SPNum, level );

            // gather index sets for B, BT
            uint_t   nPDoFs   = numberOfGlobalDoFs< P1FunctionTag >( *storage_, level );
            uint_t   nDoFs    = numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage_, level );
            uint_t   nVelDoFs = nDoFs - nPDoFs;
            PetscInt velRange[nVelDoFs];
            for ( PetscInt i = 0; i < nVelDoFs; ++i )
            {
               velRange[i] = i;
            }
            PetscInt pRange[nPDoFs];
            for ( PetscInt i = nVelDoFs; i < nDoFs; ++i )
            {
               pRange[i - nVelDoFs] = i;
            }
            IS             velIndices, pIndices;
            PetscErrorCode PetscErr;
            PetscErr = ISCreateGeneral( PETSC_COMM_WORLD, nVelDoFs, velRange, PETSC_COPY_VALUES, &velIndices );
            WALBERLA_CHECK( !PetscErr );
            PetscErr = ISCreateGeneral( PETSC_COMM_WORLD, nPDoFs, pRange, PETSC_COPY_VALUES, &pIndices );
            WALBERLA_CHECK( !PetscErr );
            IS BBTSrcIndices[2] = { velIndices, pIndices };
            IS BBTDstIndices[2] = { pIndices, velIndices };

            // get B and BT from K
            Mat BMat, BTMat, AMat;
            PetscErr = MatCreateSubMatrix( KMat.get(), velIndices, pIndices, MAT_INITIAL_MATRIX, &BTMat );
            WALBERLA_CHECK( !PetscErr );
            PetscErr = MatCreateSubMatrix( KMat.get(), pIndices, velIndices, MAT_INITIAL_MATRIX, &BMat );
            WALBERLA_CHECK( !PetscErr );

            // get A
            PetscErr = MatCreateSubMatrix( KMat.get(), velIndices, velIndices, MAT_INITIAL_MATRIX, &AMat );
            WALBERLA_CHECK( !PetscErr );
            PetscErr = MatAssemblyBegin( AMat, MAT_FINAL_ASSEMBLY );
            WALBERLA_CHECK( !PetscErr );
            PetscErr = MatAssemblyEnd( AMat, MAT_FINAL_ASSEMBLY );
            WALBERLA_CHECK( !PetscErr );
            // MatView(AMat, PETSC_VIEWER_STDOUT_WORLD);

            // extract diagonal
            Vec diagonal;
            PetscErr = VecCreateSeq( PETSC_COMM_WORLD, nVelDoFs, &diagonal );
            WALBERLA_CHECK( !PetscErr );
            std::vector< PetscScalar > vals( nVelDoFs );
            PetscErr = MatGetDiagonal( AMat, diagonal );
            WALBERLA_CHECK( !PetscErr );
            PetscErr = VecGetValues( diagonal, nVelDoFs, velRange, vals.data() );
            WALBERLA_CHECK( !PetscErr );
            //MatView(diagAMat, PETSC_VIEWER_STDOUT_WORLD);
            PetscErr = MatAssemblyBegin( AMat, MAT_FINAL_ASSEMBLY );
            WALBERLA_CHECK( !PetscErr );
            PetscErr = MatAssemblyEnd( AMat, MAT_FINAL_ASSEMBLY );
            WALBERLA_CHECK( !PetscErr );

            // invert diagonal and put back into vector
            for ( PetscInt i = 0; i < nVelDoFs; ++i )
            {
               vals[i] = 1.0 / vals[i];
            }
            PetscErr = VecSetValues( diagonal, nVelDoFs, velRange, vals.data(), INSERT_VALUES );
            WALBERLA_CHECK( !PetscErr );
            // MatView(diagAMat, PETSC_VIEWER_STDOUT_WORLD);

            // insert inverted diagonal
            Mat diagAMat;
            PetscErr = MatCreate( PETSC_COMM_WORLD, &diagAMat );
            WALBERLA_CHECK( !PetscErr );
            PetscErr = MatSetType( diagAMat, MATMPIAIJ );
            WALBERLA_CHECK( !PetscErr );
            PetscErr = MatSetSizes( diagAMat, nVelDoFs, nVelDoFs, nVelDoFs, nVelDoFs );
            WALBERLA_CHECK( !PetscErr );
            PetscErr = MatMPIAIJSetPreallocation( diagAMat, 1, NULL, 0, NULL );
            WALBERLA_CHECK( !PetscErr );
            PetscErr = MatDiagonalSet( diagAMat, diagonal, INSERT_VALUES );
            WALBERLA_CHECK( !PetscErr );
         
         
            PetscErr = MatAssemblyBegin( diagAMat, MAT_FINAL_ASSEMBLY );
            WALBERLA_CHECK( !PetscErr );
            PetscErr = MatAssemblyEnd( diagAMat, MAT_FINAL_ASSEMBLY );
            WALBERLA_CHECK( !PetscErr );
            //MatCreateSubMatrices( KMat.get(), 2, BBTDstIndices, BBTSrcIndices, MAT_INITIAL_MATRIX, BBTPtrs );

            // PetscViewer viewer;
            /*
         PetscViewer viewer;
         PetscViewerASCIIOpen( PETSC_COMM_WORLD, "FOR_MATLAB_invMMat_Sinker.m", &viewer );
         PetscViewerPushFormat( viewer, PETSC_VIEWER_ASCII_MATLAB );
         MatView(diagAMat,  viewer );
         
         PetscViewerDestroy( &viewer );
         PetscViewerASCIIOpen( PETSC_COMM_WORLD, "FOR_MATLAB_BMat_Sinker.m", &viewer );
         PetscViewerPushFormat( viewer, PETSC_VIEWER_ASCII_MATLAB );
         MatView(BMat,  viewer );
         */

            // Proxy objects for createFromMatrixProduct
            auto MinvMatProxy = std::make_shared< PETScSparseMatrixProxy >( diagAMat );
            auto BMatProxy    = std::make_shared< PETScSparseMatrixProxy >( BMat );
            auto BTMatProxy   = std::make_shared< PETScSparseMatrixProxy >( BTMat );

            // final BMinvBT matrix creation
            mat->createFromMatrixProduct( { BTMatProxy, MinvMatProxy, BMatProxy } );
        // }
        // else {}
         /*
         PetscViewerDestroy( &viewer );
           PetscViewerASCIIOpen( PETSC_COMM_WORLD, "FOR_MATLAB_BinvMBTMat_Sinker.m", &viewer );
         PetscViewerPushFormat( viewer, PETSC_VIEWER_ASCII_MATLAB );
         MatView(std::dynamic_pointer_cast< PETScSparseMatrixProxy >( mat )->mat_,  viewer );
      */
      }

      SaddlePointOpType       K;
      DivOpType               B;
      DivOpTType              BT;
      P2FunctionApplyOperator Minv;
      //VectorVelocityMassOp                                     M;
      //std::shared_ptr< PETScLUSolver< VectorVelocityMassOp > > MSolver;
      VelocityFunctionType                      tmp;
      VelocityFunctionType                      tmp2;
      const std::shared_ptr< PrimitiveStorage > storage_;
   };

 private:
   VelocityBlockOpType A;
   BFBTSubOperator     Z;

   std::shared_ptr< PETScLUSolver< BFBTSubOperator > > ZSolver;

   const std::shared_ptr< PrimitiveStorage > storage_;
   PressureFunctionType                      pTmp;
   VelocityFunctionType                      vTmp;
   VelocityFunctionType                      vTmp2;
};

//TODO ToMatrix function not ready for other template specializations
using BFBT_P2P1 = BFBTOperator< P2P1ElementwiseAffineEpsilonStokesOperator,
                                P2ElementwiseAffineEpsilonOperator,
                                P2ToP1ElementwiseDivOperator,
                                P1ToP2ElementwiseDivTOperator >;

} // namespace hyteg
