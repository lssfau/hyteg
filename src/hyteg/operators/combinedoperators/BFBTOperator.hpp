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

   // types for component functions
   using PressureFunctionType    = typename DivOpType::dstType;
   using VelocityFunctionType    = typename DivOpType::srcType;
   using VectorVelocityMassOp    = VectorToVectorOperator< walberla::real_t, P2VectorFunction, P2VectorFunction >;
   using SaddlePointFunctionType = typename SaddlePointOpType::srcType;

   // types for numerators
   using VelocityNumeratorType    = typename VelocityFunctionType::template FunctionType< idx_t >;
   using PressureNumeratorType    = typename PressureFunctionType::template FunctionType< idx_t >;
   using SaddlePointNumeratorType = typename SaddlePointFunctionType::template FunctionType< idx_t >;

 public:
   BFBTOperator( const std::shared_ptr< PrimitiveStorage >&       storage,
                 size_t                                           minLevel,
                 size_t                                           maxLevel,
                 std::function< real_t( const hyteg::Point3D& ) > viscosity,
                 const std::vector< BoundaryCondition >&          VelocitySpaceBCs )
   : Operator< PressureFunctionType, PressureFunctionType >( storage, maxLevel, maxLevel )
   , A( storage, maxLevel, maxLevel, viscosity ) // A should be a discretization of the velocity block with non constant viscosity
   , Z( storage, maxLevel, viscosity, VelocitySpaceBCs )
   , ZSolver( std::make_shared< CGSolver< BFBTSubOperator > >( storage, maxLevel, maxLevel ) )
   , storage_( storage )
   , pTmp( "pTmp", storage, maxLevel, maxLevel )
   , vTmp( "vTmp", storage, maxLevel, maxLevel )
   , vTmp2( "vTmp2", storage, maxLevel, maxLevel )
   {
      //ZSolver->setPrintInfo( true );
      ZSolver->setDoFType(hyteg::All);

      // get BCs for artificial velocity space used in BFBT
      for ( uint_t c = 0; c < VelocitySpaceBCs.size(); c++ )
      {
         vTmp.setBoundaryCondition( VelocitySpaceBCs[c], c );
         vTmp2.setBoundaryCondition( VelocitySpaceBCs[c], c );
      }
   }

   void apply( const PressureFunctionType& src,
               const PressureFunctionType& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace )
   {
      //S_{w-BFBT}^{-1} = Z^{-1}*(BM(w)^{-1}AM(w)^{-1}B^t)*Z^{-1}
      pTmp.copyBoundaryConditionFromFunction( src );
      /*vTmp.interpolate( 0, level, hyteg::Inner );
      pTmp.interpolate( 0, level, hyteg::Inner );
      vTmp2.interpolate( 0, level, hyteg::Inner );
      // random initial guess
        std::uniform_real_distribution< real_t > unif( -1, 1 );
        std::default_random_engine               re;
        re.seed( 1312412415 );
        pTmp.interpolate( [&unif, &re]( const hyteg::Point3D& xx ) { return unif( re ); }, level, hyteg::Inner );*/

      ZSolver->solve( Z, pTmp, src, level ); // pTmp = Z^{-1}*src
                                             //  pTmp.interpolate( 0, level, hyteg::DirichletBoundary );
                                             //  hyteg::vertexdof::projectMean( pTmp, level );

      Z.BT.apply( pTmp, vTmp, level, hyteg::All ); // vTmp = B^t*pTmp
      vTmp.interpolate( 0, level, hyteg::DirichletBoundary );

      vTmp2.multElementwise( { *( Z.invDiagA ), vTmp }, level, hyteg::All );
      vTmp2.interpolate( 0, level, hyteg::DirichletBoundary );

      A.apply( vTmp2, vTmp, level, hyteg::All ); // vtmp= A*vTmp2
      vTmp.interpolate( 0, level, hyteg::DirichletBoundary );

      vTmp2.multElementwise( { *( Z.invDiagA ), vTmp }, level, hyteg::All );
      vTmp2.interpolate( 0, level, hyteg::DirichletBoundary );

      Z.B.apply( vTmp2, pTmp, level, hyteg::All ); // pTmp = B*vTmp2
                                                   // pTmp.interpolate( 0, level, hyteg::DirichletBoundary );

      // dst.interpolate( [&unif, &re]( const hyteg::Point3D& xx ) { return unif( re ); }, level, hyteg::Inner );

      ZSolver->solve( Z, dst, pTmp, level ); // dst = Z^{-1}*pTmp
                                             //  dst.interpolate( 0, level, hyteg::DirichletBoundary );
                                             // hyteg::vertexdof::projectMean( dst, level );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >&                          mat,
                  const typename PressureFunctionType::template FunctionType< idx_t >& src,
                  const typename PressureFunctionType::template FunctionType< idx_t >& dst,
                  size_t                                                               level,
                  DoFType                                                              flag ) const
   {
      WALBERLA_ABORT( "Requires exact inversion of A, should not be used!" );
   }

   // this implements the suboperator Z = BM(w)^{-1}B^T
   class BFBTSubOperator : public Operator< PressureFunctionType, PressureFunctionType >
   {
    public:
      BFBTSubOperator( const std::shared_ptr< PrimitiveStorage >&       storage,
                       size_t                                           level,
                       std::function< real_t( const hyteg::Point3D& ) > viscosity,
                       const std::vector< BoundaryCondition >&          VelocitySpaceBCs )
      : Operator< PressureFunctionType, PressureFunctionType >( storage, level, level )
      , B( storage, level, level )
      , BT( storage, level, level )
      , Minv( storage, level )
      , K( storage, level, level, viscosity )
      , storage_( storage )
      , tmp( "tmp", storage, level, level )
      , tmp2( "tmp2", storage, level, level )
      , tmp3( "tmp3", storage, level, level )
      , level_( level )
      , invDiagAAssembled( false )
      {
         // get BCs for temporaries
         for ( uint_t c = 0; c < VelocitySpaceBCs.size(); c++ )
         {
            tmp.setBoundaryCondition( VelocitySpaceBCs[c], c );
            tmp2.setBoundaryCondition( VelocitySpaceBCs[c], c );
         }
         tmp3.setBoundaryCondition( VelocitySpaceBCs[0] );

         // obtain inverse diagonal of viscous block
         auto [SPNum, VelocityNum, PressureNum]                          = getNumerators();
         auto [nVelDoFs, velRange, velIndices, nPDoFs, pRange, pIndices] = getIndexSets();
         auto                 invDiagonalVec   = getVelocityBlockDiagonal( nVelDoFs, velRange, VelocityNum );
         auto                 invDiagonalProxy = std::make_shared< PETScVectorProxy >( invDiagonalVec );
         VelocityFunctionType invDiagonalFunction( "Z_invDiagVec_Hyteg", storage_, level_, level_ );
         invDiagonalFunction.fromVector( VelocityNum, invDiagonalProxy, level_, hyteg::All );
         //  PETScVector invDiagVec( invDiagonalFunction, VelocityNum, level, All, "Z_invDiagVec_Hyteg" );
         //  invDiagVec.print( "FOR_MATLAB_DEBUG_Z_invDiagVec.m", false, PETSC_VIEWER_ASCII_MATLAB );
         tmp2.multElementwise( { invDiagonalFunction, tmp }, level, hyteg::All );
         invDiagA = std::make_shared< VelocityFunctionType >( invDiagonalFunction );
      }

      void apply( const PressureFunctionType& src,
                  const PressureFunctionType& dst,
                  size_t                      level,
                  DoFType                     flag,
                  UpdateType                  updateType = Replace ) const
      {
         auto [SPNum, VelocityNum, PressureNum] = getNumerators();
         bool print                             = false;
         if ( print )
         {
            PETScVector srcVec( src, PressureNum, level, All, "Z_srcVec_Hyteg" );
            srcVec.print( "FOR_MATLAB_DEBUG_Z_srcVec.m", false, PETSC_VIEWER_ASCII_MATLAB );
         }

         BT.apply( src, tmp, level, hyteg::All ); // tmp = B^T*src
         tmp.interpolate( 0, level, hyteg::DirichletBoundary );

         tmp2.multElementwise( { *invDiagA, tmp }, level, hyteg::All );
         tmp2.interpolate( 0, level, hyteg::DirichletBoundary );

         if ( print )
         {
            PETScVector tmp2Vec( tmp2, VelocityNum, level, All, "Z_tmp2Vec_Hyteg" );
            tmp2Vec.print( "FOR_MATLAB_DEBUG_Z_tmp2Vec.m", false, PETSC_VIEWER_ASCII_MATLAB );
         }

         B.apply( tmp2, tmp3, level, hyteg::All ); // dst = B*tmp
         //Why does this destroy everything?
         // tmp3.interpolate( 0, level, hyteg::DirichletBoundary );

         if ( print )
         {
            PETScVector tmp3Vec( tmp3, PressureNum, level, All, "Z_tmp3Vec_Hyteg" );
            tmp3Vec.print( "FOR_MATLAB_DEBUG_Z_tmp3Vec.m", false, PETSC_VIEWER_ASCII_MATLAB );
            WALBERLA_ABORT( "BYE" )
         }

         dst.assign( { 1 }, { tmp3 }, level, hyteg::All );
      }

      auto getNumerators() const
      {
         VelocityNumeratorType    VelocityNum( "VelocityNum", storage_, level_, level_ );
         PressureNumeratorType    PressureNum( "PressureNum", storage_, level_, level_ );
         SaddlePointNumeratorType SPNum( "SPNum", storage_, level_, level_ );
         SPNum.enumerate( level_ );
         VelocityNum.enumerate( level_ );
         PressureNum.enumerate( level_ );
         return std::make_tuple( std::move( SPNum ), std::move( VelocityNum ), std::move( PressureNum ) );
      }

      auto getIndexSets() const
      {
         uint_t                  nPDoFs   = numberOfGlobalDoFs< P1FunctionTag >( *storage_, level_ );
         uint_t                  nDoFs    = numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage_, level_ );
         uint_t                  nVelDoFs = nDoFs - nPDoFs;
         std::vector< PetscInt > velRange( nVelDoFs );
         for ( PetscInt i = 0; i < nVelDoFs; ++i )
         {
            velRange[i] = i;
         }
         std::vector< PetscInt > pRange( nPDoFs );
         for ( PetscInt i = nVelDoFs; i < nDoFs; ++i )
         {
            pRange[i - nVelDoFs] = i;
         }
         IS             velIndices, pIndices;
         PetscErrorCode PetscErr;
         PetscErr = ISCreateGeneral( PETSC_COMM_WORLD, nVelDoFs, velRange.data(), PETSC_COPY_VALUES, &velIndices );
         WALBERLA_CHECK( !PetscErr );
         PetscErr = ISCreateGeneral( PETSC_COMM_WORLD, nPDoFs, pRange.data(), PETSC_COPY_VALUES, &pIndices );
         WALBERLA_CHECK( !PetscErr );
         return std::make_tuple(
             nVelDoFs, std::move( velRange ), std::move( velIndices ), nPDoFs, std::move( pRange ), std::move( pIndices ) );
      }

      auto
          getVelocityBlockDiagonal( uint_t nVelDoFs, std::vector< PetscInt >& velRange, VelocityNumeratorType& VelocityNum ) const
      {
         PetscErrorCode PetscErr;

         // assemble visous block
         PETScSparseMatrix< VelocityBlockOpType > AMat;
         AMat.createMatrixFromOperator( K.viscOp, level_, VelocityNum, All );
         AMat.applyDirichletBCSymmetrically( VelocityNum, level_ );

         // extract diagonal
         Vec diagonal;
         PetscErr = VecCreateSeq( PETSC_COMM_WORLD, nVelDoFs, &diagonal );
         WALBERLA_CHECK( !PetscErr );
         std::vector< PetscScalar > vals( nVelDoFs );
         PetscErr = MatGetDiagonal( AMat.get(), diagonal );
         WALBERLA_CHECK( !PetscErr );
         PetscErr = VecGetValues( diagonal, nVelDoFs, velRange.data(), vals.data() );
         WALBERLA_CHECK( !PetscErr );
         PetscErr = MatAssemblyBegin( AMat.get(), MAT_FINAL_ASSEMBLY );
         WALBERLA_CHECK( !PetscErr );
         PetscErr = MatAssemblyEnd( AMat.get(), MAT_FINAL_ASSEMBLY );
         WALBERLA_CHECK( !PetscErr );

         // invert diagonal and put back into vector
         for ( PetscInt i = 0; i < nVelDoFs; ++i )
         {
            vals[i] = 1.0 / vals[i];
         }
         PetscErr = VecSetValues( diagonal, nVelDoFs, velRange.data(), vals.data(), INSERT_VALUES );
         WALBERLA_CHECK( !PetscErr );

         return diagonal;
      }

      void toMatrix( const std::shared_ptr< SparseMatrixProxy >&                          mat,
                     const typename PressureFunctionType::template FunctionType< idx_t >& src,
                     const typename PressureFunctionType::template FunctionType< idx_t >& dst,
                     size_t                                                               level,
                     DoFType                                                              flag ) const
      {
         // Z = BM(w)^{-1}B^T to matrix
         // this assembly only happens once in the solver, so no need to have an assembledOnce mechanism

         // numerators
         auto [SPNum, VelocityNum, PressureNum] = getNumerators();

         PetscErrorCode PetscErr;
         // assemble K, apply boundary conditions
         PETScSparseMatrix< SaddlePointOpType > KMat;
         KMat.createMatrixFromOperator( K, level_, SPNum, All );
         KMat.applyDirichletBCSymmetrically( SPNum, level_ );

         // gather index sets for B, BT
         auto [nVelDoFs, velRange, velIndices, nPDoFs, pRange, pIndices] = getIndexSets();
         IS BBTSrcIndices[2]                                             = { velIndices, pIndices };
         IS BBTDstIndices[2]                                             = { pIndices, velIndices };

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

         // extract inverted diagonal
         Vec invDiagonal = getVelocityBlockDiagonal( nVelDoFs, velRange, VelocityNum );

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
         PetscErr = MatDiagonalSet( diagAMat, invDiagonal, INSERT_VALUES );
         WALBERLA_CHECK( !PetscErr );
         PetscErr = MatAssemblyBegin( diagAMat, MAT_FINAL_ASSEMBLY );
         WALBERLA_CHECK( !PetscErr );
         PetscErr = MatAssemblyEnd( diagAMat, MAT_FINAL_ASSEMBLY );
         WALBERLA_CHECK( !PetscErr );

         // Proxy objects for createFromMatrixProduct
         auto MinvMatProxy = std::make_shared< PETScSparseMatrixProxy >( diagAMat );
         auto BMatProxy    = std::make_shared< PETScSparseMatrixProxy >( BMat );
         auto BTMatProxy   = std::make_shared< PETScSparseMatrixProxy >( BTMat );

         // final BMinvBT matrix creation
         mat->createFromMatrixProduct( { BTMatProxy, MinvMatProxy, BMatProxy } );
      }

      SaddlePointOpType K;
      //TODO save K only once
      DivOpType               B;
      DivOpTType              BT;
      P2FunctionApplyOperator Minv;
      //VectorVelocityMassOp                                     M;
      //std::shared_ptr< PETScLUSolver< VectorVelocityMassOp > > MSolver;

      //TODO reuse temporaries
      VelocityFunctionType tmp;
      VelocityFunctionType tmp2;
      PressureFunctionType tmp3;

      std::shared_ptr< VelocityFunctionType >   invDiagA;
      const std::shared_ptr< PrimitiveStorage > storage_;
      const uint_t                              level_;
      bool                                      invDiagAAssembled;
   };

 private:
   VelocityBlockOpType A;
   BFBTSubOperator     Z;

   std::shared_ptr< PETScCGSolver< BFBTSubOperator > > ZSolver;

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
