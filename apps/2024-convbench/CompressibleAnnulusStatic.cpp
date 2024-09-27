#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosWriter.hpp"
#include "mixed_operator/VectorMassOperator.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesFullOperator.hpp"
// #include "hyteg_operators/operators/full_stokes/P2VectorElementwiseFullStokesAnnulusMap.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg_operators/operators/terraneo/P2VectorToP1ElementwiseFrozenVelocityAnnulusMap.hpp"

// #include "SimpleCompStokesOperator.hpp"

using walberla::real_t;
using walberla::uint_t;

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      cfg->readParameterFile( "./CompressibleAnnulusStatic.prm" );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   const real_t rMin = mainConf.getParameter<real_t>("rMin");
   const real_t rMax = mainConf.getParameter<real_t>("rMax");

   const uint_t nTan = mainConf.getParameter<uint_t>("nTan");
   const uint_t nRad = mainConf.getParameter<uint_t>("nRad");

   auto meshInfo = MeshInfo::meshAnnulus(rMin, rMax, MeshInfo::CRISS, nTan, nRad);

   auto setupStorage = std::make_shared< hyteg::SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   AnnulusMap::setMap(*setupStorage);

   auto storage = std::make_shared< hyteg::PrimitiveStorage >( *setupStorage, 1 );

   const uint_t minLevel = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel = mainConf.getParameter< uint_t >( "maxLevel" );

   BoundaryCondition bcVelocity;
   bcVelocity.createAllInnerBC();
   bcVelocity.createDirichletBC("DirichletBC", { MeshInfo::flagInnerBoundary, MeshInfo::flagOuterBoundary });

   P2P1TaylorHoodFunction<real_t> u("u", storage, minLevel, maxLevel + 1, bcVelocity);
   P2P1TaylorHoodFunction<real_t> f("f", storage, minLevel, maxLevel, bcVelocity);
   P2P1TaylorHoodFunction<real_t> fTest("fTest", storage, minLevel, maxLevel, bcVelocity);
   P2P1TaylorHoodFunction<real_t> fTestOpGen("fTestOpGen", storage, minLevel, maxLevel, bcVelocity);
   P2P1TaylorHoodFunction<real_t> uAnalytical("uAnalytical", storage, minLevel, maxLevel + 1, bcVelocity);

   P2Function<real_t> rhoP2("rhoP2", storage, minLevel, maxLevel);
   P2VectorFunction<real_t> gradRhoOverRho("gradRhoOverRho", storage, minLevel, maxLevel);

   P2VectorFunction<real_t> uError("uError", storage, maxLevel + 1, maxLevel + 1, bcVelocity);
   P2VectorFunction<real_t> uErrorTmp("uErrorTmp", storage, maxLevel + 1, maxLevel + 1, bcVelocity);

   real_t VrMin = mainConf.getParameter<real_t>("VrMin");

   real_t VpMin = mainConf.getParameter<real_t>("VpMin");
   real_t VpMax = mainConf.getParameter<real_t>("VpMax");

   real_t det = rMin * (rMax * rMax) - rMax * (rMin * rMin);

   std::function<real_t(const Point3D&)> rhoFunc = [](const Point3D& x)
   {  
    real_t r = std::sqrt(x[0] * x[0] + x[1] * x[1]);

   //  return 1.0;
   //  return r * r;
      return 1.0 / r;
   };

   std::function<real_t(const Point3D&)> uX = [&](const Point3D& x)
   {
    real_t r = std::sqrt(x[0] * x[0] + x[1] * x[1]);
    real_t theta = std::atan2(x[1], x[0]);
    
    real_t c1 = rMax * rMax * VpMin - rMin * rMin * VpMax;
    real_t c2 = -rMax * VpMin + rMin * VpMax;

    real_t rho = rhoFunc(x);

    real_t ur = rMin * VrMin / (rho*r);
    real_t uphi = c1 * r + c2 * r * r;

    return ur * std::cos(theta) - uphi * std::sin(theta);
   };

   std::function<real_t(const Point3D&)> uY = [&](const Point3D& x)
   {
    real_t r = std::sqrt(x[0] * x[0] + x[1] * x[1]);
    real_t theta = std::atan2(x[1], x[0]);
    
    real_t c1 = rMax * rMax * VpMin - rMin * rMin * VpMax;
    real_t c2 = -rMax * VpMin + rMin * VpMax;

    real_t rho = rhoFunc(x);

    real_t ur = rMin * VrMin / (rho*r);
    real_t uphi = c1 * r + c2 * r * r;

    return ur * std::sin(theta) + uphi * std::cos(theta);
   };

   uAnalytical.uvw().interpolate({uX, uY}, maxLevel, All);
   uAnalytical.uvw().interpolate({uX, uY}, maxLevel + 1, All);

   u.uvw().interpolate({uX, uY}, maxLevel, DirichletBoundary);

   std::function<real_t(const Point3D&)> viscFunc = [](const Point3D& x)
   {
      real_t r = std::sqrt(x[0] * x[0] + x[1] * x[1]);

      return 1.0 / std::pow(r, 3.0);
   };

   std::function<real_t(const Point3D&)> gradRhoOverRhoFuncX = [](const Point3D& x)
   {
      real_t r = std::sqrt(x[0] * x[0] + x[1] * x[1]);
      return -x[0] / (r * r);
      // return 0.0;
   };

   std::function<real_t(const Point3D&)> gradRhoOverRhoFuncY = [](const Point3D& x)
   {
      real_t r = std::sqrt(x[0] * x[0] + x[1] * x[1]);
      return -x[1] / (r * r);
      // return 0.0;
   };

   rhoP2.interpolate(rhoFunc, maxLevel, All);

   gradRhoOverRho.interpolate({gradRhoOverRhoFuncX, gradRhoOverRhoFuncY}, maxLevel, All);

   P0Function< real_t > viscP0("viscP0", storage, minLevel, maxLevel);
   viscP0.interpolate(viscFunc, maxLevel, All);

   P2Function< real_t > viscP2("viscP2", storage, minLevel, maxLevel);
   viscP2.interpolate(viscFunc, maxLevel, All);

   using StokesOperatorType   = operatorgeneration::P2P1StokesFullAnnulusMapOperator;
//    using StokesOperatorType = operatorgeneration::P2P1StokesFullP0ViscosityAnnulusMapOperator;

   StokesOperatorType stokesOperator(storage, minLevel, maxLevel, viscP2);
//    StokesOperatorType stokesOperator(storage, minLevel, maxLevel, viscP0);
   P2ElementwiseVectorMassOperator vecMassOperator(storage, maxLevel + 1, maxLevel + 1);

   operatorgeneration::P2VectorToP1ElementwiseFrozenVelocityAnnulusMap gradRhoRhoOp(storage, minLevel, maxLevel, rhoP2);

   gradRhoRhoOp.apply(uAnalytical.uvw(), fTestOpGen.p(), maxLevel, All);

   real_t stokesMinresRelTol = mainConf.getParameter<real_t>("stokesMinresRelTol");
   uint_t stokesMinresIter = mainConf.getParameter<uint_t>("stokesMinresIter");

   MinResSolver<StokesOperatorType> minresSolver(storage, minLevel, maxLevel, stokesMinresIter, stokesMinresRelTol);
   minresSolver.setPrintInfo(true);

   minresSolver.solve(stokesOperator, u, fTestOpGen, maxLevel);

   // PETScLUSolver<P2P1THCompStokesOperator> directSolver(storage, maxLevel);

   // directSolver.solve(stokesOperator, u, fTest, maxLevel);

   P2P1StokesToP2P1StokesProlongation prolongationStokes;

   prolongationStokes.prolongate(u, maxLevel, All);

   uError.assign({1.0, -1.0}, {u.uvw(), uAnalytical.uvw()}, maxLevel + 1, All);
   
   vecMassOperator.apply(uError, uErrorTmp, maxLevel + 1, All);

   real_t l2Error = uError.dotGlobal(uErrorTmp, maxLevel + 1, All);

   WALBERLA_LOG_INFO_ON_ROOT(walberla::format("\n\nl2Error = %4.7e\n\n", l2Error));

   std::string outputPath = mainConf.getParameter<std::string>("outputPath");
   std::string outputFilename = mainConf.getParameter<std::string>("outputFilename");

   std::shared_ptr< AdiosWriter > adios2Output = std::make_shared<AdiosWriter>(outputPath, outputFilename, storage);

   adios2Output->add(u);
   adios2Output->add(f);
   adios2Output->add(fTest);
   adios2Output->add(fTestOpGen);
   // adios2Output->add(uError);
   adios2Output->add(uAnalytical);

//    adios2Output->write(maxLevel, 0U);

   return 0;
}