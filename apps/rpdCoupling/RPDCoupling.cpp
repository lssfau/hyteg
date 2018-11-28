
#if 0
#include "core/Environment.h"

#include "tinyhhg_core/mesh/MeshInfo.hpp"

#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"

#include "tinyhhg_core/rpd_coupling/SetupPrimitiveStorageInterface.hpp"

namespace hhg {

void rpdCoupling()
{
  auto meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_24el.msh" );
  auto setupPrimitiveStorage = std::make_shared< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  rpd::SetupPrimitiveStorageInterface rpdSetupStorage( setupPrimitiveStorage );

}

}

int main( int argc, char* argv[] )
{
  walberla::Environment walberlaEnv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();
  hhg::rpdCoupling();
}
#endif

#include <rpd/vtk/ParticleVtkOutput.h>

#include <rpd/data/LinkedCells.h>
#include <rpd/data/ParticleStorage.h>
#include <rpd/data/ShapeStorage.h>

#include <rpd/domain/BlockForestDomain.h>

#include <rpd/kernel/ExplicitEulerWithShape.h>
#include <rpd/kernel/GenerateLinkedCells.h>
#include <rpd/kernel/SyncProperty.h>
#include <rpd/kernel/SpringDashpot.h>
#include <rpd/kernel/SyncNextNeighbors.h>

#include <blockforest/BlockForest.h>
#include <core/Abort.h>
#include <core/Environment.h>
#include <core/math/Random.h>
#include <core/mpi/Reduce.h>
#include <core/grid_generator/SCIterator.h>
#include <core/logging/Logging.h>
#include <core/timing/Timer.h>
#include <core/waLBerlaBuildInfo.h>
#include <pe/utility/CreateWorld.h>
#include <vtk/VTKOutput.h>

#include "tinyhhg_core/mesh/MeshInfo.hpp"

#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"

#include "tinyhhg_core/rpd_coupling/SetupPrimitiveStorageInterface.hpp"

#include <functional>
#include <memory>

namespace walberla {

using namespace walberla::rpd;

int main( int argc, char ** argv )
{
   using namespace walberla::timing;

   Environment env(argc, argv);
   mpi::MPIManager::instance()->useWorldComm();

   logging::Logging::instance()->setStreamLogLevel(logging::Logging::INFO);
   logging::Logging::instance()->setFileLogLevel(logging::Logging::INFO);

   WALBERLA_LOG_INFO_ON_ROOT( "config file: " << argv[1] );
   WALBERLA_LOG_INFO_ON_ROOT( "waLBerla Revision: " << WALBERLA_GIT_SHA1 );

   math::seedRandomGenerator( static_cast<unsigned int>(1337 * mpi::MPIManager::instance()->worldRank()) );

   WALBERLA_LOG_INFO_ON_ROOT("*** READING CONFIG FILE ***");
   auto cfg = env.config();
   if (cfg == nullptr) WALBERLA_ABORT("No config specified!");
   const Config::BlockHandle mainConf  = cfg->getBlock( "PeriodicGranularGas" );

   const real_t spacing = mainConf.getParameter<real_t>("spacing", real_t(1.0) );
   WALBERLA_LOG_INFO_ON_ROOT("spacing: " << spacing);

   const real_t radius = mainConf.getParameter<real_t>("radius", real_t(0.5) );
   WALBERLA_LOG_INFO_ON_ROOT("radius: " << radius);

   int simulationSteps = mainConf.getParameter<int>("simulationSteps", 10 );
   WALBERLA_LOG_INFO_ON_ROOT("simulationSteps: " << simulationSteps);

   real_t dt = mainConf.getParameter<real_t>("dt", real_c(0.01) );
   WALBERLA_LOG_INFO_ON_ROOT("dt: " << dt);

   const int visSpacing = mainConf.getParameter<int>("visSpacing",  1000 );
   WALBERLA_LOG_INFO_ON_ROOT("visSpacing: " << visSpacing);
   const std::string path = mainConf.getParameter<std::string>("path",  "vtk_out" );
   WALBERLA_LOG_INFO_ON_ROOT("path: " << path);

   WALBERLA_LOG_INFO_ON_ROOT("*** BLOCKFOREST ***");
#if 0
   // create forest
   shared_ptr< BlockForest > forest = pe::createBlockForestFromConfig( mainConf );
   if (!forest)
   {
      WALBERLA_LOG_INFO_ON_ROOT( "No BlockForest created ... exiting!");
      return EXIT_SUCCESS;
   }
   BlockForestDomain domain(forest);

   auto simulationDomain = forest->getDomain();
   WALBERLA_CHECK_EQUAL(forest->size(), 1);
   auto localDomain      = forest->begin()->getAABB();
#endif

   auto meshInfo = hhg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_6el.msh" );
   math::AABB simulationDomain( Vec3( 0, 0, 0 ), Vec3( 1, 1, 1 ) );
   auto setupPrimitiveStorage = std::make_shared< hhg::SetupPrimitiveStorage >( meshInfo, uint_c( mpi::MPIManager::instance()->numProcesses() ) );
   hhg::rpd::SetupPrimitiveStorageInterface rpdSetupStorage( setupPrimitiveStorage );
   auto storage = std::make_shared< hhg::PrimitiveStorage >( *setupPrimitiveStorage );
   hhg::writeDomainPartitioningVTK( storage, "../../output", "Domain" );

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - START ***");

   //init data structures
   data::ParticleStorage ps(5000);
   data::ShapeStorage    ss;

   auto  smallSphere = ss.create<data::Sphere>( radius );
   for (auto pt : grid_generator::SCGrid(simulationDomain, Vector3<real_t>(spacing, spacing, spacing) * real_c(0.5), spacing))
   {
      const auto rank = mpi::MPIManager::instance()->rank();
      if ( rpdSetupStorage.isContainedInProcessSubdomain( rank, pt ) )
      {
        data::Particle&& p = ps.create();
        p.getPosition()          = pt;
        p.getLinearVelocity()    = Vec3(math::realRandom(-1.,1.), math::realRandom(-1.,1.), math::realRandom(-1.,1.));
        p.getInteractionRadius() = radius;
        p.getShapeID()           = smallSphere;
        p.getOwner()             = rank;
      }
   }
   uint_t numParticles = ps.size();
   mpi::reduceInplace(numParticles, mpi::SUM);
   WALBERLA_LOG_INFO_ON_ROOT("#particles created: " << numParticles);

   data::Particle&& p0 = ps.create(UniqueID<data::Particle>::createGlobal());
   p0.getPosition() = simulationDomain.minCorner();
   p0.getShapeID()  = ss.create<data::HalfSpace>( Vec3(0,0,1) );
   p0.getOwner()    = mpi::MPIManager::instance()->rank();
   data::particle_flags::set(p0.getFlags(), data::particle_flags::INFINITE);
   data::particle_flags::set(p0.getFlags(), data::particle_flags::FIXED);
   data::particle_flags::set(p0.getFlags(), data::particle_flags::NON_COMMUNICATING);

   data::Particle&& p1 = ps.create(UniqueID<data::Particle>::createGlobal());
   p1.getPosition() = simulationDomain.maxCorner();
   p1.getShapeID()  = ss.create<data::HalfSpace>( Vec3(0,0,-1) );
   p1.getOwner()    = mpi::MPIManager::instance()->rank();
   data::particle_flags::set(p1.getFlags(), data::particle_flags::INFINITE);
   data::particle_flags::set(p1.getFlags(), data::particle_flags::FIXED);
   data::particle_flags::set(p1.getFlags(), data::particle_flags::NON_COMMUNICATING);

   data::Particle&& p2 = ps.create(UniqueID<data::Particle>::createGlobal());
   p2.getPosition() = simulationDomain.minCorner();
   p2.getShapeID()  = ss.create<data::HalfSpace>( Vec3(1,0,0) );
   p2.getOwner()    = mpi::MPIManager::instance()->rank();
   data::particle_flags::set(p2.getFlags(), data::particle_flags::INFINITE);
   data::particle_flags::set(p2.getFlags(), data::particle_flags::FIXED);
   data::particle_flags::set(p2.getFlags(), data::particle_flags::NON_COMMUNICATING);

   data::Particle&& p3 = ps.create(UniqueID<data::Particle>::createGlobal());
   p3.getPosition() = simulationDomain.maxCorner();
   p3.getShapeID()  = ss.create<data::HalfSpace>( Vec3(-1,0,0) );
   p3.getOwner()    = mpi::MPIManager::instance()->rank();
   data::particle_flags::set(p3.getFlags(), data::particle_flags::INFINITE);
   data::particle_flags::set(p3.getFlags(), data::particle_flags::FIXED);
   data::particle_flags::set(p3.getFlags(), data::particle_flags::NON_COMMUNICATING);

   data::Particle&& p4 = ps.create(UniqueID<data::Particle>::createGlobal());
   p4.getPosition() = simulationDomain.minCorner();
   p4.getShapeID()  = ss.create<data::HalfSpace>( Vec3(0,1,0) );
   p4.getOwner()    = mpi::MPIManager::instance()->rank();
   data::particle_flags::set(p4.getFlags(), data::particle_flags::INFINITE);
   data::particle_flags::set(p4.getFlags(), data::particle_flags::FIXED);
   data::particle_flags::set(p4.getFlags(), data::particle_flags::NON_COMMUNICATING);

   data::Particle&& p5 = ps.create(UniqueID<data::Particle>::createGlobal());
   p5.getPosition() = simulationDomain.maxCorner();
   p5.getShapeID()  = ss.create<data::HalfSpace>( Vec3(0,-1,0) );
   p5.getOwner()    = mpi::MPIManager::instance()->rank();
   data::particle_flags::set(p5.getFlags(), data::particle_flags::INFINITE);
   data::particle_flags::set(p5.getFlags(), data::particle_flags::FIXED);
   data::particle_flags::set(p5.getFlags(), data::particle_flags::NON_COMMUNICATING);

   WALBERLA_LOG_INFO_ON_ROOT("*** SETUP - END ***");

   WALBERLA_LOG_INFO_ON_ROOT("*** VTK ***");
   auto vtkOutput       = make_shared<rpd::vtk::ParticleVtkOutput>(ps) ;
   auto vtkWriter       = walberla::vtk::createVTKOutput_PointData(vtkOutput, "Bodies", 1, "vtk", "simulation_step", false, false);

   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - START ***");
   // Init kernels
   kernel::GenerateLinkedCells    generateLinkedCells;
   kernel::ExplicitEulerWithShape explicitEulerWithShape( dt );
   kernel::SpringDashpot          dem;
   kernel::SyncNextNeighbors      SNN;
   kernel::ReduceProperty         RP;
   WcTimer timer;
   WcTimingPool tp;
   SNN(ps, rpdSetupStorage);
   for (int i=0; i < simulationSteps; ++i)
   {
      if( i % 1 == 0 )
      {
         WALBERLA_LOG_DEVEL_ON_ROOT( "Timestep " << i << " / " << simulationSteps );
      }

      if (i % visSpacing == 0)
      {
         vtkWriter->write();
      }
#if 0
      tp["GenerateLinkedCells"].start();
      generateLinkedCells(ps, lc);
      tp["GenerateLinkedCells"].end();
#endif
      tp["DEM"].start();
      dem(ps, ss);
      tp["DEM"].end();

      tp["ReduceForce"].start();
      RP(ps, [](data::Particle p) -> auto& {return p.getForce();});
      tp["ReduceForce"].end();

      tp["Euler"].start();
      explicitEulerWithShape(ps, ss);
      tp["Euler"].end();

      tp["SNN"].start();
      SNN(ps, rpdSetupStorage);
      tp["SNN"].end();
   }
   timer.end();
   auto tp_reduced = tp.getReduced();
   WALBERLA_LOG_INFO_ON_ROOT(*tp_reduced);
   WALBERLA_LOG_INFO_ON_ROOT("runtime: " << timer.average());
   WALBERLA_LOG_INFO_ON_ROOT("PUpS: " << real_c(numParticles) * real_c(simulationSteps) / timer.average());
   WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - END ***");

   WALBERLA_LOG_INFO_ON_ROOT("*** CHECKING RESULT - START ***");
   auto pIt = ps.begin();
   for (auto it = grid_generator::SCIterator(simulationDomain, Vector3<real_t>(spacing, spacing, spacing) * real_c(0.5), spacing);
        it != grid_generator::SCIterator();
        ++it, ++pIt)
   {
     const auto rank = mpi::MPIManager::instance()->rank();
     if ( rpdSetupStorage.isContainedInProcessSubdomain( rank, pIt->getPosition() ) )
     {
       WALBERLA_CHECK_UNEQUAL(pIt, ps.end());
       WALBERLA_CHECK_FLOAT_EQUAL((*pIt).getPosition(), *it);
     }
   }
   WALBERLA_LOG_INFO_ON_ROOT("*** CHECKING RESULT - END ***");

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
   return walberla::main( argc, argv );
}
