
#include <rpd/vtk/ParticleVtkOutput.h>
#include <rpd/data/ParticleStorage.h>
#include <rpd/domain/BlockForestDomain.h>
#include <rpd/kernel/SyncNextNeighbors.h>
#include <rpd/kernel/ExplicitEuler.h>
#include <rpd/kernel/GenerateLinkedCells.h>
#include <rpd/kernel/SyncProperty.h>
#include <rpd/kernel/SpringDashpot.h>
#include <vtk/VTKOutput.h>

#include <blockforest/BlockForest.h>
#include <core/Environment.h>
#include <core/logging/Logging.h>
#include <core/mpi/Reduce.h>
#include <pe/utility/CreateWorld.h>

#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"

#include "tinyhhg_core/rpd_coupling/SetupPrimitiveStorageInterface.hpp"

#include <iostream>
#include <memory>

using namespace walberla;
using namespace walberla::rpd;


walberla::id_t createSphere( const Vec3 & position, const Vec3 & velocity, const real_t & radius,
                             data::ParticleStorage & ps, IDomain & domain )
{
   walberla::id_t uid = 0;
   auto owned = domain.isContainedInProcessSubdomain( uint_c(mpi::MPIManager::instance()->rank()), position );
   if (owned)
   {
      data::Particle&& p    = ps.create();
      p.getPosition()          = position;
      p.getInteractionRadius() = radius;
      p.getRotation()          = Rot3(Quat());
      p.getLinearVelocity()    = velocity;
      p.getAngularVelocity()   = Vec3(4,5,6);
      p.getOwner()             = mpi::MPIManager::instance()->rank();
      uid = p.getUid();
   }

   WALBERLA_CHECK_NOT_IDENTICAL( domain.findContainingProcessRank( position ), -1 );

   mpi::allReduceInplace(uid, mpi::SUM);
   return uid;
}

int main( int argc, char ** argv )
{
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   mpi::MPIManager::instance()->useWorldComm();

   logging::Logging::instance()->setStreamLogLevel(logging::Logging::DETAIL);
   logging::Logging::instance()->includeLoggingToFile("RPD_Kernel_SyncNextNeighbor");
   logging::Logging::instance()->setFileLogLevel(logging::Logging::DETAIL);

   const real_t dt              = 0.01;
   const uint_t simulationSteps = 8000;
   const uint_t visSpacing      = 10;

   auto meshInfo = hhg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_24el.msh" );
   math::AABB simulationDomain( Vec3( 0, 0, 0 ), Vec3( 1, 1, 1 ) );
   auto setupPrimitiveStorage = std::make_shared< hhg::SetupPrimitiveStorage >( meshInfo, uint_c( mpi::MPIManager::instance()->numProcesses() ) );
   hhg::rpd::SetupPrimitiveStorageInterface rpdSetupStorage( setupPrimitiveStorage );
   auto storage = std::make_shared< hhg::PrimitiveStorage >( *setupPrimitiveStorage );
   hhg::writeDomainPartitioningVTK( storage, "../../output", "SyncNextNeighborDomain" );

   //init data structures
   data::ParticleStorage ps(100);

   for ( uint_t y = 1; y < 10; y++ )
   {
     for ( uint_t z = 1; z < 10; z++ )
     {
       //initialize particle
       const Vec3 position( 0.1, y * 0.1, z * 0.1 );
       const Vec3 velocity( 0.01, 0.0, 0.0 );
       const real_t radius = 0.05;
       auto uid = createSphere(position, velocity, radius, ps, rpdSetupStorage);
       WALBERLA_LOG_DEVEL_ON_ROOT("uid: " << uid);
     }
   }

   //init kernels
   kernel::SyncNextNeighbors SNN;

    WALBERLA_LOG_INFO_ON_ROOT("*** VTK ***");
    auto vtkOutput       = make_shared<rpd::vtk::ParticleVtkOutput>(ps) ;
    vtkOutput->addOutput< rpd::data::SelectParticleOwner >( "particleRank" );
    auto vtkWriter       = walberla::vtk::createVTKOutput_PointData(vtkOutput, "Bodies", 1, "../../output", "simulation_step", false, false);

    WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - START ***");
    // Init kernels
    kernel::ExplicitEuler          explicitEuler( dt );
    kernel::ReduceProperty         RP;

    WcTimer timer;
    WcTimingPool tp;
    SNN(ps, rpdSetupStorage);
    for (uint_t i=0; i < simulationSteps; ++i)
    {
      if( i % 1 == 0 )
      {
        WALBERLA_LOG_DEVEL_ON_ROOT( "Timestep " << i << " / " << simulationSteps );
      }

      if (i % visSpacing == 0)
      {
        vtkWriter->write();
      }

      for ( auto p : ps )
      {
        if ( rpdSetupStorage.isContainedInProcessSubdomain( mpi::MPIManager::instance()->rank(), p.getPosition()))
        {
          WALBERLA_CHECK( !data::particle_flags::isSet( p->getFlags(), data::particle_flags::GHOST ));
        }
        else
        {
          WALBERLA_CHECK( data::particle_flags::isSet( p->getFlags(), data::particle_flags::GHOST ));
        }
      }

      tp["ReduceForce"].start();
      RP(ps, [](data::Particle p) -> auto& {return p.getForce();});
      tp["ReduceForce"].end();

      tp["Euler"].start();
      explicitEuler(ps);
      tp["Euler"].end();

      tp["SNN"].start();
      SNN(ps, rpdSetupStorage, 0.001);
      tp["SNN"].end();
    }
    timer.end();
    auto tp_reduced = tp.getReduced();
    WALBERLA_LOG_INFO_ON_ROOT(*tp_reduced);
    WALBERLA_LOG_INFO_ON_ROOT("runtime: " << timer.average());
    WALBERLA_LOG_INFO_ON_ROOT("PUpS: " << real_c(1 /* == num particles */ ) * real_c(simulationSteps) / timer.average());
    WALBERLA_LOG_INFO_ON_ROOT("*** SIMULATION - END ***");

    for ( auto p : ps )
    {
      if ( !data::particle_flags::isSet( p->getFlags(), data::particle_flags::GHOST) )
      {
        WALBERLA_CHECK_FLOAT_EQUAL( p.getPosition()[0], 0.9 );
      }
    }

   return EXIT_SUCCESS;
}
