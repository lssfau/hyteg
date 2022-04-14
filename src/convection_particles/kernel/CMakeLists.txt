target_sources( convection_particles
    PRIVATE
    DetectAndStoreContacts.h
    DoubleCast.h
    ForceLJ.h
    ExplicitEuler.h
    AssocToBlock.h
    NonLinearSpringDashpot.h
    ExplicitEulerWithShape.h
    TemperatureIntegration.h
    InitParticlesForHCSITS.h
    InsertParticleIntoSparseLinkedCells.h
    VelocityVerlet.h
    InsertParticleIntoLinkedCells.h
    LinearSpringDashpot.h
    SpringDashpotSpring.h
    SingleCast.h
    SpringDashpot.h
    HeatConduction.h
    InitContactsForHCSITS.h
    ParticleSelector.h
    IntegrateParticlesHCSITS.h
    HCSITSRelaxationStep.h     
)

