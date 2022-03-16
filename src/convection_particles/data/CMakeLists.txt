target_sources( convection_particles
    PRIVATE
    ShapeStorage.h
    ParticleStorage.h
    ParticleAccessor.h
    STLOverloads.h
    Flags.h
    ContactStorage.h
    SparseLinkedCells.h
    ContactHistory.h
    LinkedCells.h
    DataTypes.h
    ContactAccessor.h
    IAccessor.h
    IContactAccessor.h
    ParticleAccessorWithShape.h     
)
add_subdirectory( shape )
