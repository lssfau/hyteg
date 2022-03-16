target_sources( convection_particles
    PRIVATE
    SyncNextNeighborsBlockForest.h
    ShapePackUnpack.h
    RemoveAndNotify.h
    SyncNextNeighbors.cpp
    SyncNextNeighbors.h
    ClearGhostOwnerSync.h
    ClearNextNeighborSync.h
    ContactFilter.h
    SyncNextNeighborsNoGhosts.h
    SyncNextNeighborsBlockForest.cpp
    SyncGhostOwners.cpp
    RemoveAndNotify.cpp
    SyncNextNeighborsNoGhosts.cpp
    SyncGhostOwners.h
    ReduceProperty.h
    BroadcastProperty.h
    ReduceContactHistory.h     
)
add_subdirectory( notifications )
