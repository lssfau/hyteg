target_sources( hyteg
    PRIVATE
    JacobiPreconditioner.hpp
    IdentityPreconditioner.hpp     
)
add_subdirectory( stokes )
