target_sources( hyteg
    PRIVATE
    DataTools.hpp
    DataLoader.hpp
)

add_subdirectory(DensityProfiles)
add_subdirectory(GravityProfiles)
add_subdirectory(PressureProfiles)
add_subdirectory(TemperatureProfiles)
add_subdirectory(ViscosityProfiles)