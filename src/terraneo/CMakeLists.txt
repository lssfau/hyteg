################################################################################
#
# Module for mantle convection simulations based on HyTeG
#
################################################################################

add_library( terraneo )
target_link_libraries( terraneo PUBLIC core hyteg )
# target_sources( terraneo )

add_subdirectory( sphericalharmonics )
