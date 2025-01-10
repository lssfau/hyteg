waLBerla_link_files_to_builddir( data )

waLBerla_add_executable(NAME adaptiveRefinement
      FILES adaptiveRefinement.cpp
      DEPENDS hyteg opgen-mass constant_stencil_operator)

waLBerla_add_executable(NAME fullAMR
      FILES fullAMR.cpp
      DEPENDS hyteg opgen-mass constant_stencil_operator)