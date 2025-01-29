set( CURRENT_LIST_DIR ${CMAKE_CURRENT_LIST_DIR} )

if ( NOT DEFINED hyteg_git_input )
   set( hyteg_git_input ${hyteg_SOURCE_DIR}/src/hyteg/Git.in.hpp )
endif ()

if ( NOT DEFINED hyteg_git_output )
   set( hyteg_git_output ${hyteg_BINARY_DIR}/src/hyteg/Git.hpp )
endif ()

function( CheckGitWrite GIT_CHECKSUM )
   file( WRITE ${CMAKE_BINARY_DIR}/git-state.txt ${GIT_CHECKSUM} )
endfunction()

function( CheckGitRead GIT_CHECKSUM_CACHE )
   if ( EXISTS ${CMAKE_BINARY_DIR}/git-state.txt )
      file( READ ${CMAKE_BINARY_DIR}/git-state.txt file_content )
      set( GIT_CHECKSUM_CACHE ${file_content} PARENT_SCOPE )
   endif ()
endfunction()

function( CheckGitVersion )

   execute_process(
           COMMAND ${GIT_EXECUTABLE} log -1 --format=%h
           WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
           OUTPUT_VARIABLE GIT_COMMIT_HASH
           OUTPUT_STRIP_TRAILING_WHITESPACE
   )

   execute_process(
           COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
           OUTPUT_VARIABLE GIT_BRANCH
           OUTPUT_STRIP_TRAILING_WHITESPACE
   )

   execute_process(
           COMMAND ${GIT_EXECUTABLE} diff --stat
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
           OUTPUT_VARIABLE GIT_DIFF
           OUTPUT_STRIP_TRAILING_WHITESPACE
   )

   string(JOIN "\n" GIT_CHECKSUM ${GIT_COMMIT_HASH} ${GIT_BRANCH} ${GIT_DIFF})

   CheckGitRead( GIT_CHECKSUM_CACHE )

   if ( NOT DEFINED GIT_CHECKSUM_CACHE )
      set( GIT_CHECKSUM_CACHE "INVALID" )
   endif ()

   # Only update the git_version.cpp if the hash has changed. This will
   # prevent us from rebuilding the project more than we need to.
   if ( NOT ${GIT_CHECKSUM} STREQUAL ${GIT_CHECKSUM_CACHE} OR NOT EXISTS ${hyteg_git_output} )
      CheckGitWrite( ${GIT_CHECKSUM} )
      configure_file( ${hyteg_git_input} ${hyteg_git_output} @ONLY )
   endif ()

endfunction()

function( CheckGitSetup )

   add_custom_target( AlwaysCheckGit COMMAND ${CMAKE_COMMAND}
           -DRUN_CHECK_GIT_VERSION=1
           -Dhyteg_git_input=${hyteg_git_input}
           -Dhyteg_git_output=${hyteg_git_output}
           -DGIT_EXECUTABLE=${GIT_EXECUTABLE}
           -P ${CURRENT_LIST_DIR}/CheckGit.cmake
           BYPRODUCTS ${hyteg_git_output}
   )

   add_library( git_version INTERFACE ${hyteg_git_output} ${CURRENT_LIST_DIR}/CheckGit.cmake )
   add_dependencies( git_version AlwaysCheckGit )

   CheckGitVersion()
endfunction()

# This is used to run this function from an external cmake process.
if ( RUN_CHECK_GIT_VERSION )
   CheckGitVersion()
endif ()
