if( HYTEG_BUILD_WITH_ADIOS2 )
  waLBerla_compile_test(FILES BasicADIOSTest.cpp DEPENDS hyteg core)
  waLBerla_execute_test(NAME BasicADIOSTest)

  if( NOT EXISTS ./AdiosWriterTest-Output )
    file( MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/AdiosWriterTest-Output" )
  endif()

  waLBerla_compile_test(FILES AdiosWriterTest.cpp DEPENDS hyteg core)
  waLBerla_execute_test(NAME AdiosWriterTest-Sequential COMMAND $<TARGET_FILE:AdiosWriterTest> PROCESSES 1)

  if( HYTEG_BUILD_WITH_PYTHON3 )
    add_subdirectory( AdiosWriterTest-References )

    if( WALBERLA_DOUBLE_ACCURACY )
      set ( WALBERLA_FP_BITS 64 )
    else ()
      set ( WALBERLA_FP_BITS 32 )
    endif ()

    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/AdiosWriterTest-verify.py
      ${CMAKE_CURRENT_BINARY_DIR}/AdiosWriterTest-verify.py )

    waLBerla_execute_test( NAME AdiosWriterVerifyTest-Sequential
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/AdiosWriterTest-verify.py --num-writers 1 --fp-bits ${WALBERLA_FP_BITS} )
    set_tests_properties( AdiosWriterTest-Sequential PROPERTIES FIXTURES_SETUP AdiosOutput-Sequential )
    set_tests_properties( AdiosWriterVerifyTest-Sequential PROPERTIES FIXTURES_REQUIRED AdiosOutput-Sequential )

    if( WALBERLA_BUILD_WITH_MPI )
      waLBerla_execute_test(NAME AdiosWriterTest-Parallel COMMAND $<TARGET_FILE:AdiosWriterTest> PROCESSES 3)
      waLBerla_execute_test( NAME AdiosWriterVerifyTest-Parallel
	COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/AdiosWriterTest-verify.py --num-writers 3 --fp-bits ${WALBERLA_FP_BITS} )
      set_tests_properties( AdiosWriterTest-Parallel PROPERTIES FIXTURES_SETUP AdiosOutput-Parallel )
      set_tests_properties( AdiosWriterVerifyTest-Parallel PROPERTIES FIXTURES_REQUIRED AdiosOutput-Parallel )
    endif()

  endif()

  waLBerla_compile_test(FILES CheckpointingTest.cpp DEPENDS hyteg core)
  waLBerla_execute_test(NAME CheckpointingTest)

  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${CMAKE_CURRENT_SOURCE_DIR}/CheckpointRestoreTest.prm
    ${CMAKE_CURRENT_BINARY_DIR}/CheckpointRestoreTest.prm )
  waLBerla_compile_test(FILES CheckpointRestoreTest.cpp DEPENDS hyteg core)

  waLBerla_execute_test(NAME CheckpointRestoreTest-Sequential COMMAND $<TARGET_FILE:CheckpointRestoreTest> PROCESSES 1)

  if( WALBERLA_BUILD_WITH_MPI )
    waLBerla_execute_test(NAME CheckpointRestoreTest-Parallel-Same-Communicator COMMAND $<TARGET_FILE:CheckpointRestoreTest> PROCESSES 3)

    set_tests_properties( CheckpointRestoreTest-Sequential                 PROPERTIES FIXTURES_SETUP CheckPointRestore )
    set_tests_properties( CheckpointRestoreTest-Parallel-Same-Communicator PROPERTIES FIXTURES_SETUP CheckPointRestore )

    waLBerla_execute_test(NAME CheckpointRestoreTest-Sequential-Import-Only COMMAND $<TARGET_FILE:CheckpointRestoreTest> CheckpointRestoreTest.prm -Parameters.onlyImport=true -Parameters.restoreFromFileWithName=CheckpointRestoreTest-np3.bp PROCESSES 1)
    waLBerla_execute_test(NAME CheckpointRestoreTest-Parallel-Import-Only1 COMMAND $<TARGET_FILE:CheckpointRestoreTest> CheckpointRestoreTest.prm -Parameters.onlyImport=true -Parameters.restoreFromFileWithName=CheckpointRestoreTest-np3.bp PROCESSES 3)
    waLBerla_execute_test(NAME CheckpointRestoreTest-Parallel-Import-Only2 COMMAND $<TARGET_FILE:CheckpointRestoreTest> CheckpointRestoreTest.prm -Parameters.onlyImport=true -Parameters.restoreFromFileWithName=CheckpointRestoreTest-np1.bp PROCESSES 3)

    set_tests_properties( CheckpointRestoreTest-Sequential-Import-Only PROPERTIES FIXTURES_REQUIRED CheckPointRestore )
    set_tests_properties( CheckpointRestoreTest-Parallel-Import-Only1  PROPERTIES FIXTURES_REQUIRED CheckPointRestore )
    set_tests_properties( CheckpointRestoreTest-Parallel-Import-Only2  PROPERTIES FIXTURES_REQUIRED CheckPointRestore )

  endif ()

endif ()

waLBerla_compile_test(FILES TableTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME TableTest)
set_tests_properties(TableTest PROPERTIES FIXTURES_SETUP TableTest)

if(WALBERLA_DOUBLE_ACCURACY)
   set(TABLE_REFERENCE_FILE "TableTest.dat.ref.dp")
else()
   set(TABLE_REFERENCE_FILE "TableTest.dat.ref.sp")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/${TABLE_REFERENCE_FILE}
      ${CMAKE_CURRENT_BINARY_DIR}/${TABLE_REFERENCE_FILE})
ADD_TEST(TableTestCheckOutput ${CMAKE_COMMAND} -E compare_files TableTest.dat ${TABLE_REFERENCE_FILE})
set_tests_properties(TableTestCheckOutput PROPERTIES FIXTURES_REQUIRED TableTest)

waLBerla_compile_test(FILES KeyValueStoreTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME KeyValueStoreTest)
set_tests_properties(KeyValueStoreTest PROPERTIES FIXTURES_SETUP KeyValueStoreTest)

execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/KeyValueStoreTest.tex.ref
      ${CMAKE_CURRENT_BINARY_DIR}/KeyValueStoreTest.tex.ref)
ADD_TEST(KeyValueStoreTestCheckOutput ${CMAKE_COMMAND} -E compare_files KeyValueStoreTest.tex KeyValueStoreTest.tex.ref)
set_tests_properties(KeyValueStoreTestCheckOutput PROPERTIES FIXTURES_REQUIRED KeyValueStoreTest)

waLBerla_compile_test(FILES VTKOutputTest.cpp DEPENDS hyteg core mixed_operator)
waLBerla_execute_test(NAME VTKOutputTest)

waLBerla_compile_test(FILES VTKDuplicatePointCountingTest.cpp DEPENDS hyteg core)
waLBerla_execute_test(NAME VTKDuplicatePointCountingTest)

waLBerla_compile_test(FILES VTKOutputPictureNormGenerator.cpp DEPENDS hyteg core)

execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/VTKOutputPictureNormTest.py
        ${CMAKE_CURRENT_BINARY_DIR}/VTKOutputPictureNormTest.py )

execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/VTKOutputPictureNormTest-References
        ${CMAKE_CURRENT_BINARY_DIR}/VTKOutputPictureNormTest-References )

if ( NOT EXISTS ./VTKOutputPictureNormTest-Fresh )
  file( MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/VTKOutputPictureNormTest-Fresh" )
endif ()

if ( HYTEG_BUILD_WITH_PYTHON3 )
  waLBerla_execute_test( NAME VTKOutputPictureNormTest
          COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/VTKOutputPictureNormTest.py --with-xvfb )
endif ()