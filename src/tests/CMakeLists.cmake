
#add_custom_command( TARGET DiscreteModel
#    POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E echo "Copying tests data.."
#    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/src/tests ${CMAKE_BINARY_DIR}/tests
#    )

#file(COPY ${CMAKE_SOURCE_DIR}/src/tests DESTINATION ${CMAKE_BINARY_DIR})

# add_test(
#   NAME setup
#   COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/test/setup.py
#   )
# set_tests_properties(
#   setup
#   PROPERTIES
#     FIXTURES_SETUP my-fixture
#   )

add_test(
   NAME setup
   #COMMAND ${CMAKE_COMMAND} -E echo "Copying tests data.."
   COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/src/tests ${CMAKE_BINARY_DIR}/tests
   )
set_tests_properties(
   setup
   PROPERTIES
     FIXTURES_SETUP my-fixture
   )


add_test(
  NAME SpringMechElastic
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/SpringMechElastic
  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_BINARY_DIR}/tests/check.py $<TARGET_FILE:DiscreteModel>
  )

add_test(
  NAME SpringMechElastic_3D
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/SpringMechElastic_3D
  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_BINARY_DIR}/tests/check.py $<TARGET_FILE:DiscreteModel>
  )

set_tests_properties(
  SpringMechElastic
  SpringMechElastic_3D
  PROPERTIES
    LABELS "tests"
    TIMEOUT 10
    FIXTURES_REQUIRED my-fixture
  )

# add_test(
#   NAME cleanup
#   COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/test/cleanup.py
#   )
# set_tests_properties(
#   cleanup
#   PROPERTIES
#     FIXTURES_CLEANUP my-fixture
#   )
