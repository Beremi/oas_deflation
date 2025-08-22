add_custom_command( TARGET OAS
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "Copying benchmark data.."
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/src/benchmark ${CMAKE_BINARY_DIR}/benchmark
    )

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
  NAME CantileverCoupledLinear
  COMMAND $<TARGET_FILE:OAS> benchmark/CantileverCoupledLinear/master.inp
  )

add_test(
  NAME CantileverMechanicsNonlinear
  COMMAND $<TARGET_FILE:OAS> benchmark/CantileverMechanicsNonlinear/master.inp
  )

add_test(
  NAME 2DUniPressConfined
  COMMAND $<TARGET_FILE:OAS> benchmark/2DUniPressConfined/master.inp
  )

add_test(
  NAME 2DUniPressFree
  COMMAND $<TARGET_FILE:OAS> benchmark/2DUniPressFree/master.inp
  )

add_test(
  NAME 3DUniPressConfined
  COMMAND $<TARGET_FILE:OAS> benchmark/3DUniPressConfined/master.inp
  )

add_test(
  NAME 3DUniPressFree
  COMMAND $<TARGET_FILE:OAS> benchmark/3DUniPressFree/master.inp
  )

add_test(
  NAME DiamondTest
  COMMAND $<TARGET_FILE:OAS> benchmark/DiamondTest/master.inp
  )

add_test(
  NAME RectangleTest
  COMMAND $<TARGET_FILE:OAS> benchmark/RectangleTest/master.inp
  )

add_test(
  NAME SpringMechanicsShear
  COMMAND $<TARGET_FILE:OAS> benchmark/SpringMechanicsShear/master.inp
  )


set_tests_properties(
  2DUniPressConfined
  2DUniPressFree
  3DUniPressConfined
  3DUniPressFree
  CantileverCoupledLinear
  CantileverMechanicsNonlinear
  DiamondTest
  RectangleTest
  SpringMechanicsShear
  PROPERTIES
    LABELS "benchmark"
    TIMEOUT 10
    #FIXTURES_REQUIRED my-fixture
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
