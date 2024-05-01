# in case Git is not available, we default to "unknown"
set(GIT_HASH "unknown" CACHE STRING "Description" FORCE)

# find Git and if available set GIT_HASH variable
find_package(Git QUIET)
if(GIT_FOUND)
  execute_process(
    WORKING_DIRECTORY ${SOURCE_DIR}
    #COMMAND ${GIT_EXECUTABLE} log -1 --pretty=format:%h
    COMMAND ${GIT_EXECUTABLE} describe --abbrev=10 --long --always --dirty --tags
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    )
  if (WIN32)
    execute_process(
      WORKING_DIRECTORY ${SOURCE_DIR}
      #COMMAND ${GIT_EXECUTABLE} log -1 --pretty=format:%h
      COMMAND cmd /C ${GIT_EXECUTABLE} diff ':!*version.*' ${SOURCE_DIR}/src/solver/src
      OUTPUT_VARIABLE GIT_DIFF
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
      )
  else ()
    execute_process(
      WORKING_DIRECTORY ${SOURCE_DIR}
      #COMMAND ${GIT_EXECUTABLE} log -1 --pretty=format:%h
      COMMAND bash -c ${GIT_EXECUTABLE} diff -- ':!*version.*' ${SOURCE_DIR}/src/solver/src
      OUTPUT_VARIABLE GIT_DIFF
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
      )
  endif ()
endif()

if(GIT_HASH STREQUAL "")
    # Set a default value if the output is empty
    set(GIT_HASH "nogit")
endif()

if(GIT_DIFF STREQUAL "")
    # Set a default value if the output is empty
    set(GIT_DIFF "It was not possible to get differences between the source code and the git version.")
endif()

# git diff -- src/solver/src/* CMakeLists.txt
# git diff -- *.{cpp,h}

message(STATUS "Git hash is ${GIT_HASH}")

foreach(key
  IN ITEMS
    OS_NAME
    OS_RELEASE
    OS_VERSION
    OS_PLATFORM
  )
  cmake_host_system_information(RESULT _${key} QUERY ${key})
endforeach()

string(TIMESTAMP TIME_STRING "%Y-%m-%d %H:%M")
string(TIMESTAMP DATE_STRING "%Y-%m-%d")

# generate file version.h based on version.h.in
configure_file(
  ${CMAKE_CURRENT_LIST_DIR}/src/solver/src/version.h.in
  ${TARGET_DIR}/generated/version.h
  @ONLY
  )

set(my_output ${TARGET_DIR}/generated/hash.txt)
file(WRITE ${my_output} ${GIT_HASH})

set(my_output ${TARGET_DIR}/generated/date.txt)
file(WRITE ${my_output} ${DATE_STRING})

set(my_output ${TARGET_DIR}/generated/target_name.txt)
file(WRITE ${my_output} ${TARGET_NAME}_${DATE_STRING}_${GIT_HASH})
