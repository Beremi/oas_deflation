# in case Git is not available, we default to "unknown"
set(GIT_HASH "unknown" CACHE STRING "Description" FORCE)

# find Git and if available set GIT_HASH variable
find_package(Git QUIET)
if(GIT_FOUND)
  execute_process(
    #COMMAND ${GIT_EXECUTABLE} log -1 --pretty=format:%h
    COMMAND ${GIT_EXECUTABLE} describe --abbrev=10 --long --always --dirty --tags
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    )
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
file(WRITE ${my_output} ${TARGET_NAME}_${GIT_HASH}_${DATE_STRING})
