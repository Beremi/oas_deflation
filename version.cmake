# in case Git is not available, we default to "unknown"
set(GIT_HASH "unknown")

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

# generate file version.h based on version.h.in
configure_file(
  ${CMAKE_CURRENT_LIST_DIR}/src/solver/src/version.h.in
  ${TARGET_DIR}/generated/version.h
  @ONLY
  )
