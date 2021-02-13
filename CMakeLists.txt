cmake_minimum_required(VERSION 3.10.0 FATAL_ERROR)

message(STATUS "project source directory -> ${CMAKE_SOURCE_DIR}")
message(STATUS "project build directory -> ${CMAKE_BINARY_DIR}")
# Do this before project() to skip the compiler config process
#if ("${CMAKE_BINARY_DIR}" STREQUAL "${CMAKE_SOURCE_DIR}")
if (("${CMAKE_BINARY_DIR}" MATCHES "${CMAKE_SOURCE_DIR}[\/\\]") OR ("${CMAKE_BINARY_DIR}" STREQUAL "${CMAKE_SOURCE_DIR}"))
  message(FATAL_ERROR "In-source builds are disabled.
    Please create a subfolder outside of repository
    and use `cmake PATH_TO_SRC` inside it.
    NOTE: cmake will now create CMakeCache.txt and CMakeFiles/*.
          You must delete them, or cmake will refuse to work.")
endif() # yes, end-markers and even else() need empty parens

include (CheckSymbolExists)
include (CheckCXXSymbolExists)
include (CheckIncludeFiles)
include (CheckIncludeFileCXX)
include (CheckCXXSourceCompiles)
include (GenerateExportHeader)
include (ProcessorCount)

set(CMAKE_COLOR_MAKEFILE ON)
set(CMAKE_VERBOSE_MAKEFILE OFF)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})


project(PartMod CXX)
set (PartMod_VERSION_MAJOR 0)
set (PartMod_VERSION_MINOR 1)
set (PartMod_VERSION_PATCH 0) # Unused
set (PartMod_HOMEPAGE "https://kelidas.gitlab.io/partmod-site")
#set (PartMod_COPYRIGHT "Copyright (C) 2019")
#set (PartMod_EMAIL "elias.j@fce.vutbr.cz")

set(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

message(STATUS "project source directory -> ${PartMod_SOURCE_DIR}")
message(STATUS "project build directory -> ${PartMod_BINARY_DIR}")

include(GNUInstallDirs)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})


# Default is release build
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release)
endif ()

if (MSVC)
    add_definitions (-D_USE_MATH_DEFINES)
    set (CMAKE_REQUIRED_DEFINITIONS "-D_USE_MATH_DEFINES")
    # I'm not sure we should block these warnings or not:
    # C4800 Conversion from int to bool performance warning (we should probably be using bool in most of these places though).
    # C4127 Conditional expression is constant: Most of these uses are intentional (debugging code etc.)
    # C4267 size_t to int conversion. Most of these are indices, but many places can be replaced by a for-each loop as well.
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W4 -wd4786 -wd4100 -wd4996 -wd4706 -wd4512 -wd4251 -wd4800")# todo: find how to set multiprocessor compilation /MP, debug information format /Zi etc...
    #set(CMAKE_CXX_FLAGS_RELEASE "/MT")
	  #set(CMAKE_CXX_FLAGS_DEBUG "/MTd /Zi")
else ()
    #set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fPIC") # -Wextra
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fPIC -Wall")# -Wextra -Wshadow -fstack-protector-strong -fsanitize=address -fsanitize=undefined -fsanitize=leak -fno-omit-frame-pointer")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -fPIC") # Wextra
endif ()

check_cxx_compiler_flag("-march=native" _march_native_works)
check_cxx_compiler_flag("-xHost" _xhost_works)

set(_CXX_FLAGS)
if(_march_native_works)
  message(STATUS "Using processor's vector instructions (-march=native compiler flag set)")
  set(_CXX_FLAGS "-march=native")
elseif(_xhost_works)
  message(STATUS "Using processor's vector instructions (-xHost compiler flag set)")
  set(_CXX_FLAGS "-xHost")
else()
  message(STATUS "No suitable compiler flag found for vectorization")
endif()

set (CMAKE_CXX_STANDARD 17)

#set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
#set (CMAKE_CXX_FLAGS_PROFILING "${CMAKE_CXX_FLAGS_PROFILING} -O2 -pg")

# detect
find_package(PythonInterp REQUIRED)
find_package (Eigen3 3.3 NO_MODULE)
find_package(OpenMP REQUIRED)

if(TARGET Eigen3::Eigen)
  message(STATUS "Eigen3 v${EIGEN3_VERSION_STRING} found in ${EIGEN3_INCLUDE_DIR}")
  include_directories(${EIGEN3_INCLUDE_DIR})
else()
  message(STATUS "Eigen3 not found. Version from the thirdparty directory will be used.")
  add_custom_target(unpack-eigen
  ALL
  COMMAND
    ${CMAKE_COMMAND} -E tar xzf ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/eigen-3.3.7.tar.gz
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT
    "Unpacking Eigen3 in ${CMAKE_CURRENT_BINARY_DIR}/eigen-3.3.7"
  )
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/eigen-3.3.7)
endif()

#######################################################################
######################## Options ######################################
#######################################################################
option (USE_PREPROCESSOR "Enable preprocessor" ON)
option (USE_POSTPROCESSOR "Enable postprocessor" ON)
option (USE_SOLVER "Enable solver" ON)
option (USE_DOCS "Enable documentation" ON)

include_directories ("${PartMod_BINARY_DIR}")
set (LIB_TYPE OBJECT)

if (USE_SOLVER)
  include_directories ("${PartMod_SOURCE_DIR}/src/solver/src")
  add_subdirectory ("${PartMod_SOURCE_DIR}/src/solver/src")
  #list (APPEND MODULE_LIST "solver")
  list (APPEND LIBS $<TARGET_OBJECTS:solver>)
  target_compile_options(solver PRIVATE "${OpenMP_CXX_FLAGS}")
endif ()


#file(GLOB SOURCE "${PartMod_SOURCE_DIR}/src/solver/src/*.h" "${PartMod_SOURCE_DIR}/src/*.cpp")

#add_executable(DiscreteModel ${SOURCE})
#add_executable(DiscreteModel ${PartMod_SOURCE_DIR}/src/solver/src/main.cpp ${LIBS})
#add_executable(DiscreteModel ${PartMod_SOURCE_DIR}/src/solver/src/main.cpp ${solver_files})

add_executable(DiscreteModel 
                ${PartMod_SOURCE_DIR}/src/solver/src/main.cpp
                ${PartMod_SOURCE_DIR}/src/solver/src/version.cpp 
                ${LIBS}
                )
#target_link_libraries (DiscreteModel libDiscreteModel Eigen3::Eigen)
#find_package(Boost ${Boost_MINIMUM_REQUIRED} QUIET COMPONENTS "${BOOST_COMPONENTS_REQUIRED}")

if (CMAKE_CXX_COMPILER_ID STREQUAL GNU)
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
        set(CXX_FILESYSTEM_LIBRARIES stdc++fs) # Boost::filesystem
    endif()
else()
    set(CXX_FILESYSTEM_LIBRARIES)
endif()

if (NOT TARGET Eigen3::Eigen)
  add_dependencies(solver unpack-eigen)
endif()

target_link_libraries(DiscreteModel
                      PRIVATE
                      ${CXX_FILESYSTEM_LIBRARIES}
                      #Eigen3::Eigen
                     )
if(OpenMP_CXX_FOUND)
    target_link_libraries(DiscreteModel PUBLIC OpenMP::OpenMP_CXX)
endif()
target_compile_options(DiscreteModel
  PRIVATE
    ${_CXX_FLAGS}
)

install (
  TARGETS
    DiscreteModel
  DESTINATION
    ${CMAKE_INSTALL_BINDIR})

find_package(BLAS)

if(BLAS_FOUND)
  message(STATUS "Eigen will use some subroutines from BLAS.")
  message(STATUS "See: http://eigen.tuxfamily.org/dox-devel/TopicUsingBlasLapack.html")
  target_compile_definitions(DiscreteModel
    PRIVATE
      EIGEN_USE_BLAS
    )
  target_link_libraries(DiscreteModel
    PUBLIC
      ${BLAS_LIBRARIES}
    )
else()
  message(STATUS "BLAS not found. Using Eigen own functions")
endif()

find_package(LAPACK)
find_package(LAPACKE)

if(LAPACK_FOUND AND LAPACKE_FOUND)
  message(STATUS "Eigen will use some subroutines from LAPACK.")
  message(STATUS "See: http://eigen.tuxfamily.org/dox-devel/TopicUsingBlasLapack.html")
  target_compile_definitions(DiscreteModel
    PRIVATE
      EIGEN_USE_LAPACKE
    )
  target_link_libraries(DiscreteModel
    PUBLIC
      ${LAPACK_LIBRARIES}
      ${LAPACKE_LIBRARIES}
    )
else()
  message(STATUS "LAPACK not found. Using Eigen own functions")
endif()



set_target_properties(
  DiscreteModel
  PROPERTIES
    COMPILE_FLAGS
      ${OpenMP_CXX_FLAGS}
    LINK_FLAGS
      ${OpenMP_CXX_FLAGS}
  )
  
  
#######################################################################
######################## VERSION ######################################
#######################################################################

# needs to find the generated header file
target_include_directories(DiscreteModel
  PRIVATE
    ${PartMod_BINARY_DIR}/generated
  )
  
add_custom_command(
  OUTPUT
    ${PartMod_BINARY_DIR}/generated/version.h
    ALL
  COMMAND
    ${CMAKE_COMMAND} -D TARGET_DIR=${PartMod_BINARY_DIR} -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -P ${PartMod_SOURCE_DIR}/version.cmake
  WORKING_DIRECTORY
    ${PartMod_SOURCE_DIR}
  )

# rebuild version.hpp every time
add_custom_target(
  get_git_hash
  ALL
  DEPENDS
    ${PartMod_BINARY_DIR}/generated/version.h
  )

# version.hpp has to be generated
# before we start building DiscreteModel
add_dependencies(DiscreteModel get_git_hash)

#######################################################################
######################## Documentation ################################
#######################################################################
find_package (Doxygen)
if (DOXYGEN_FOUND)
    set (DOXYGEN_INPUT ${PartMod_BINARY_DIR}/Doxyfile)
    set (DOXYGEN_OUTPUT ${PartMod_BINARY_DIR}/partmod_docs/)
    #file (MAKE_DIRECTORY ${DOXYGEN_OUTPUT})

    #configure_file (
    #    "${PartMod_SOURCE_DIR}/doc/refman/footer.html.in"
  #    "${PartMod_BINARY_DIR}/footer.html")
    configure_file (
        "${PartMod_SOURCE_DIR}/docs/Doxyfile"
      "${PartMod_BINARY_DIR}/Doxyfile")


    add_custom_target (
        docs
        COMMAND ${CMAKE_COMMAND} -E echo_append "Building Documentation..."
        COMMAND ${CMAKE_COMMAND} -E make_directory ${DOXYGEN_OUTPUT}
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}
        COMMAND ${CMAKE_COMMAND} -E echo "Done."
        WORKING_DIRECTORY ${PartMod_BINARY_DIR}
        DEPENDS)

    get_filename_component(PARTMOD_SITE_SOURCE ../partmod-site ABSOLUTE $(PartMod_SOURCE_DIR))
    message(STATUS "PARTMOD_SITE_SOURCE -> ${PARTMOD_SITE_SOURCE}")
    find_package(Git)
    IF(EXISTS "${PARTMOD_SITE_SOURCE}" AND IS_DIRECTORY "${PARTMOD_SITE_SOURCE}" AND GIT_FOUND)
        add_custom_target (
            docs_publish
            COMMAND ${CMAKE_COMMAND} -E echo_append "Deploying Documentation..."
            COMMAND ${GIT_EXECUTABLE} merge master
            COMMAND ${CMAKE_COMMAND} -E remove_directory ${PARTMOD_SITE_SOURCE}/public
            COMMAND ${CMAKE_COMMAND} -E copy_directory ${DOXYGEN_OUTPUT}/html ${PARTMOD_SITE_SOURCE}/public
            COMMAND ${GIT_EXECUTABLE} add -A ${PARTMOD_SITE_SOURCE}
            COMMAND ${GIT_EXECUTABLE} commit -m "documentation updated"
            COMMAND ${GIT_EXECUTABLE} push origin master
            COMMAND ${CMAKE_COMMAND} -E echo "Done."
            WORKING_DIRECTORY ${PARTMOD_SITE_SOURCE}
            DEPENDS)
    endif ()

    add_custom_target (
        docs_publish_rsync
        COMMAND ${CMAKE_COMMAND} -E echo_append "Deploying Documentation by RSync..."
        #COMMAND ${GIT_EXECUTABLE} merge master
        #COMMAND ${CMAKE_COMMAND} -E remove_directory ${PARTMOD_SITE_SOURCE}/public
        #COMMAND ${CMAKE_COMMAND} -E copy_directory ${DOXYGEN_OUTPUT}/html ${PARTMOD_SITE_SOURCE}/public
        #COMMAND ${GIT_EXECUTABLE} add -A ${PARTMOD_SITE_SOURCE}
        #COMMAND ${GIT_EXECUTABLE} commit -m "documentation updated"
        #COMMAND ${GIT_EXECUTABLE} push
        COMMAND ${CMAKE_COMMAND} -E echo "Done."
        WORKING_DIRECTORY ${PARTMOD_SITE_SOURCE}
        DEPENDS)
endif (DOXYGEN_FOUND)


#######################################################################
######################## Packing ######################################
#######################################################################
# https://stackoverflow.com/questions/7050997/zip-files-using-cmake
### versions
set (CPACK_PACKAGE_VERSION_MAJOR ${PartMod_VERSION_MAJOR})
set (CPACK_PACKAGE_VERSION_MINOR ${PartMod_VERSION_MINOR})
set (CPACK_PACKAGE_VERSION_PATCH ${PartMod_VERSION_PATCH})
set (CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}")

### general settings
set (CPACK_PACKAGE_NAME "partmod")
#set (CPACK_PACKAGE_VENDOR "OOFEM development team")
#set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Object Oriented Finite Element")
set (CPACK_RESOURCE_FILE_README "${PartMod_SOURCE_DIR}/README.md")
#set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING.LGPLv2.1")
set (CPACK_PACKAGE_INSTALL_DIRECTORY ${CPACK_PACKAGE_NAME})
set (CPACK_PACKAGE_EXECUTABLES "DiscreteModel")
set (CPACK_PACKAGE_FILE_NAME ${CPACK_PACKAGE_NAME}_${CPACK_PACKAGE_VERSION}_${CMAKE_SYSTEM_PROCESSOR})

### default generators
if (WIN32)
    set(CPACK_GENERATOR "ZIP")
else ()
    set (CPACK_GENERATOR "TGZ;DEB")
endif ()

### source package (type "make source_package" to produce a zip)
set (CPACK_SOURCE_GENERATOR "ZIP")
set (CPACK_SOURCE_PACKAGE_FILE_NAME "partmod-${CPACK_PACKAGE_VERSION}")
set (CPACK_SOURCE_IGNORE_FILES "~$;/build/;tags;cscope.*;.*\\\\.out$;\\\\.out\\\\.;/\\\\..*;\\\\.kdev4$;do_release;release_filter\\\\.pl") # Note: All hidden files are ignored

### deb package options
set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Name <e-mail>")
set (CPACK_DEBIAN_PACKAGE_SECTION "Mathematics")
# autogenerate dependency information
set (CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)

set(CPACK_DEBIAN_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION}.${CPACK_PACKAGE_VERSION_PATCH}+sid1")

include (CPack)



#######################################################################
######################## Testing ######################################
#######################################################################
find_program(MEMORYCHECK_COMMAND NAMES valgrind)
set(MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full")

# add memcheck test action
include(CTest)

# define tests
enable_testing()

# BENCHMARKS
include(src/benchmark/CMakeLists.cmake)
# TESTS
include(src/tests/CMakeLists.cmake)
