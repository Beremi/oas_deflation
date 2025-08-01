cmake_minimum_required(VERSION 3.15.0 FATAL_ERROR)

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


project(OAS C CXX)
set (OAS_VERSION_MAJOR 0)
set (OAS_VERSION_MINOR 1)
set (OAS_VERSION_PATCH 0) # Unused
set (OAS_HOMEPAGE "https://kelidas.gitlab.io/OAS-site")
#set (OAS_COPYRIGHT "Copyright (C) 2019")
#set (OAS_EMAIL "elias.j@fce.vutbr.cz")

set(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

message(STATUS "project source directory -> ${OAS_SOURCE_DIR}")
message(STATUS "project build directory -> ${OAS_BINARY_DIR}")

include(GNUInstallDirs)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})


# Default is release build
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release)
endif ()

#######################################################################
######################## Options ######################################
#######################################################################
option (USE_PREPROCESSOR "Enable preprocessor" ON)
option (USE_POSTPROCESSOR "Enable postprocessor" ON)
option (USE_SOLVER "Enable solver" ON)
option (USE_ML_TORCH "Enable ML Torch" OFF)
option (ARCHIVE_BINARIES "Archive binary files" OFF)
option (USE_DOCS "Enable documentation" ON)
option (USE_VTK "Enable VTK (for exporting binary VTU-files)" ON)
option (USE_EXPRTK "Enable EXPRTK" OFF)
option (USE_SANITIZERS "Enable sanitizers" OFF)
option (CHECK_CONVERSION "Enable conversion warnings" OFF)

if (MSVC)
    add_definitions (-D_USE_MATH_DEFINES)
    # Build with multiple processes
    add_definitions(/MP)
    # Allow big object
    add_definitions(/bigobj)
    set (CMAKE_REQUIRED_DEFINITIONS "-D_USE_MATH_DEFINES")
    # I'm not sure we should block these warnings or not:
    # C4800 Conversion from int to bool performance warning (we should probably be using bool in most of these places though).
    # C4127 Conditional expression is constant: Most of these uses are intentional (debugging code etc.)
    # C4267 size_t to int conversion. Most of these are indices, but many places can be replaced by a for-each loop as well.
    # /wd4018 # 'expression' : signed/unsigned mismatch
    # /wd4065 # switch statement contains 'default' but no 'case' labels
    # /wd4146 # unary minus operator applied to unsigned type, result still unsigned
    # /wd4244 # 'conversion' conversion from 'type1' to 'type2', possible loss of data
    # /wd4251 # 'identifier' : class 'type' needs to have dll-interface to be used by clients of class 'type2'
    # /wd4267 # 'var' : conversion from 'size_t' to 'type', possible loss of data
    # /wd4305 # 'identifier' : truncation from 'type1' to 'type2'
    # /wd4307 # 'operator' : integral constant overflow
    # /wd4309 # 'conversion' : truncation of constant value
    # /wd4334 # 'operator' : result of 32-bit shift implicitly converted to 64 bits (was 64-bit shift intended?)
    # /wd4355 # 'this' : used in base member initializer list
    # /wd4506 # no definition for inline function 'function'
    # /wd4800 # 'type' : forcing value to bool 'true' or 'false' (performance warning)
    # /wd4996 # The compiler encountered a deprecated declaration.
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W4 -wd4786 -wd4100 -wd4996 -wd4706 -wd4512 -wd4251 -wd4800")
    #set(CMAKE_CXX_FLAGS_RELEASE "/MT")
	  #set(CMAKE_CXX_FLAGS_DEBUG "/MTd /Zi") # todo: find how to set multiprocessor compilation /MP, debug information format /Zi etc...
else ()
    #-ffast-math disables all NaN and Inf checks
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fPIC -Wall -Wextra  -Wshadow -ffast-math -fno-math-errno -march=native -DNDEBUG") # -pedantic-errorsc -ggdb  -fopenmp
    #set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fPIC -Wall -Wextra  -Wshadow -fno-math-errno -march=native -DNDEBUG") # -pedantic-errorsc -ggdb  -fopenmp
    if (CHECK_CONVERSION)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsign-conversion -Wconversion")
    endif()
    if (USE_SANITIZERS)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fstack-protector-strong -fsanitize=address -fsanitize=undefined -fsanitize=leak -fno-omit-frame-pointer") #
    endif ()
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -fPIC") #  -Wextra
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

set (CMAKE_CXX_STANDARD 20)

#set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
#set (CMAKE_CXX_FLAGS_PROFILING "${CMAKE_CXX_FLAGS_PROFILING} -O2 -pg")
set(TARGETNAME OAS)

# detect
#find_package (Python3 COMPONENTS Interpreter)
find_package (Eigen3 3.3.9 NO_MODULE)
find_package(OpenMP REQUIRED)

if (USE_MPI)
    if (MPI_DIR)
        include_directories ("${MPI_DIR}/include")
        find_library (MPI_LIBRARY NAMES mpi PATHS ${MPI_DIR}/lib64 ${MPI_DIR}/lib)
        find_library (MPI_CXX_LIBRARY NAMES mpi_cxx PATHS ${MPI_DIR}/lib64 ${MPI_DIR}/lib)
        list (APPEND EXT_LIBS ${MPI_LIBRARY})
        if (MPI_CXX_LIBRARY)
            list (APPEND test_eigen_LIBRARIES ${MPI_CXX_LIBRARY})
        endif ()
    else ()
        include (FindMPI)
        find_package (MPI REQUIRED)
        set (CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILES_FLAGS} ${MPI_COMPILE_FLAGS})
        set (CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
        add_definitions (-DMPICH_IGNORE_CXX_SEEK)
        include_directories (${MPI_INCLUDE_PATH})
        list (APPEND test_eigen_LIBRARIES ${MPI_LIBRARIES})
    endif (MPI_DIR)
    list (APPEND MODULE_LIST "MPI")
    # Compile and link flags?
endif ()

if(TARGET Eigen3::Eigen)
  message(STATUS "Eigen3 v${EIGEN3_VERSION_STRING} found in ${EIGEN3_INCLUDE_DIR}")
  include_directories(${EIGEN3_INCLUDE_DIR})
else()
  message(STATUS "Eigen3 not found. Version from the thirdparty directory will be used.")
  add_custom_target(unpack-eigen
  ALL
  COMMAND
    ${CMAKE_COMMAND} -E tar xzf ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/eigen-3.4.0.tar.gz
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT
    "Unpacking Eigen3 in ${CMAKE_CURRENT_BINARY_DIR}/eigen-3.4.0"
  )
  include_directories(SYSTEM ${CMAKE_CURRENT_BINARY_DIR}/eigen-3.4.0)
endif()

add_custom_target(unpack-spectra
  ALL
  COMMAND
    ${CMAKE_COMMAND} -E tar xzf ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/spectra-1.1.0.tar.gz
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT
    "Unpacking Spectra in ${CMAKE_CURRENT_BINARY_DIR}/spectra-1.1.0"
  )
include_directories(SYSTEM ${CMAKE_CURRENT_BINARY_DIR}/spectra-1.1.0/include)

include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/exprtk)

include_directories ("${OAS_BINARY_DIR}")
set (LIB_TYPE OBJECT) # set (LIB_TYPE SHARED)


include_directories ("${OAS_SOURCE_DIR}/src/solver/src")
add_subdirectory ("${OAS_SOURCE_DIR}/src/solver/src")

if (USE_SOLVER)
  #list (APPEND MODULE_LIST "solver")
  list (APPEND LIBS $<TARGET_OBJECTS:solver>)
  target_compile_options(solver PRIVATE "${OpenMP_CXX_FLAGS}")
endif ()

if (USE_ML_TORCH)
  find_package(Torch REQUIRED PATHS "${CMAKE_BINARY_DIR}/libtorch")
  if (NOT TORCH_FOUND)
    message("-- Torch not found")
  else ()
    message("-- Torch (${Torch_VERSION}) found in ${Torch_DIR}.")
    list (APPEND LIBS $<TARGET_OBJECTS:ml_torch> )
    target_compile_options(solver PRIVATE "${OpenMP_CXX_FLAGS}")
    list (APPEND EXT_LIBS ${TORCH_LIBRARIES})
    add_definitions(-DML_TORCH_FOUND)
  endif ()
endif ()

#file(GLOB SOURCE "${OAS_SOURCE_DIR}/src/solver/src/*.h" "${OAS_SOURCE_DIR}/src/*.cpp")

#add_executable(${TARGETNAME} ${SOURCE})
#add_executable(${TARGETNAME} ${OAS_SOURCE_DIR}/src/solver/src/main.cpp ${LIBS})
#add_executable(${TARGETNAME} ${OAS_SOURCE_DIR}/src/solver/src/main.cpp ${solver_files})
add_executable(${TARGETNAME}
                ${OAS_SOURCE_DIR}/src/solver/src/main.cpp
                ${OAS_SOURCE_DIR}/src/solver/src/version.cpp
                ${LIBS}
                )

#add_library (l${TARGETNAME}
#        SHARED
#        ${LIBS}
#        )
#target_link_libraries (lib${TARGETNAME} ${EXT_LIBS})

if (USE_VTK)
    find_package (VTK REQUIRED NO_MODULE)
    if (NOT VTK_FOUND)
        message("-- Skipping WriteVTU: ${VTK_NOT_FOUND_MESSAGE}")
    endif()
    if (VTK_VERSION VERSION_LESS "8.90.0")
        message("-- VTK is older than 8.90.0 (${VTK_VERSION})")
        # old system
        include(${VTK_USE_FILE})
        set (CMAKE_LIBRARY_DIR ${VTK_DIR})
        target_compile_definitions(solver PUBLIC __VTK_MODULE)
        list (APPEND EXT_LIBS ${VTK_LIBRARIES})
        list (APPEND MODULE_LIST "VTK")
        add_definitions(-DVTK_FOUND)
    else ()
        message("-- VTK found: ${VTK_VERSION}")
        target_compile_definitions(solver PUBLIC __VTK_MODULE)
        list (APPEND EXT_LIBS ${VTK_LIBRARIES})
        list (APPEND MODULE_LIST "VTK")
        # vtk_module_autoinit is needed
        vtk_module_autoinit(
            TARGETS ${TARGETNAME} solver
            MODULES ${VTK_LIBRARIES}
            )
        add_definitions(-DVTK_FOUND)
    endif ()
endif ()

if (USE_EXPRTK)
    message("-- Adding Exprtk")
    target_compile_definitions (solver PUBLIC -D__EXPRTK_MODULE)
endif ()

# precompiled headers - no speedup?
#target_precompile_headers(solver
#                PUBLIC
#                  ${OAS_SOURCE_DIR}/src/solver/src/exprtk.hpp)

#target_link_libraries (${TARGETNAME} libOAS Eigen3::Eigen)
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
add_dependencies(solver unpack-spectra)

target_link_libraries(solver
                      PRIVATE
                      ${CXX_FILESYSTEM_LIBRARIES}
                      ${EXT_LIBS}
                      #Eigen3::Eigen
                     )

if (USE_ML_TORCH)
  target_compile_definitions (solver PUBLIC -DML_TORCH_FOUND)
  target_link_libraries(ml_torch
                        PRIVATE
                        ${CXX_FILESYSTEM_LIBRARIES}
                        ${EXT_LIBS}
                        #Eigen3::Eigen
                      )
endif ()

target_link_libraries(${TARGETNAME}
                      PRIVATE
                      ${CXX_FILESYSTEM_LIBRARIES}
                      ${EXT_LIBS}
                      #Eigen3::Eigen
                     )

if(OpenMP_CXX_FOUND)
    target_link_libraries(solver PUBLIC OpenMP::OpenMP_CXX)
endif()
target_compile_options(${TARGETNAME}
  PRIVATE
    ${_CXX_FLAGS}
)

install (
  TARGETS
    ${TARGETNAME}
  DESTINATION
    ${CMAKE_INSTALL_BINDIR})

find_package(BLAS)

if(BLAS_FOUND)
  message(STATUS "Eigen will use some subroutines from BLAS.")
  message(STATUS "See: http://eigen.tuxfamily.org/dox-devel/TopicUsingBlasLapack.html")
  target_compile_definitions(${TARGETNAME}
    PRIVATE
      EIGEN_USE_BLAS
    )
  target_link_libraries(${TARGETNAME}
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
  target_compile_definitions(${TARGETNAME}
    PRIVATE
      EIGEN_USE_LAPACKE
    )
  target_link_libraries(${TARGETNAME}
    PUBLIC
      ${LAPACK_LIBRARIES}
      ${LAPACKE_LIBRARIES}
    )
else()
  message(STATUS "LAPACK not found. Using Eigen own functions")
endif()



set_target_properties(
  ${TARGETNAME}
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
target_include_directories(${TARGETNAME}
  PRIVATE
    ${OAS_BINARY_DIR}/generated
  )

add_custom_command(
  OUTPUT
    ${OAS_BINARY_DIR}/generated/version.h
    ALL
  COMMAND
    ${CMAKE_COMMAND} -D TARGET_NAME=${TARGETNAME} -D TARGET_DIR=${OAS_BINARY_DIR} -D SOURCE_DIR=${OAS_SOURCE_DIR} -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                     -D VTK_VERSION=${VTK_VERSION} -P ${OAS_SOURCE_DIR}/version.cmake
  WORKING_DIRECTORY
    ${OAS_SOURCE_DIR}
  )

# rebuild version.hpp every time
add_custom_target(
  get_git_hash
  ALL
  DEPENDS
    ${OAS_BINARY_DIR}/generated/version.h
  )

# version.hpp has to be generated
# before we start building OAS
add_dependencies(${TARGETNAME} get_git_hash)


# Copy binary
if (ARCHIVE_BINARIES)
  #add_custom_command(
  #        TARGET ${TARGETNAME} POST_BUILD
  #        COMMAND ${CMAKE_COMMAND} -E copy
  #                ${CMAKE_CURRENT_BINARY_DIR}/bin/OAS
  #                ${CMAKE_CURRENT_BINARY_DIR}/bin/OAS_${GIT_HASH}_${DATE_STRING})
  add_custom_command(
          TARGET ${TARGETNAME} POST_BUILD
          COMMAND bash ${OAS_SOURCE_DIR}/archive.sh
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif ()

#######################################################################
######################## Documentation ################################
#######################################################################
find_package (Doxygen)
if (DOXYGEN_FOUND)
    set (DOXYGEN_INPUT ${OAS_BINARY_DIR}/Doxyfile)
    set (DOXYGEN_OUTPUT ${OAS_BINARY_DIR}/OAS_docs/)
    #file (MAKE_DIRECTORY ${DOXYGEN_OUTPUT})

    #configure_file (
    #    "${OAS_SOURCE_DIR}/doc/refman/footer.html.in"
  #    "${OAS_BINARY_DIR}/footer.html")
    configure_file (
        "${OAS_SOURCE_DIR}/docs/Doxyfile"
      "${OAS_BINARY_DIR}/Doxyfile")


    add_custom_target (
        docs
        COMMAND ${CMAKE_COMMAND} -E echo_append "Building Documentation..."
        COMMAND ${CMAKE_COMMAND} -E make_directory ${DOXYGEN_OUTPUT}
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}
        COMMAND ${CMAKE_COMMAND} -E echo "Done."
        WORKING_DIRECTORY ${OAS_BINARY_DIR}
        DEPENDS)

    get_filename_component(OAS_SITE_SOURCE ../OAS-site ABSOLUTE $(OAS_SOURCE_DIR))
    message(STATUS "OAS_SITE_SOURCE -> ${OAS_SITE_SOURCE}")
    find_package(Git)
    IF(EXISTS "${OAS_SITE_SOURCE}" AND IS_DIRECTORY "${OAS_SITE_SOURCE}" AND GIT_FOUND)
        add_custom_target (
            docs_publish
            COMMAND ${CMAKE_COMMAND} -E echo_append "Deploying Documentation..."
            COMMAND ${GIT_EXECUTABLE} merge master
            COMMAND ${CMAKE_COMMAND} -E remove_directory ${OAS_SITE_SOURCE}/public
            COMMAND ${CMAKE_COMMAND} -E copy_directory ${DOXYGEN_OUTPUT}/html ${OAS_SITE_SOURCE}/public
            COMMAND ${GIT_EXECUTABLE} add -A ${OAS_SITE_SOURCE}
            COMMAND ${GIT_EXECUTABLE} commit -m "documentation updated"
            COMMAND ${GIT_EXECUTABLE} push origin master
            COMMAND ${CMAKE_COMMAND} -E echo "Done."
            WORKING_DIRECTORY ${OAS_SITE_SOURCE}
            DEPENDS)
    endif ()

    add_custom_target (
        docs_publish_rsync
        COMMAND ${CMAKE_COMMAND} -E echo_append "Deploying Documentation by RSync..."
        #COMMAND ${GIT_EXECUTABLE} merge master
        #COMMAND ${CMAKE_COMMAND} -E remove_directory ${OAS_SITE_SOURCE}/public
        #COMMAND ${CMAKE_COMMAND} -E copy_directory ${DOXYGEN_OUTPUT}/html ${OAS_SITE_SOURCE}/public
        #COMMAND ${GIT_EXECUTABLE} add -A ${OAS_SITE_SOURCE}
        #COMMAND ${GIT_EXECUTABLE} commit -m "documentation updated"
        #COMMAND ${GIT_EXECUTABLE} push
        COMMAND ${CMAKE_COMMAND} -E echo "Done."
        WORKING_DIRECTORY ${OAS_SITE_SOURCE}
        DEPENDS)
endif (DOXYGEN_FOUND)


#######################################################################
######################## Packing ######################################
#######################################################################
# https://stackoverflow.com/questions/7050997/zip-files-using-cmake
### versions
set (CPACK_PACKAGE_VERSION_MAJOR ${OAS_VERSION_MAJOR})
set (CPACK_PACKAGE_VERSION_MINOR ${OAS_VERSION_MINOR})
set (CPACK_PACKAGE_VERSION_PATCH ${OAS_VERSION_PATCH})
set (CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}")

### general settings
set (CPACK_PACKAGE_NAME "OAS")
#set (CPACK_PACKAGE_VENDOR "OOFEM development team")
#set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Object Oriented Finite Element")
set (CPACK_RESOURCE_FILE_README "${OAS_SOURCE_DIR}/README.md")
#set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING.LGPLv2.1")
set (CPACK_PACKAGE_INSTALL_DIRECTORY ${CPACK_PACKAGE_NAME})
set (CPACK_PACKAGE_EXECUTABLES ${TARGETNAME})
set (CPACK_PACKAGE_FILE_NAME ${CPACK_PACKAGE_NAME}_${CPACK_PACKAGE_VERSION}_${CMAKE_SYSTEM_PROCESSOR})

### default generators
if (WIN32)
    set(CPACK_GENERATOR "ZIP")
else ()
    set (CPACK_GENERATOR "TGZ;DEB")
endif ()

### source package (type "make source_package" to produce a zip)
set (CPACK_SOURCE_GENERATOR "ZIP")
set (CPACK_SOURCE_PACKAGE_FILE_NAME "OAS-${CPACK_PACKAGE_VERSION}")
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
