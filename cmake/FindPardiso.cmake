# FindPardiso.cmake - Find Intel MKL Pardiso library
# This module defines:
#  Pardiso_FOUND - System has Pardiso (Intel MKL)
#  Pardiso_INCLUDE_DIRS - Pardiso/MKL include directories
#  Pardiso_LIBRARIES - Libraries needed to use Pardiso
#  Pardiso_VERSION - The version of MKL found

# Pardiso is part of Intel MKL, so we need to find MKL
# Try to find MKL using pkg-config first
find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
    pkg_check_modules(PC_MKL QUIET mkl)
endif()

# Find the MKL include directory
find_path(Pardiso_INCLUDE_DIR
    NAMES mkl_pardiso.h mkl.h
    PATHS
        ${PC_MKL_INCLUDE_DIRS}
        $ENV{MKLROOT}/include
        $ENV{MKL_DIR}/include
        $ENV{CONDA_PREFIX}/include
        /opt/intel/mkl/include
        /opt/intel/oneapi/mkl/latest/include
        /usr/include/mkl
        /usr/local/include
        /usr/include
    PATH_SUFFIXES
        mkl
)

# Find MKL libraries
# Pardiso requires: mkl_core, mkl_intel_lp64 (or ilp64), mkl_sequential (or intel_thread)
# For threaded version also need iomp5

# Try to find the core MKL library first
find_library(MKL_CORE_LIBRARY
    NAMES mkl_core
    PATHS
        ${PC_MKL_LIBRARY_DIRS}
        $ENV{MKLROOT}/lib/intel64
        $ENV{MKLROOT}/lib
        $ENV{MKL_DIR}/lib/intel64
        $ENV{MKL_DIR}/lib
        $ENV{CONDA_PREFIX}/lib
        /opt/intel/mkl/lib/intel64
        /opt/intel/oneapi/mkl/latest/lib/intel64
        /usr/lib/x86_64-linux-gnu
        /usr/lib
        /usr/local/lib
)

# Find interface layer (LP64 for 32-bit integers)
find_library(MKL_INTERFACE_LIBRARY
    NAMES mkl_intel_lp64
    PATHS
        ${PC_MKL_LIBRARY_DIRS}
        $ENV{MKLROOT}/lib/intel64
        $ENV{MKLROOT}/lib
        $ENV{MKL_DIR}/lib/intel64
        $ENV{MKL_DIR}/lib
        $ENV{CONDA_PREFIX}/lib
        /opt/intel/mkl/lib/intel64
        /opt/intel/oneapi/mkl/latest/lib/intel64
        /usr/lib/x86_64-linux-gnu
        /usr/lib
        /usr/local/lib
)

# Find threading layer - IMPORTANT: Prefer threaded versions for parallel Pardiso!
# Priority: gnu_thread > intel_thread > sequential (single-threaded)
find_library(MKL_THREADING_LIBRARY
    NAMES mkl_gnu_thread mkl_intel_thread mkl_sequential
    PATHS
        ${PC_MKL_LIBRARY_DIRS}
        $ENV{MKLROOT}/lib/intel64
        $ENV{MKLROOT}/lib
        $ENV{MKL_DIR}/lib/intel64
        $ENV{MKL_DIR}/lib
        $ENV{CONDA_PREFIX}/lib
        /opt/intel/mkl/lib/intel64
        /opt/intel/oneapi/mkl/latest/lib/intel64
        /usr/lib/x86_64-linux-gnu
        /usr/lib
        /usr/local/lib
)

# Try to extract MKL version
if(Pardiso_INCLUDE_DIR)
    if(EXISTS "${Pardiso_INCLUDE_DIR}/mkl_version.h")
        file(STRINGS "${Pardiso_INCLUDE_DIR}/mkl_version.h" MKL_VERSION_STRING
             REGEX "#define __INTEL_MKL(_MINOR|_UPDATE)?__")
        
        if(MKL_VERSION_STRING)
            string(REGEX REPLACE ".*#define __INTEL_MKL__[ \t]+([0-9]+).*" "\\1" 
                   MKL_MAJOR_VERSION "${MKL_VERSION_STRING}")
            string(REGEX REPLACE ".*#define __INTEL_MKL_MINOR__[ \t]+([0-9]+).*" "\\1" 
                   MKL_MINOR_VERSION "${MKL_VERSION_STRING}")
            string(REGEX REPLACE ".*#define __INTEL_MKL_UPDATE__[ \t]+([0-9]+).*" "\\1" 
                   MKL_UPDATE_VERSION "${MKL_VERSION_STRING}")
            
            set(Pardiso_VERSION "${MKL_MAJOR_VERSION}.${MKL_MINOR_VERSION}.${MKL_UPDATE_VERSION}")
        endif()
    endif()
endif()

# Handle the QUIETLY and REQUIRED arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Pardiso
    FOUND_VAR Pardiso_FOUND
    REQUIRED_VARS MKL_CORE_LIBRARY MKL_INTERFACE_LIBRARY MKL_THREADING_LIBRARY Pardiso_INCLUDE_DIR
    VERSION_VAR Pardiso_VERSION
)

if(Pardiso_FOUND)
    # Set the libraries in the correct order for linking
    # The order matters for MKL: interface -> threading -> core
    set(Pardiso_LIBRARIES 
        ${MKL_INTERFACE_LIBRARY}
        ${MKL_THREADING_LIBRARY}
        ${MKL_CORE_LIBRARY}
    )
    set(Pardiso_INCLUDE_DIRS ${Pardiso_INCLUDE_DIR})
    
    # Create imported target
    if(NOT TARGET Pardiso::Pardiso)
        add_library(Pardiso::Pardiso INTERFACE IMPORTED)
        set_target_properties(Pardiso::Pardiso PROPERTIES
            INTERFACE_LINK_LIBRARIES "${Pardiso_LIBRARIES}"
            INTERFACE_INCLUDE_DIRECTORIES "${Pardiso_INCLUDE_DIR}"
        )
    endif()
    
    # Detect which threading layer was found and add appropriate OpenMP library
    get_filename_component(MKL_THREADING_NAME "${MKL_THREADING_LIBRARY}" NAME_WE)
    if(MKL_THREADING_NAME MATCHES "gnu_thread")
        # Using GNU OpenMP - need libgomp
        find_library(GOMP_LIBRARY 
            NAMES gomp
            PATHS
                /usr/lib/gcc/x86_64-linux-gnu/13
                /usr/lib/gcc/x86_64-linux-gnu/12
                /usr/lib/gcc/x86_64-linux-gnu/11
                /usr/lib/gcc/x86_64-linux-gnu/10
                /usr/lib/gcc/x86_64-linux-gnu/9
                /usr/lib/x86_64-linux-gnu
                /usr/lib
        )
        if(GOMP_LIBRARY)
            list(APPEND Pardiso_LIBRARIES ${GOMP_LIBRARY})
            message(STATUS "Pardiso: Using GNU threading with libgomp: ${GOMP_LIBRARY}")
        else()
            message(WARNING "Pardiso: Using GNU threading but libgomp not found - parallel execution may not work")
            message(WARNING "  Try: sudo apt-get install libgomp1")
        endif()
    elseif(MKL_THREADING_NAME MATCHES "intel_thread")
        # Using Intel OpenMP - need libiomp5
        find_library(IOMP5_LIBRARY NAMES iomp5)
        if(IOMP5_LIBRARY)
            list(APPEND Pardiso_LIBRARIES ${IOMP5_LIBRARY})
            message(STATUS "Pardiso: Using Intel threading (libiomp5)")
        else()
            message(WARNING "Pardiso: Using Intel threading but libiomp5 not found - parallel execution may not work")
        endif()
    elseif(MKL_THREADING_NAME MATCHES "sequential")
        message(WARNING "Pardiso: Using sequential (single-threaded) MKL - no parallel execution!")
        message(WARNING "  For parallel Pardiso, install libmkl-gnu-thread or libmkl-intel-thread")
    endif()
    
    # Note: Pardiso also needs libpthread and libm on Linux
    if(UNIX AND NOT APPLE)
        find_library(PTHREAD_LIBRARY pthread)
        find_library(M_LIBRARY m)
        find_library(DL_LIBRARY dl)
        if(PTHREAD_LIBRARY)
            list(APPEND Pardiso_LIBRARIES ${PTHREAD_LIBRARY})
        endif()
        if(M_LIBRARY)
            list(APPEND Pardiso_LIBRARIES ${M_LIBRARY})
        endif()
        if(DL_LIBRARY)
            list(APPEND Pardiso_LIBRARIES ${DL_LIBRARY})
        endif()
    endif()
endif()

mark_as_advanced(
    Pardiso_INCLUDE_DIR
    MKL_CORE_LIBRARY
    MKL_INTERFACE_LIBRARY
    MKL_THREADING_LIBRARY
)
