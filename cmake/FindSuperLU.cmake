# FindSuperLU.cmake - Find SuperLU library
# This module defines:
#  SuperLU_FOUND - System has SuperLU
#  SuperLU_INCLUDE_DIRS - SuperLU include directories
#  SuperLU_LIBRARIES - Libraries needed to use SuperLU
#  SuperLU_VERSION - The version of SuperLU found

# Try to find SuperLU using pkg-config first
find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
    pkg_check_modules(PC_SuperLU QUIET superlu)
endif()

# Find the include directory
find_path(SuperLU_INCLUDE_DIR
    NAMES supermatrix.h slu_ddefs.h
    PATHS
        ${PC_SuperLU_INCLUDE_DIRS}
        $ENV{SUPERLU_DIR}/include
        $ENV{CONDA_PREFIX}/include
        $ENV{CONDA_PREFIX}/include/superlu
        /usr/include/superlu
        /usr/local/include/superlu
        /usr/include
        /usr/local/include
    PATH_SUFFIXES
        superlu
        SuperLU
)

# Find the library
find_library(SuperLU_LIBRARY
    NAMES superlu superlu_5.0 superlu_5.1 superlu_5.2 superlu_5.3 superlu_6.0
    PATHS
        ${PC_SuperLU_LIBRARY_DIRS}
        $ENV{SUPERLU_DIR}/lib
        $ENV{CONDA_PREFIX}/lib
        /usr/lib
        /usr/local/lib
        /usr/lib/x86_64-linux-gnu
)

# Try to extract version from slu_ddefs.h
if(SuperLU_INCLUDE_DIR)
    file(STRINGS "${SuperLU_INCLUDE_DIR}/slu_ddefs.h" SuperLU_VERSION_STRING
         REGEX "#define SUPERLU_[A-Z]+_VERSION")
    
    if(SuperLU_VERSION_STRING)
        string(REGEX REPLACE ".*#define SUPERLU_MAJOR_VERSION[ \t]+([0-9]+).*" "\\1" 
               SuperLU_MAJOR_VERSION "${SuperLU_VERSION_STRING}")
        string(REGEX REPLACE ".*#define SUPERLU_MINOR_VERSION[ \t]+([0-9]+).*" "\\1" 
               SuperLU_MINOR_VERSION "${SuperLU_VERSION_STRING}")
        string(REGEX REPLACE ".*#define SUPERLU_PATCH_VERSION[ \t]+([0-9]+).*" "\\1" 
               SuperLU_PATCH_VERSION "${SuperLU_VERSION_STRING}")
        
        set(SuperLU_VERSION "${SuperLU_MAJOR_VERSION}.${SuperLU_MINOR_VERSION}.${SuperLU_PATCH_VERSION}")
    endif()
endif()

# Handle the QUIETLY and REQUIRED arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SuperLU
    FOUND_VAR SuperLU_FOUND
    REQUIRED_VARS SuperLU_LIBRARY SuperLU_INCLUDE_DIR
    VERSION_VAR SuperLU_VERSION
)

if(SuperLU_FOUND)
    set(SuperLU_LIBRARIES ${SuperLU_LIBRARY})
    set(SuperLU_INCLUDE_DIRS ${SuperLU_INCLUDE_DIR})
    
    if(NOT TARGET SuperLU::SuperLU)
        add_library(SuperLU::SuperLU UNKNOWN IMPORTED)
        set_target_properties(SuperLU::SuperLU PROPERTIES
            IMPORTED_LOCATION "${SuperLU_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${SuperLU_INCLUDE_DIR}"
        )
    endif()
endif()

mark_as_advanced(
    SuperLU_INCLUDE_DIR
    SuperLU_LIBRARY
)
