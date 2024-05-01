#include "version.h"

std :: string version_info(bool diff) {
    std :: string s = "This code has been built from version " + GIT_HASH + " : " + TIME_STRING + "\n";
    s += "OS name: " + OS_NAME + "\n";
    s += "OS sub-type: " + OS_RELEASE + "\n";
    s += "OS build ID: " + OS_VERSION + "\n";
    s += "OS platform: " + OS_PLATFORM + "\n";
    
    s += "Build type: ";
    #ifdef _MSC_VER
        #ifdef _DEBUG
            // Code specific to Debug build
            s += "DEBUG";
            s += "\n";
        #elif NDEBUG
            // Code specific to Release build
            s += "RELEASE";
            s += "\n";
        #endif
    #else
        s += BUILD_TYPE;
        s += "\n";
    #endif
    
    s += "Build contains: | ";
    #ifdef _OPENMP
        s += "OPENMP | ";
    #endif
        
    #ifdef VTK_FOUND
        s += "VTK (" + VTK_VERSION + ") | ";
    #endif

    #ifdef EIGEN_USE_LAPACKE
        s += "LAPACK | ";
    #endif
        
    #ifdef EIGEN_USE_BLAS
        s += "BLAS | ";
    #endif
    
    if (diff) s += "\nGit diff: " + GIT_DIFF + "\n";
    
    return s;
}
