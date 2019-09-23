#ifndef GLOBAL_H
 #define GLOBAL_H

// disable uncrustify
// *INDENT-OFF*
#ifdef __has_include                           // Check if __has_include is present
 #if __has_include(<filesystem>)          // Check for a standard library
  #include <filesystem>
namespace fs = std :: filesystem;
 #elif __has_include(<experimental/filesystem>)// Check for an experimental version
  #include <experimental/filesystem>
namespace fs = std :: experimental :: filesystem;
 #elif __has_include(<boost/filesystem.hpp>)// Try with an external library
  #include <boost/filesystem.hpp>
namespace fs = boost :: filesystem;
 #else                                      // Not found at all
  #error "Missing <filesystem>"
 #endif
#endif
// continue uncrustify
//*INDENT-ON*

namespace GlobPaths {
extern fs :: path INPUT;
extern fs :: path INPUTFILENAME;
extern fs :: path BASEDIR;
extern fs :: path RESULTDIR;
}


/*
 * struct GlobPaths {
 *  static const fs::path INPUTFILE;
 *  static const fs::path BASEDIR;
 *  static const fs::path RESULTDIR;
 * };*/


#endif
