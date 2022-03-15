#ifndef _GLOBALS_H
#define _GLOBALS_H

 #include <iostream>
 #include <fstream>
 #include <stdio.h>
 #include <string>
 #include <cmath>
 #include <sstream>

// time management:
 #include <chrono>
 #include <ctime>
 #include <iomanip>

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

#define PRINT_TIME true
#define PRINT_DEBUG_TIME false
# define M_PI  3.14159265358979323846  /* pi */

namespace GlobPaths {
//extern fs :: path INPUT;
//extern fs :: path INPUTFILENAME;
extern fs :: path BASEDIR;
//extern fs :: path RESULTDIR;
}
/*
 * struct GlobPaths {
 *  static const fs::path INPUTFILE;
 *  static const fs::path BASEDIR;
 *  static const fs::path RESULTDIR;
 * };*/

//master model accessible from anywhere;
class Model;
extern Model* masterModel;

std :: string convertTimeToString(std :: chrono :: duration< double >);

#endif /* _GLOBALS_H */
