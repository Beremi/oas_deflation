#ifndef _MISC_H
#define _MISC_H
// file with help functions with logical operations like "isInVector" and transformations like PointToVector etc
#include <string>
#include <vector>
#include <algorithm>


template< typename Typo >
bool isInVect(const Typo &val, const std :: vector< Typo > &vect) {
    return std :: find(vect.begin(), vect.end(), val) != vect.end();
}

// above method does not work for strings
bool isStringInVect(const std :: string &val, const std :: vector< std :: string > &vect) {
    return std :: find(vect.begin(), vect.end(), val) != vect.end();
}





#endif /* _MISC_H */
