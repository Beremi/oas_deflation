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

/*
 * Case Sensitive Implementation of endsWith()
 * It checks if the string 'mainStr' ends with given string 'toMatch'
 */
bool endsWith(const std::string &mainStr, const std::string &toMatch)
{
    if(mainStr.size() >= toMatch.size() &&
            mainStr.compare(mainStr.size() - toMatch.size(), toMatch.size(), toMatch) == 0)
            return true;
        else
            return false;
}



#endif
