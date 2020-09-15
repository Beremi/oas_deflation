#include "version.h"

std :: string version_info() {
    std :: string s = "This code has been built from version " + GIT_HASH + " : " + TIME_STRING + "\n";
    s += "OS name: " + OS_NAME + "\n";
    s += "OS sub-type: " + OS_RELEASE + "\n";
    s += "OS build ID: " + OS_VERSION + "\n";
    s += "OS platform: " + OS_PLATFORM + "\n";
    s += "Build type: " + BUILD_TYPE;
    return s;
}
