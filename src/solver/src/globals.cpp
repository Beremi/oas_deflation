#include "globals.h"

using namespace std;

string convertTimeToString(std :: chrono :: duration< double >time_interval) {
    if (time_interval.count() <= 1e-4){
        return std :: to_string(time_interval.count());
    }
    int hours = time_interval.count() / 3600;
    int minutes = time_interval.count() / 60 - hours * 60;
    int seconds = time_interval.count() - hours * 3600 - minutes * 60;
    int miliseconds = ( time_interval.count() - hours * 3600 - minutes * 60 - seconds ) * 1000;
    stringstream ss;
    ss << std :: setw(2) << std :: setfill('0') << hours << ":"
       << std :: setw(2) << std :: setfill('0') << minutes << ":"
       << std :: setw(2) << std :: setfill('0') << seconds << "."
       << std :: setw(3) << std :: setfill('0') << miliseconds;
    return ss.str();
}

