/**
 * @Author: jose
 * @Date:   2019-04-05T15:17:19+02:00
 * @Last modified by:   jose
 * @Last modified time: 2019-04-05T20:03:25+02:00
 */



#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <string>
#include <vector>
#include <valarray>
//#include <io.h>

// time management:
#include <chrono>
#include <ctime>

#include "prepro/elements/element_master.h"
#include "prepro/particles/particle.h"
#include "prepro/loader/loadGeometry.h"

using namespace std;

#define PRINT_TIME true
#define PRINT_ERRORS true
#define CHECK_AGGs_GRID false


int main(int argc, char **argv) {
    auto start = std::chrono::system_clock::now();
    auto now = start;
    auto start_part = start;
    std::chrono::duration<double> elapsed_seconds;
    if (PRINT_TIME){
      std::time_t time_now = std::chrono::system_clock::to_time_t(start);
      std::cout << "##################################################" << '\n';
      std::cout << "#### start of calculation on: " << std::ctime(&time_now);
      std::cout << "##################################################" << '\n';
    }

    int a = 0;
    pokus(a);

    if (PRINT_TIME){
      now = std::chrono::system_clock::now();
      elapsed_seconds = now-start;
      int hours = elapsed_seconds.count() / 3600;
      int minutes = elapsed_seconds.count() / 60 - hours * 60;
      int seconds = elapsed_seconds.count() - hours * 3600 - minutes * 60;
      std::time_t time_now = std::chrono::system_clock::to_time_t(now);
      std::cout << "##################################################" << '\n';
      std::cout << "#### end of calculation on: " << std::ctime(&time_now);
      std::cout << "#### total duration: " << hours << " : " << minutes << " : " << seconds << "\n";
      std::cout << "##################################################" << '\n';
    }

    return 0;
}
