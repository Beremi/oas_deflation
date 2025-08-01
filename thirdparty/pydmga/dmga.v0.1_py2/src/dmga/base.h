/*
 * base.h
 *
 *  Created on: 24-04-2013
 *      Author: robson
 */

#ifndef BASE_H_
#define BASE_H_

#ifdef DEBUG
#include <iostream>
#endif
#include <string>

#ifndef DEBUG
#define DEB(s)
#else
#define DEB(s) std::cerr << "level(0): " << s << std::endl << std::flush
#endif

#if DEBUG < 1
#define DEB1(s)
#else
#define DEB1(s) std::cerr << "level(1): " << s << std::endl << std::flush
#endif

#if DEBUG < 2
#define DEB2(s)
#else
#define DEB2(s) std::cerr << "level(2): " << s << std::endl << std::flush
#endif

#if DEBUG < 3
#define DEB3(s)
#else
#define DEB3(s) std::cerr << "level(3): " << s << std::endl << std::flush
#endif

namespace dmga{

class DMGAObject{
public:
	virtual std::string show(){
		return "DMGAObject";
	}
	virtual ~DMGAObject(){
		DEB3("DMGAObject::__destruct__(): " << show());
	}
};

}//namespace dmga


#endif /* BASE_H_ */
