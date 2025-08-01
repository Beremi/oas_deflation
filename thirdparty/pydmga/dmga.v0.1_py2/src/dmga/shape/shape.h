/*
 * shape.h
 *
 *  Created on: 12-11-2012
 *      Author: Robson
 */

#ifndef SHAPE_H_
#define SHAPE_H_

#include <base.h>
#include <diagram/diagram.h>
#include <vector>
#include <model/model.hpp>
#include <utils/utils.hpp>
#include <algorithm>
#include <algebra/algebra.hpp>

namespace dmga{

namespace shape{


class GenericShape : public DMGAObject{
public:
	std::string show(){ return "dmga::shape::GenericShape"; }

	~GenericShape(){ 			DEB3("GenericShape::__destruct__():");
		//destructor...
	}
};

class GenericKineticShape : public DMGAObject{
public:
	std::string show(){ return "dmga::shape::GenericKineticShape"; }

	~GenericKineticShape(){ 			DEB3("GenericKineticShape::__destruct__():");
		//destructor...
	}
};

} //namespace shape

} //namespace dmga


#endif /* SHAPE_H_ */
