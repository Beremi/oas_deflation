/*
 * trajectory.h
 *
 *  Created on: 16-11-2012
 *      Author: Robson
 */

#ifndef TRAJECTORY_H_
#define TRAJECTORY_H_

#include <model/model.hpp>

namespace dmga{

namespace trajectory{

class Trajectory{
public:
	double current_time;

	virtual bool advance() = 0;

	virtual ~Trajectory(){
	}

	virtual double* getRaw(int i) = 0;

	virtual dmga::model::Point3D getPos(int i) = 0;
};

} //namespace trajectory

} //namespace dmga


#endif /* TRAJECTORY_H_ */
