/*
 * particle.h
 *
 *  Created on: 10-11-2012
 *      Author: Robson
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <model/primitives.hpp>

namespace dmga{

namespace model{

template<typename KeySpec>
class Particle : public Point3D{
public:
	typedef KeySpec KeyType;

	const static int RAW_SIZE = 3;

	KeyType key;

	Particle(const KeyType& key, double x, double y, double z) : Point3D(x, y, z), key(key){
	}

	Particle(const KeyType& key, double* raw) : Point3D(raw), key(key){
	}
	/**
	 * for compatibility with poly kind of diagrams
	 */
	inline double getR(){
		return 1.0;
	}

	const KeyType& getKey() const{
		return key;
	}
};

template<typename KeySpec>
class WeightedParticle : public Point4D{
public:
	typedef KeySpec KeyType;

	const static int RAW_SIZE = 4;

	KeyType key;

	WeightedParticle(const KeyType& key, double x, double y, double z, double r) : Point4D(x, y, z, r), key(key){
	}

	WeightedParticle(const KeyType& key, double* raw) : Point4D(raw), key(key){
	}

	const KeyType& getKey() const{
		return key;
	}

	double getR() const{
		return data[3];
	}
};

}//namespace model

}//namespace dmga


#endif /* PARTICLE_H_ */
