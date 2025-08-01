/*
 * primitives.h
 *
 *  Created on: 10-11-2012
 *      Author: Robson
 */

#ifndef PRIMITIVES_H_
#define PRIMITIVES_H_

namespace dmga{

namespace model{

template<int N>
class Point{
public:
//	double x;
//	double y;
//	double z;
	double data[N];

	Point(){
	}

	Point(bool clear){
		if (clear){
			for (int i = 0; i < N; ++i){
				data[i] = 0.0;
			}
		}
	}

	Point(double x, double y, double z){
		data[0] = x;
		data[1] = y;
		data[2] = z;
		for (int i = 3; i < N; ++i){
			data[i] = 0.0;
		}
	}

	Point(double* raw){
		for (int i = 0; i < N; ++i){
			data[i] = raw[i];
		}
	}

	inline double getX() const{
		return data[0];
	}

	inline double getY() const{
		return data[1];
	}
	inline double getZ() const{
		return data[2];
	}
	inline double getI(int i) const{
		return data[i];
	}
	inline const double* getRaw() const{
		return &data[0];
	}
	inline double* getRaw(){
		return &data[0];
	}

	friend std::ostream& operator<<(std::ostream& out, const Point<N>& point){
		out << "(";
		for (int i = 0; i < N-1; i++){
			out << point.data[i] << ", ";
		}
		out << point.data[N-1] << ")";
		return out;
	}
};

class Point3D : public Point<3>{
public:
	Point3D(double x, double y, double z){
		data[0] = x;
		data[1] = y;
		data[2] = z;
	}
	Point3D(double* raw){
		data[0] = raw[0];
		data[1] = raw[1];
		data[2] = raw[2];
	}
};

class Point4D : public Point<4>{
public:
	Point4D(double x, double y, double z, double r = 1.0){
		data[0] = x;
		data[1] = y;
		data[2] = z;
		data[3] = r;
	}
	Point4D(double* raw){
		data[0] = raw[0];
		data[1] = raw[1];
		data[2] = raw[2];
		data[3] = raw[3];
	}
};

}//namespace model

}//namespace dmga


#endif /* PRIMITIVES_H_ */
