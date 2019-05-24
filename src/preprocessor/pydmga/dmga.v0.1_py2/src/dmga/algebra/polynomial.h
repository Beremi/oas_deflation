/*
 * polynomial.h
 *
 *  Created on: 14-01-2013
 *      Author: Robson
 */

#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#include <iostream>

namespace dmga{

namespace algebra{

namespace polynomial{

/**
 * represents standard polynomial
 * by its coefficients
 */
class Polynomial{
public:
	/**
	 * creates a polynomial of (theoretical) degree n.
	 * If coeffs is null then all coefficients are set to 0.0
	 * otherwise we copy (n+1) coefficients from data (data[j] is a_j i.e. the coefficient by x^j)
	 */
	Polynomial(int n, double* data): n(n){
		coeffs = new double[n + 1];
		for (int i = 0; i <= n; i++){
			coeffs[i] = data[i];
		}
	}
	/**
	 * creates a polynomial of (theoretical) degree n.
	 * that represents constant polynomial with a_0 = a0.
	 */
	Polynomial(int n, double a0 = 0.0): n(n){
		coeffs = new double[n + 1];
		coeffs[0] = a0;
		for (int i = 1; i <= n; i++){
			coeffs[i] = 0.0;
		}
	}
	/** output polynomial */
	friend std::ostream& operator<<(std::ostream& out, const Polynomial& item){
		out << item.coeffs[0];
		for (int i = 1; i <= item.n; i++){
			out << " + " << item.coeffs[i] << " x^" << i;
		}
		return out;
	}

	/** degree */
	int n;
	/** coefficients coeffs[i] = a_i where a_i stands by x^i */
	double* coeffs;
private:
};

/**
 * this class allows for holding the precomputed
 * base for a family of Newton Polynomials for multiple functions
 * with the same spaced points a_1, .., a_N
 *
 * this base allows for fast retrival of normal Polynomial from
 * given NewtonPolynomial
 *
 * we hold the information about f(x) - the polinomial
 * and the f(x+1) - the polinomial computed at the x+1
 * this is because we will use it in root finding
 * for the interval [0,1] in our Kinetic Voronoi computation
 */
class NewtonBase{
public:
	void init(int N, double* pts){
		conv = new double[(n+1) * (n+1)];
		points = new double[n+1];
		for (int i = 0; i <= n; i++){
			points[i] = pts[i];
		}
		conv[0] = 1.0;
		for (int i = 1; i <= n; i++){
			conv[index(i, 0)] = conv[index(i-1, 0)] * (-points[i-1]);
			conv[xindex(i, 0)] = conv[xindex(i-1, 0)] * (1.0 - points[i-1]);
			conv[index(i, i)] = 1.0;
		}
		for (int k = 2; k <= n; k++){
			for (int j = 1; j < k; j++){
				conv[index(k, j)] = conv[index(k-1, j-1)] - conv[index(k-1, j)] * points[k-1];
				conv[xindex(k, j)] = conv[index(k-1, j-1)] + conv[xindex(k-1, j)] * (1.0 - points[k-1]);
			}
		}
	}
	NewtonBase(): n(0), conv(0), points(0){
	}
	/**
	 * @param int N degree of the base
	 * @param double* pts the points of the interpolation
	 */
	NewtonBase(int N, double* pts): n(N){
		init(N, pts);
	}
	/** index in the conv where j-th coefficient of k-th base polynomial f(x) is stored */
	inline int index(int k, int j){
		return k * (n+1) + j;
	}
	/** index in the conv where j-th coefficient of k-th base polynomial f(x+1) is stored */
	inline int xindex(int k, int j){
		return j * (n+1) + k;
	}

	friend std::ostream& operator<<(std::ostream& out, const NewtonBase& item){
		for (int k = 0; k <= item.n; k++){
			out << item.conv[k * (item.n+1)];
			for (int i = 1; i <= item.n; i++){
				out << "\t" << item.conv[k * (item.n+1) + i];
			} out << "\n";
		}
		return out;
	}

	/** degree of polynomial */
	int n;
	/** conversion table from Neton to Normal */
	double* conv;
	/** points on which the interpolation polynomial is based */
	double* points;
};
/**
 * Does the newton interpolation for given data points
 */
class NewtonPolynomial{
public:
	NewtonBase& base;
	double* f;

	/**
	 * computes polynomial interpolation in points given by base
	 * and values given in values
	 */
	NewtonPolynomial(NewtonBase& base, double* values): base(base){
		double* divdiff = new double[base.n + 1];
		f = new double[base.n + 1];
		for (int i = 0; i <= base.n; i++){
			divdiff[i] = values[i];
		}
		for (int i = 0; i < base.n; i++){
			f[i] = divdiff[0];
			for (int j = 0; j < base.n - i; j++){
				divdiff[j] = (divdiff[j+1] - divdiff[j]) / (base.points[j + 1 + i] - base.points[j]);
			}
		}
		f[base.n] = divdiff[0];
		delete[] divdiff;
	}

	operator Polynomial(){
		return normal_polynomial();
	}

	Polynomial normal_polynomial(){
		Polynomial p(base.n);
		for (int j = 0; j <= base.n; j++){ // for each coefficient, we want to compute...
			for (int k = j; k <= base.n; k++){ //for each basic polynomial, which can contain a coefficient of order x^j
				p.coeffs[j] += f[k] * base.conv[base.index(k, j)];
			}
		}
		return p;
	}

	Polynomial shifted_polynomial(){
		Polynomial p(base.n);
		for (int j = 0; j <= base.n; j++){ // for each coefficient, we want to compute...
			for (int k = j; k <= base.n; k++){ //for each basic polynomial, which can contain a coefficient of order x^j
				p.coeffs[j] += f[k] * base.conv[base.xindex(k, j)];
			}
		}
		return p;
	}

	friend std::ostream& operator<<(std::ostream& out, const NewtonPolynomial& item){
		out << item.f[0];
		for (int i = 1; i <= item.base.n; i++){
			out << " " << item.f[i];
		}
		return out;
	}
private:
};

} //namespace polynomial

} //namespace algebra

} //namespace dmga


#endif /* POLYNOMIAL_H_ */
