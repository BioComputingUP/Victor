/*
 * Lin2dSolver.h
 *
 *  Created on: Jan 29, 2013
 *      Author: Arvind A de Menezes Pereira
 *
 * 	A simple solver for solving a system of 2D linear equations
 */

#ifndef LIN2DSOLVER_H_
#define LIN2DSOLVER_H_

#include <exception>
#include <cmath>

using namespace std;

typedef const double CDBL;

/** Invalid argument exception class. Specializing the general exception class because
 *  the std::exception specializations such as logic_error, runtime_error etc are not
 *  always defined on all platforms.
 */
class CramerUnSolvable : public exception {
public:
	explicit CramerUnSolvable( const string _what_arg, const double val ) : what_arg( _what_arg ), value( val )
	{}
	const string what_arg;
	const double value;
	virtual ~CramerUnSolvable() throw() {}
};


typedef struct Lin2dSolution {
	double x;
	double y;
	bool isValid;

	Lin2dSolution() : x(0.0), y(0.0), isValid(false) {}

	bool operator==(Lin2dSolution &ls) const {
		return (x==ls.x && y==ls.y && isValid==ls.isValid);
	}

} Lin2dSolution;

class Lin2dSolver
{
public:
	Lin2dSolver() {}

	/** Solves the 2d linear system of equations using Cramer's law.
	 *
	 * Reference: 	http://en.wikipedia.org/wiki/Cramer's_rule
	 *
	 * ax+by = e and cx+dy=f which in matrix format is: [ a b; c d ] * [ x ; y ] = [ e; f ]
	 *
	 * @param a
	 * @param b
	 * @param e
	 * @param c
	 * @param d
	 * @param f
	 * @return {Lin2dSolution} : structure containing the solution of the equation. Fields are: x,y,isValid.
	 * 			x, y contain the solution of x, y if found. isValid indicates whether a valid solution was found.
	 */
	static Lin2dSolution  Solve( CDBL a, CDBL b, CDBL e, CDBL c, CDBL d, CDBL f ) {
		Lin2dSolution l2dsol;

		const double l2dEPS = 1e-20;

		double adMinusBc = a*d-b*c;
		if( fabs( adMinusBc ) <= l2dEPS ) {
			//throw CramerUnSolvable("Error. Cannot solve using Cramer's law.",(adMinusBc));
			l2dsol.x = 0; l2dsol.y = 0; l2dsol.isValid = false;
		}
		else {
			l2dsol.x = (e*d-b*f)/adMinusBc;
			l2dsol.y = (a*f-e*c)/adMinusBc;
			l2dsol.isValid = true;
		}
		return l2dsol;
	}

	virtual ~Lin2dSolver() {}
};



#endif /* LIN2DSOLVER_H_ */
