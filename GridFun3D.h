//######################################################################
//   Class GridFun3D
//######################################################################
/*
Class GridFun3D: A class which represents a scalar valued function on a
3 dimentional rectangular surface.

The function values itself are to be specified in the values object,
whereas everything else is not used in this class, but is useful for
RelaxOp2D



*/
// Math 270A
// Mar. 1, 2015
//######################################################################
//
//
#ifndef _GridFun3D_h_
#define _GridFun3D_h_

#include "Math270A_DoubleArray3D.h"

#include <assert.h>
#include <iostream>

class GridFun3D
{
public:
	GridFun3D();
	GridFun3D(const GridFun3D& G);
	GridFun3D(long xPanel, double xMin, double xMax,
		long yPanel, double yMin, double yMax,
		long zPanel, double zMin, double zMax);
	~GridFun3D();
	void initialize();
	void initialize(const GridFun3D& G);
	void initialize(long xPanel, double xMin, double xMax, 
		long yPanel, double yMin, double yMax,
		long zyPanel, double zMin, double zMax);
	void operator =(const GridFun3D& G); // Equates 2 function
	void operator +=(const GridFun3D& G); // Adds two function, equiv to (f+g)(x)
	void operator -=(const GridFun3D& G); // (f-g)(x)
	void operator *=(double alpha); // Scalar multiplication
	void operator /=(double alpha); // Scalar division
	void setToValue(double d);
	double normInf();
	friend ostream& operator <<(ostream& outStream, const GridFun3D& V);
	Math270A::DoubleArray3D values;
	double hx;		// Spatial resolution
	double xMin;	// Where the x values start
	double xMax;	// and end
	long xPanel;	// Number of panels total. Number of grid points = panels + 1
	double hy;
	double yMin;
	double yMax;
	long yPanel;
	double hz;
	double zMin;
	double zMax;
	long zPanel;
};

#endif