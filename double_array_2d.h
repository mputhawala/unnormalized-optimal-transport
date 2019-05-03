#ifndef _DOUBLE_ARRAY_2D_H_
#define _DOUBLE_ARRAY_2D_H_

#include <iostream>
#include <iomanip>
#include <cmath>
#include <assert.h>
// Toggle off asserts if not in debug mode
#ifndef  _DEBUG
#define _NDEBUG
#endif

namespace Math270A
{
  //
  //####################################################################
  //                    Math270A_double_array_2d.h
  //####################################################################
  /**
  Provides a very "light weight" two dimensional array structure
  with initialization capabilities and optional bounds checking.

  <pre>
  The beginning index is 0               : (C convention)
  Data for the array is stored by ROWS   : (C convention)
  Access using (*,*), e.g. A(i,j) for (i,j)th element.: (NOT C convention)
  </pre>
  Values are stored as:

  (0,0)	(0,1)	(0,2)	...
  (1,0)	(1,1)	(1,2)	...
  (2,0)	(2,1)	(2,2)	...
  .		  .		  .		.
  .		  .		  .		 .
  .		  .		  .		  .

  The copy constructor creates a duplicate instance.<p>

  Bounds checking is toggled on by specifying the compiler pre-processor
  define _DEBUG
  <p>
  Created for use in Math 270A<p>

  @author Chris Anderson (C) UCLA
  @version  Jan. 26, 2015
  */
  //#####################################################################
  //
  //#####################################################################
  //

  class DoubleArray2D
  {

  public:
    //
    //###################################################################
    //                 Constructors/Initialization
    //###################################################################
    //
    DoubleArray2D();
    DoubleArray2D(size_t m, size_t n);

    DoubleArray2D(const DoubleArray2D& D);

    virtual ~DoubleArray2D();

    void Initialize();
    void Initialize(const DoubleArray2D& D);
    void Initialize(size_t m, size_t n);


    void operator=(const DoubleArray2D& D);
    void operator+=(const DoubleArray2D& G);
    void operator-=(const DoubleArray2D& G);
    void operator*=(double alpha);
    void operator/=(double alpha);

    void SetToValue(double d);

    double NormInf() const;
    //
    //###################################################################
    //      Element Access with bounds checking toggled
    //      using _DEBUG
    //###################################################################
    //
    double&  operator()(size_t i1, size_t i2);

    const double&  operator()(size_t i1, size_t i2) const;
    //
    //###################################################################
    //                Array Structure Access Functions
    //###################################################################
    //
    /// Returns the number of indices in the 1st coordinate direction 

    size_t GetIndex1Size()  const;

    /// Returns the number of indices in the 2nd coordinate direction 

    size_t GetIndex2Size()  const;

    /// Returns the total number of data values 

    size_t GetDataSize()    const;

    /// Returns a pointer to the double array containing the data values 

    double* const GetDataPointer() const;

    //  Input/Output
    //
    //  Prints out values as as if they were in the first Cartesian 
    //  quadrant --- not in matrix indexing. 
    //
    //
    friend std::ostream&  operator <<(std::ostream& out_stream, const DoubleArray2D& A);
    //
    //###################################################################
    //                      Class Data Members
    //###################################################################
    //
  protected:

    double*      data_ptr;     // data pointer
    size_t      index_1_size;     // coordinate 1 size
    size_t      index_2_size;     // coordinate 2 size
    //
    //###################################################################
    //                      Bounds Checking
    //###################################################################
    //
    bool BoundsCheck(size_t i, size_t begin, size_t end, int coordinate) const;
  };


}       /* Namespace Math270A */
#endif  /* _DOUBLE_ARRAY_2D_H_ */
