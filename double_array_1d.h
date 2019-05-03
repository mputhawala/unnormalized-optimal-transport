#ifndef _DOUBLE_ARRAY_1D_H_
#define _DOUBLE_ARRAY_1D_H_

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
  //                    Math270A_double_array_1d.h
  //####################################################################
  /**
  Provides a very "light weight" one dimensional array structure
  with initialization capabilities and optional bounds checking.

  <pre>
  The beginning index is 0               : (C convention)
  Access using (*), e.g. A(i) for (i)th element.
  </pre>
  Values are stored as:

  (1), (2), ...

  The copy constructor creates a duplicate instance.<p>

  Bounds checking is toggled on by specifying the compiler pre-processor
  define _DEBUG
  <p>
  Created for use in Math 270A<p>

  @author Chris Anderson (C) UCLA
  @and Michael Puthawala
  @version  3/6/2019
  */
  //#####################################################################
  //
  //#####################################################################
  //

  class DoubleArray1D
  {

  public:
    //
    //###################################################################
    //                 Constructors/Initialization
    //###################################################################
    //
    DoubleArray1D();
    DoubleArray1D(size_t m);
    DoubleArray1D(const DoubleArray1D& D);

    virtual ~DoubleArray1D();

    void Initialize();
    void Initialize(const DoubleArray1D& D);
    void Initialize(size_t m);

    void operator=(const DoubleArray1D& D);
    void operator+=(const DoubleArray1D& G);
    void operator-=(const DoubleArray1D& G);

    void operator*=(double alpha);
    void operator/=(double alpha);

    void SetToValue(double d);

    double NormInf();
    //
    //###################################################################
    //      Element Access with bounds checking toggled
    //      using _DEBUG
    //###################################################################
    //
#ifdef _DEBUG 
    double&  operator()(size_t i1);

    const double&  operator()(size_t i1) const;
#else
    double&  operator()(size_t i1);

    const double&  operator()(size_t i1) const;
#endif

    //
    //###################################################################
    //                Array Structure Access Functions
    //###################################################################
    //
    /// Returns the number of indices in the 1st coordinate direction 

    size_t GetIndex1Size() const;

    /// Returns the number of indices in the 2nd coordinate direction 

    /// Returns the total number of data values 

    size_t GetDataSize() const;

    /// Returns a pointer to the double array containing the data values 

    double* const GetDataPointer() const;

    //  Input/Output
    //
    //  Prints out values as as if they were in the first Cartesian 
    //  quadrant --- not in matrix indexing. 
    //
    //
    friend std::ostream&  operator <<(std::ostream& out_stream, const DoubleArray1D& A);
    //
    //###################################################################
    //                      Class Data Members
    //###################################################################
    //
  protected:

    double*      data_ptr;     // data pointer
    size_t      index_1_size;     // coordinate 1 size
    //
    //###################################################################
    //                      Bounds Checking
    //###################################################################
    //
    bool BoundsCheck(size_t i, size_t begin, size_t end, int coordinate) const;
  };


}       /* Namespace Math270A */
#endif  /* _DoubleArray1D_    */
