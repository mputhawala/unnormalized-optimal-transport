#ifndef _DOUBLE_ARRAY_3D_H_
#define _DOUBLE_ARRAY_3D_H_

#include <iostream>
#include <cassert>

// Toggle off asserts if not in debug mode
#ifndef  _DEBUG
#define _NDEBUG
#endif

namespace Math270A
{
  //
  //####################################################################
  //                    Math270A_double_array_3d.h
  //####################################################################
  /**
  Provides a very "light weight" two dimensional array structure
  with initialization capabilities and optional bounds checking.

  <pre>
  The beginning index is 0               : (C convention)
  Data for the array is stored by ROWS   : (C convention)
  Access using (*,*,*), e.g. A(i,j,k) for (i,j,k)th element.: (NOT C convention)
  </pre>
  Values are stored as:

  (0,0,0)	(0,0,1)	(0,0,2)	... 	(1,0,0)	(1,0,1)	(1,0,2)	... (2,0,0)	(2,0,1)	(2,0,2)	...
  (0,1,0)	(0,1,1)	(0,1,2)	... 	(1,1,0)	(1,1,1)	(1,1,2)	... (2,1,0)	(2,1,1)	(2,1,2)	...
  (0,2,0)	(0,2,1)	(0,2,2)	... 	(1,2,0)	(1,2,1)	(1,2,2)	... (2,2,0)	(2,2,1)	(2,2,2)	...
  .		  .		  .		.   	.		  .		  .		.   .		  .		  .		.
  .		  .		  .		 .  	.		  .		  .		 .  .		  .		  .		 .
  .		  .		  .		  . 	.		  .		  .		  . .		  .		  .		  .

  The copy constructor creates a duplicate instance.<p>

  Bounds checking is toggled on by specifying the compiler pre-processor
  define _DEBUG
  <p>

  @author Chris Anderson (C) UCLA
  @ and Michael Puthawala
  @version  Dec. 24 2018
  */
  //#####################################################################
  //
  //#####################################################################
  //

  class DoubleArray3D {
  private:
    void Initialize();


    void Initialize(const DoubleArray3D& D);


    void Initialize(size_t m, size_t n, size_t o);

  public:
    //
    //###################################################################
    //                 Constructors/Initialization
    //###################################################################
    //
    DoubleArray3D();

    DoubleArray3D(size_t m, size_t n, size_t o);

    DoubleArray3D(const DoubleArray3D& D);

    virtual ~DoubleArray3D();

    void Resize(size_t m, size_t n, size_t o);

    DoubleArray3D& operator=(const DoubleArray3D& D);

    void operator+=(const DoubleArray3D& rhs);

    void operator-=(const DoubleArray3D& rhs);

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
    double& operator()(size_t i1, size_t i2, size_t i3);

    const double& operator()(size_t i1, size_t i2, size_t i3) const;
#else
    double& operator()(size_t i1, size_t i2, size_t i3);

    const double& operator()(size_t i1, size_t i2, size_t i3) const;
#endif

    bool operator==(const DoubleArray3D& rhs);

    //
    //###################################################################
    //                Array Structure Access Functions
    //###################################################################
    //
    /// Returns the number of indices in the 1st coordinate direction 

    inline size_t GetIndex1Size()  const { return index_1_size; }

    /// Returns the number of indices in the 2nd coordinate direction 

    inline size_t GetIndex2Size()  const { return index_2_size; }

    /// Returns the number of indices in the 3rd coordinate direction 

    inline size_t GetIndex3Size()  const { return index_3_size; }

    /// Returns the total number of data values 

    inline size_t GetDataSize()    const { return index_1_size*index_2_size*index_3_size; }

    /// Returns a pointer to the double array containing the data values 

    inline double* const GetDataPointer() const { return data_ptr; }

    //  Input/Output
    //
    //  Prints out values as as if they were in the first Cartesian 
    //  quadrant --- not in matrix indexing. 
    //
    //
    friend std::ostream&  operator <<(std::ostream& outStream, const DoubleArray3D& A);

    //
    //###################################################################
    //                      Class Data Members
    //###################################################################
    //
  protected:

    double*     data_ptr;     // data pointer
    size_t      index_1_size;     // coordinate 1 size
    size_t      index_2_size;     // coordinate 2 size
    size_t      index_3_size;     // coordinate 3 size
    //
    //###################################################################
    //                      Bounds Checking
    //###################################################################
    //
#ifdef _DEBUG
    bool BoundsCheck(size_t i, size_t begin, size_t end, int coordinate) const;
#else
    inline bool BoundsCheck(size_t, size_t, size_t, int) const;
#endif
  };
}       /* Namespace Math270A */
#endif  /* _DOUBLE_ARRAY_3D_H_*/
