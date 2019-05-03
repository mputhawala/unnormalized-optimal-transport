#include "double_array_1d.h"

namespace Math270A {
  DoubleArray1D::DoubleArray1D() {
    data_ptr = 0;
    index_1_size = 0;
  }

  DoubleArray1D::DoubleArray1D(size_t m) {
    data_ptr = 0;
    index_1_size = 0;
    Initialize(m);
  };

  DoubleArray1D::DoubleArray1D(const DoubleArray1D& D) {
    data_ptr = 0;
    index_1_size = 0;
    Initialize(D);
  };

  DoubleArray1D::~DoubleArray1D() {
    if (data_ptr != 0) delete[] data_ptr;
  }

  void DoubleArray1D::Initialize() {
    if (data_ptr != 0) delete[] data_ptr;
    data_ptr = 0;
    index_1_size = 0;
  };


  void DoubleArray1D::Initialize(const DoubleArray1D& D) {
    // Initialization of null object

    if (D.data_ptr == 0)
    {
      Initialize();
      return;
    }

    if ((index_1_size != D.index_1_size))
    {
      if (data_ptr != 0) delete[] data_ptr;
      index_1_size = D.index_1_size;
      data_ptr = new double[index_1_size];
    }
    size_t i;
    for (i = 0; i < index_1_size; i++)
    {
      data_ptr[i] = D.data_ptr[i];
    }
  };


  void DoubleArray1D::Initialize(size_t m) {
    if ((index_1_size != m))
    {
      if (data_ptr != 0) delete[] data_ptr;
      index_1_size = m;
      data_ptr = new double[index_1_size];
    }
    size_t i;
    for (i = 0; i < index_1_size; i++)
    {
      data_ptr[i] = 0.0;
    }
  };

  void DoubleArray1D::operator=(const DoubleArray1D& D) {
    Initialize(D);
  }
  size_t i2;
  void DoubleArray1D::operator+=(const DoubleArray1D& G)
  {
    assert(index_1_size == G.index_1_size);
    for (size_t i1 = 0; i1 < index_1_size; i1++)
    {
      *(data_ptr + i1) = *(data_ptr + i1) + *(G.data_ptr + i1);
    }
  }

  void DoubleArray1D::operator-=(const DoubleArray1D& G) {
    assert(index_1_size == G.index_1_size);
    for (size_t i1 = 0; i1 < index_1_size; i1++)
      *(data_ptr + i1) = *(data_ptr + i1) - *(G.data_ptr + i1);
  }

  void DoubleArray1D::operator*=(double alpha) {
    for (size_t i1 = 0; i1 < index_1_size; i1++)
      *(data_ptr + i1) = *(data_ptr + i1)*alpha;
  }

  void DoubleArray1D::operator/=(double alpha) {
    assert(alpha != 0);
    for (size_t i1 = 0; i1 < index_1_size; i1++)
      *(data_ptr + i1) = *(data_ptr + i1) / alpha;
  }

  void DoubleArray1D::SetToValue(double d) {
    for (size_t i1 = 0; i1 < index_1_size; i1++)
      *(data_ptr + i1) = d;
  }

  double DoubleArray1D::NormInf() {
    double norm = 0;
    for (size_t i1 = 0; i1 < index_1_size; i1++)
      if (abs(*(data_ptr + i1)) > norm)
        norm = abs(*(data_ptr + i1));
    return norm;
  }
  //
  //###################################################################
  //      Element Access with bounds checking toggled
  //      using _DEBUG
  //###################################################################
  //
#ifdef _DEBUG 
  double&  DoubleArray1D::operator()(size_t i1) {
    assert(BoundsCheck(i1, 0, index_1_size - 1, 1));
    return *(data_ptr + i1);
  };

  const double&  DoubleArray1D::operator()(size_t i1) const {
    assert(BoundsCheck(i1, 0, index_1_size - 1, 1));
    return *(data_ptr + i1);
  };
#else
  double&  DoubleArray1D::operator()(size_t i1) {
    return *(data_ptr + i1);
  };

  const double&  DoubleArray1D::operator()(size_t i1) const {
    return *(data_ptr + i1);
  };
#endif

  //
  //###################################################################
  //                Array Structure Access Functions
  //###################################################################
  //
  /// Returns the number of indices in the 1st coordinate direction 

  size_t DoubleArray1D::GetIndex1Size()  const { return index_1_size; }

  /// Returns the number of indices in the 2nd coordinate direction 

  /// Returns the total number of data values 

  size_t DoubleArray1D::GetDataSize()    const { return index_1_size; }

  /// Returns a pointer to the double array containing the data values 

  double* const DoubleArray1D::GetDataPointer() const { return data_ptr; }

  //  Input/Output
  //
  //  Prints out values as as if they were in the first Cartesian 
  //  quadrant --- not in matrix indexing. 
  //
  //
  std::ostream&  operator <<(std::ostream& out_stream, const DoubleArray1D& A)
  {
    for (size_t i = 0; i < A.index_1_size - 1; i++)
    {
      // out_stream << setw(12) << A(i) << " ";
      out_stream << A(i) << '\t';
    }
    out_stream << A(A.index_1_size - 1);
    return out_stream;
  }

#ifdef _DEBUG
  bool DoubleArray1D::BoundsCheck(size_t i, size_t begin, size_t end, int coordinate) const {
    if ((i < begin) || (i > end))
    {
      std::cerr << "Array index " << coordinate << " out of bounds " << std::endl;
      std::cerr << "Offending index value : " << i << " Acceptable Range [" << begin << "," << end << "]" << std::endl;
      return false;
    }
    return true;
  }
#else
  bool DoubleArray1D::BoundsCheck(size_t, size_t, size_t, int) const { return true; }
#endif

}