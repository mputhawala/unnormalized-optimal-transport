#include <iostream>
#include <iomanip>
#include <cmath>
#include "double_array_3d.h"

// Toggle off asserts if not in debug mode
#ifndef  _DEBUG
#define _NDEBUG
#endif

namespace Math270A
{
  void DoubleArray3D::Initialize() {
    if (data_ptr != nullptr)
      delete[] data_ptr;
    data_ptr = nullptr;
    index_1_size = 0;
    index_2_size = 0;
    index_3_size = 0;
  };


  void DoubleArray3D::Initialize(const DoubleArray3D& D) {
    // Initialization of null object

    if (D.data_ptr == nullptr) {
      Initialize();
      return;
    }

    if ((index_1_size != D.index_1_size) || (index_2_size != D.index_2_size) || (index_3_size != D.index_3_size)) {
      if (data_ptr != nullptr) {
        delete[] data_ptr;
        data_ptr = nullptr;
      }
      index_1_size = D.index_1_size;
      index_2_size = D.index_2_size;
      index_3_size = D.index_3_size;
    }
    if (data_ptr == nullptr) {
      data_ptr = new double[index_1_size*index_2_size*index_3_size];
    }
    for (size_t i = 0; i < index_1_size*index_2_size*index_3_size; ++i) {
      data_ptr[i] = D.data_ptr[i];
    }
  };


  void DoubleArray3D::Initialize(size_t m, size_t n, size_t o) {
    if ((index_1_size != m) || (index_2_size != n) || (index_3_size != o)) {
      if (data_ptr != nullptr) {
        delete[] data_ptr;
        data_ptr = nullptr;
      }
      index_1_size = m;
      index_2_size = n;
      index_3_size = o;
    }
    if (data_ptr == nullptr) {
      data_ptr = new double[index_1_size*index_2_size*index_3_size];
    }
    for (size_t i = 0; i < index_1_size*index_2_size*index_3_size; ++i) {
      data_ptr[i] = 0.0;
    }
  };

  //
  //###################################################################
  //                 Constructors/Initialization
  //###################################################################
  //
  DoubleArray3D::DoubleArray3D() : data_ptr(nullptr), index_1_size(0), index_2_size(0), index_3_size(0) {
  };

  DoubleArray3D::DoubleArray3D(size_t m, size_t n, size_t o) : data_ptr(nullptr), index_1_size(m), index_2_size(n), index_3_size(o) {
    Initialize(m, n, o);
  };

  DoubleArray3D::DoubleArray3D(const DoubleArray3D& D) : data_ptr(nullptr), index_1_size(D.index_1_size), index_2_size(D.index_2_size), index_3_size(D.index_3_size) {
    Initialize(D);
  };

  DoubleArray3D::~DoubleArray3D() {
    std::cout << "DoubleArray3D::~DoubleArray3D" << std::endl;
    if (data_ptr != nullptr) {
      delete[] data_ptr;
      data_ptr = nullptr;
    }
  }

  void DoubleArray3D::Resize(size_t m, size_t n, size_t o) {
    Initialize(m, n, o);
  }

  DoubleArray3D& DoubleArray3D::operator=(const DoubleArray3D& D) {
    Initialize(D);
    return *this;
  }

  void DoubleArray3D::operator+=(const DoubleArray3D& rhs) {
    assert(index_1_size == rhs.index_1_size);
    assert(index_2_size == rhs.index_2_size);
    assert(index_3_size == rhs.index_3_size);
    for (size_t i1 = 0; i1 < index_1_size; ++i1) {
      for (size_t i2 = 0; i2 < index_2_size; ++i2) {
        for (size_t i3 = 0; i3 < index_3_size; ++i3) {
          *(data_ptr + i1*index_2_size*index_3_size + i2*index_3_size + i3) += *(rhs.data_ptr + i1*index_2_size*index_3_size + i2*index_3_size + i3);
        }
      }
    }
  }

  void DoubleArray3D::operator-=(const DoubleArray3D& rhs) {
    assert(index_1_size == rhs.index_1_size);
    assert(index_2_size == rhs.index_2_size);
    assert(index_3_size == rhs.index_3_size);
    for (size_t i1 = 0; i1 < index_1_size; ++i1) {
      for (size_t i2 = 0; i2 < index_2_size; ++i2) {
        for (size_t i3 = 0; i3 < index_3_size; ++i3) {
          *(data_ptr + i1*index_2_size*index_3_size + i2*index_3_size + i3) -= *(rhs.data_ptr + i1*index_2_size*index_3_size + i2*index_3_size + i3);
        }
      }
    }
  }

  void DoubleArray3D::operator*=(double alpha) {
    for (size_t i1 = 0; i1 < index_1_size; ++i1) {
      for (size_t i2 = 0; i2 < index_2_size; ++i2) {
        for (size_t i3 = 0; i3 < index_3_size; ++i3) {
          *(data_ptr + i1*index_2_size*index_3_size + i2*index_3_size + i3) *= alpha;
        }
      }
    }
  }

  void DoubleArray3D::operator/=(double alpha) {
    for (size_t i1 = 0; i1 < index_1_size; ++i1) {
      for (size_t i2 = 0; i2 < index_2_size; ++i2) {
        for (size_t i3 = 0; i3 < index_3_size; ++i3) {
          *(data_ptr + i1*index_2_size*index_3_size + i2*index_3_size + i3) /= alpha;
        }
      }
    }
  }

  void DoubleArray3D::SetToValue(double d) {
    for (size_t i1 = 0; i1 < index_1_size; ++i1) {
      for (size_t i2 = 0; i2 < index_2_size; ++i2) {
        for (size_t i3 = 0; i3 < index_3_size; ++i3) {
          *(data_ptr + i1*index_2_size*index_3_size + i2*index_3_size + i3) = d;
        }
      }
    }
  }

  double DoubleArray3D::NormInf() {
    double norm = 0;
    for (size_t i1 = 0; i1 < index_1_size; ++i1) {
      for (size_t i2 = 0; i2 < index_2_size; ++i2) {
        for (size_t i3 = 0; i3 < index_3_size; ++i3) {
          if (abs(*(data_ptr + i1*index_2_size*index_3_size + i2*index_3_size + i3)) > norm)
            norm = abs(*(data_ptr + i1*index_2_size*index_3_size + i2*index_3_size + i3));
        }
      }
    }
    return norm;
  }
  //
  //###################################################################
  //      Element Access with bounds checking toggled
  //      using _DEBUG
  //###################################################################
  //
#ifdef _DEBUG 
  double& DoubleArray3D::operator()(size_t i1, size_t i2, size_t i3) {
    assert(BoundsCheck(i1, 0, index_1_size - 1, 1));
    assert(BoundsCheck(i2, 0, index_2_size - 1, 2));
    assert(BoundsCheck(i3, 0, index_3_size - 1, 2));
    return *(data_ptr + i1*index_2_size*index_3_size + i2*index_3_size + i3);
  };

  const double& DoubleArray3D::operator()(size_t i1, size_t i2, size_t i3) const {
    assert(BoundsCheck(i1, 0, index_1_size - 1, 1));
    assert(BoundsCheck(i2, 0, index_2_size - 1, 2));
    assert(BoundsCheck(i3, 0, index_3_size - 1, 2));
    return *(data_ptr + i1*index_2_size*index_3_size + i2*index_3_size + i3);
  };
#else
  double& DoubleArray3D::operator()(size_t i1, size_t i2, size_t i3) {
    return *(data_ptr + i1*index_2_size*index_3_size + i2*index_3_size + i3);
  };

  const double& DoubleArray3D::operator()(size_t i1, size_t i2, size_t i3) const {
    return *(data_ptr + i1*index_2_size*index_3_size + i2*index_3_size + i3);
  };
#endif

  bool DoubleArray3D::operator==(const DoubleArray3D& rhs) {
    if (&rhs == this)
      return true;
    if (index_1_size != rhs.index_1_size || index_2_size != rhs.index_2_size || index_3_size != rhs.index_3_size)
      return false;
    for (size_t i = 0; i < index_1_size*index_2_size*index_3_size; ++i) {
      if (*(data_ptr + i) != *(rhs.data_ptr + i))
        return false;
    }
    return true;
  }

  //  Input/Output
  //
  //  Prints out values as as if they were in the first Cartesian 
  //  quadrant --- not in matrix indexing.
  //
  std::ostream&  operator <<(std::ostream& outStream, const DoubleArray3D& A) {
    for (size_t i1 = 0; i1 < A.index_1_size; ++i1) {
      for (size_t i2 = 0; i2 < A.index_2_size; ++i2) {
        for (size_t i3 = 0; i3 < A.index_3_size - 1; ++i3) {
          // outStream << setw(12) << A(i1, i2, i3) << " ";
          outStream << A(i1, i2, i3) << '\t';
        }
        outStream << A(i1, i2, A.index_3_size - 1) << '\n';
      }
      outStream << '\n';
    }
    return outStream;
  }


#ifdef _DEBUG
  bool DoubleArray3D::BoundsCheck(size_t i, size_t begin, size_t end, int coordinate) const {
    if ((i < begin) || (i > end)) {
      std::cerr << "Array index " << coordinate << " out of bounds \n";
      std::cerr << "Offending index value : " << i << " Acceptable Range [" << begin << "," << end << "]\n";
      return false;
    }
    return true;
  }
#else
  inline bool DoubleArray3D::BoundsCheck(size_t, size_t, size_t, int) const { return true; }
#endif
}       /* Namespace Math270A */
