#include "double_array_2d.h"


namespace Math270A
{
    //
    //###################################################################
    //                 Constructors/Initialization
    //###################################################################
    //
  DoubleArray2D::DoubleArray2D()
    {
      data_ptr = 0;
      index_1_size = 0;
      index_2_size = 0;
    };

  DoubleArray2D::DoubleArray2D(size_t m, size_t n)
    {
      data_ptr = 0;
      index_1_size = 0;
      index_2_size = 0;
      Initialize(m, n);
    };

  DoubleArray2D::DoubleArray2D(const DoubleArray2D& D)
    {
      data_ptr = 0;
      index_1_size = 0;
      index_2_size = 0;
      Initialize(D);
    };

    DoubleArray2D::~DoubleArray2D()
    {
      if (data_ptr != 0) delete[] data_ptr;
    }


    void DoubleArray2D::Initialize()
    {
      if (data_ptr != 0) delete[] data_ptr;
      data_ptr = 0;
      index_1_size = 0;
      index_2_size = 0;
    };


    void DoubleArray2D::Initialize(const DoubleArray2D& D)
    {
      // Initialization of null object

      if (D.data_ptr == 0)
      {
        Initialize();
        return;
      }

      if ((index_1_size != D.index_1_size) || (index_2_size != D.index_2_size))
      {
        if (data_ptr != 0) delete[] data_ptr;
        index_1_size = D.index_1_size;
        index_2_size = D.index_2_size;
        data_ptr = new double[index_1_size*index_2_size];
      }
      size_t i;
      for (i = 0; i < index_1_size*index_2_size; i++)
      {
        data_ptr[i] = D.data_ptr[i];
      }
    };


    void DoubleArray2D::Initialize(size_t m, size_t n)
    {
      if ((index_1_size != m) || (index_2_size != n))
      {
        if (data_ptr != 0) delete[] data_ptr;
        index_1_size = m;
        index_2_size = n;
        data_ptr = new double[index_1_size*index_2_size];
      }
      size_t i;
      for (i = 0; i < index_1_size*index_2_size; i++)
      {
        data_ptr[i] = 0.0;
      }
    };


    void DoubleArray2D::operator=(const DoubleArray2D& D)
    {
      Initialize(D);
    }

    void DoubleArray2D::operator+=(const DoubleArray2D& G)
    {
      assert(index_1_size == G.index_1_size);
      assert(index_2_size == G.index_2_size);
      for (size_t i1 = 0; i1 < index_1_size; i1++)
      {
        /*#ifdef _OPENMP
        #pragma omp parallel for default(shared) private(i2) schedule(static)
        for (i2 = 0; i2 < index_2_size; i2++)
        {
        *(data_ptr + i2 + i1*index_2_size) = *(data_ptr + i2 + i1*index_2_size) + *(G.data_ptr + i2 + i1*G.index_2_size);
        }
        #else*/
        for (size_t i2 = 0; i2 < index_2_size; i2++)
        {
          *(data_ptr + i2 + i1*index_2_size) = *(data_ptr + i2 + i1*index_2_size) + *(G.data_ptr + i2 + i1*G.index_2_size);
        }
        //#endif
      }
    }

    void DoubleArray2D::operator-=(const DoubleArray2D& G)
    {
      assert(index_1_size == G.index_1_size);
      assert(index_2_size == G.index_2_size);
      for (size_t i1 = 0; i1 < index_1_size; i1++)
        for (size_t i2 = 0; i2 < index_2_size; i2++)
        {
          *(data_ptr + i2 + i1*index_2_size) = *(data_ptr + i2 + i1*index_2_size) - *(G.data_ptr + i2 + i1*G.index_2_size);
        }
    }

    void DoubleArray2D::operator*=(double alpha)
    {
      for (size_t i1 = 0; i1 < index_1_size; i1++)
        for (size_t i2 = 0; i2 < index_2_size; i2++)
        {
          *(data_ptr + i2 + i1*index_2_size) = *(data_ptr + i2 + i1*index_2_size)*alpha;
        }
    }

    void DoubleArray2D::operator/=(double alpha)
    {
      assert(alpha != 0);
      for (size_t i1 = 0; i1 < index_1_size; i1++)
        for (size_t i2 = 0; i2 < index_2_size; i2++)
        {
          *(data_ptr + i2 + i1*index_2_size) = *(data_ptr + i2 + i1*index_2_size) / alpha;
        }
    }

    void DoubleArray2D::SetToValue(double d)
    {
      for (size_t i1 = 0; i1 < index_1_size; i1++)
        for (size_t i2 = 0; i2 < index_2_size; i2++)
        {
          *(data_ptr + i2 + i1*index_2_size) = d;
        }
    }

    double DoubleArray2D::NormInf() const
    {
      double norm = 0;
      for (size_t i1 = 0; i1 < index_1_size; i1++)
        for (size_t i2 = 0; i2 < index_2_size; i2++)
        {
          if (abs(*(data_ptr + i2 + i1*index_2_size)) > norm)
            norm = abs(*(data_ptr + i2 + i1*index_2_size));
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
    double&  DoubleArray2D::operator()(size_t i1, size_t i2)
    {
      assert(BoundsCheck(i1, 0, index_1_size - 1, 1));
      assert(BoundsCheck(i2, 0, index_2_size - 1, 2));
      return *(data_ptr + i2 + i1*index_2_size);
    };

    const double& DoubleArray2D::operator()(size_t i1, size_t i2) const
    {
      assert(BoundsCheck(i1, 0, index_1_size - 1, 1));
      assert(BoundsCheck(i2, 0, index_2_size - 1, 2));
      return *(data_ptr + i2 + i1*index_2_size);
    };
#else
    double&  DoubleArray2D::operator()(size_t i1, size_t i2)
    {
      return *(data_ptr + i2 + i1*index_2_size);
    };

    const double&  DoubleArray2D::operator()(size_t i1, size_t i2) const
    {
      return *(data_ptr + i2 + i1*index_2_size);
    };
#endif

    //
    //###################################################################
    //                Array Structure Access Functions
    //###################################################################
    //
    /// Returns the number of indices in the 1st coordinate direction 

    size_t DoubleArray2D::GetIndex1Size()  const { return index_1_size; }

    /// Returns the number of indices in the 2nd coordinate direction 

    size_t DoubleArray2D::GetIndex2Size()  const { return index_2_size; }

    /// Returns the total number of data values 

    size_t DoubleArray2D::GetDataSize()    const { return index_1_size*index_2_size; }

    /// Returns a pointer to the double array containing the data values 

    double* const DoubleArray2D::GetDataPointer() const { return data_ptr; }

    //  Input/Output
    //
    //  Prints out values as as if they were in the first Cartesian 
    //  quadrant --- not in matrix indexing. 
    //
    //
    std::ostream&  operator <<(std::ostream& out_stream, const DoubleArray2D& A)
    {
      for (size_t i = 0; i < A.index_1_size; ++i) {
        for (size_t j = 0; j < A.index_2_size - 1; ++j) {
          //out_stream << setw(12) << A(i, j) << " ";
          out_stream << std::setw(12) << A(i, j) << '\t';
        }
        out_stream << std::setw(12) << A(i, A.index_2_size - 1) << '\n';
      }
      return out_stream;
    }

    //###################################################################
    //                      Bounds Checking
    //###################################################################
    //
#ifdef _DEBUG
    bool DoubleArray2D::BoundsCheck(size_t i, size_t begin, size_t end, int coordinate) const
    {
      if ((i < begin) || (i  > end))
      {
        std::cerr << "Array index " << coordinate << " out of bounds " << std::endl;
        std::cerr << "Offending index value : " << i << " Acceptable Range [" << begin << "," << end << "]" << std::endl;
        return false;
      }
      return true;
    }
#else
    bool DoubleArray2D::BoundsCheck(size_t, size_t, size_t, int) const { return true; }
#endif

}       /* Namespace Math270A */
