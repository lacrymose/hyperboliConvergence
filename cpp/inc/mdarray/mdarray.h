
# pragma once

# include <utils/type-traits.h>

# include <array>
# include <vector>

# include <cstddef>

/*
 * data type to hold an n-dimensional index
 */
   template<int nDim>
   struct Idx
  {
      const std::array<size_t,nDim> idxs;
      const size_t& operator[]( const int i ) const { return idxs[i]; }
  };


/*
 * Offset to increment an n-dimensional index. Used for defining stencils
 */
   template<int nDim>
   struct Offset
  {
      const std::array<size_t,nDim> off;
      const size_t& operator[]( const int i ) const { return off[i]; }
  };


/*
 * return index in position 'offset' relative to 'idx'
 */
   template<int nDim>
   Idx<nDim> operator+( const Idx<nDim> idx, const Offset<nDim> off )
  {
      Idx<nDim> result(idx);
      for( int i=0; i<nDim; i++ ){ result[i]+=off[i]; }
      return result;
  }


/*
 * data type to hold the dimensions of an n-dimensional array
 */
   template<int nDim>
   struct Dims
  {
      const std::array<size_t,nDim> dims;
      const size_t& operator[]( const int i ) const { return dims[i]; }
  };


/*
 * data type to hold the strides in memory needed to flatten an n-dimensional array
 *    last index changes fastest (row-major)
 *    eg in 3D, the element in memory corresponding to index i,j,k is: stride[0]*i + stride[1]*j + stride[2]*k
 */
   template<int nDim>
   struct Stride
  {
      const std::array<size_t,nDim> stride;

      Stride( const Dims<nDim>& dims )
     {
         stride[nDim-1]=1;
         for( int i=nDim-2; i>=0; i-- )
        {
            stride[i]=stride[i+1]*dims[i+1];
        }
     }

      const size_t& operator[]( const int i ) const { return stride[i]; }
  };


/*
 * calculate the index in the flattened array corresponding n-dimensional index 'idx', in an n-dimensional array with shape 'stride'
 */
   template<int nDim>
   size_t operator*( const Stride<nDim>& stride, const Idx<nDim>& idx )
  {
      size_t i=0;
      for( int m=0; m<nDim; m++ ){ i+= stride[m]*idx[m]; }
      return i;
  }


/*
 * Multi-dimensional array type. Last index changes fastest (row-major)
 *    Has 'nDim' dimensions, and elements of type 'ElemT'
 *    The shape (length in each dimension) must be specified at construction e.g.:
 *       MDArray<double,2> myArray(5,25);
 *    Elements are accessed using an Idx type of the same dimension, e.g.:
 *       Idx<2> idx{2,3};
 *       auto x = myArray[ idx ];
 *       auto y = myArray[ {4,11} ];
 */
   template<typename ElemT, int nDim>
   struct MDArray
  {
   // shape of multi-dimensional array
      Dims<  nDim>   dims;

   // place-values for flattening a multi-dimensional index to a 1D index
      Stride<nDim> stride;

   // the flattened 1D array in memory
      std::vector<ElemT> elems;

   // constructors (no blank construction or copying, only move)
      MDArray() = delete;
      MDArray( const MDArray&  ) = delete; /**/
      MDArray(       MDArray&& ) = default;

   // assignment (no copying, only move)
      MDArray& operator=( const MDArray&  ) = delete; /**/
      MDArray& operator=(       MDArray&& ) = default;

   // construct with specified dimensions, calculate place-value strides, and resize flattened array
      template<typename... Ints>
         requires   (sizeof...(Ints)==nDim)
                 && (is_integer_v<Ints> && ...)
      MDArray( Ints... Is ) : dims{Is...}, stride({Is...}), elems( (Is*...) ) {}

   // accessors
      const ElemT& operator[]( const Idx<nDim>& idx ) const { return elems[ stride*idx ]; }
            ElemT& operator[]( const Idx<nDim>& idx )       { return elems[ stride*idx ]; }
  };

