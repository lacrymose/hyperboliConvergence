
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
      std::array<size_t,nDim> idxs;
      const size_t& operator[]( const int i ) const { return idxs[i]; }
  };


/*
 * Offset to increment an n-dimensional index. Used for defining stencils
 */
   template<int nDim>
   struct Offset
  {
      std::array<int,nDim> offsets;
      const int& operator[]( const int i ) const { return offsets[i]; }
  };


/*
 * return index in position 'offset' relative to 'idx'
 *    Offload to helper function with integer pack to access idx & off elements
 *       Idx is a const class, so cannot create new Idx object then add offset
 *       Cannot hard-code brace initialisation because do not know nDim a-priori
 */
   template<int nDim, typename Indices=std::make_index_sequence<nDim>>
   Idx<nDim> operator+( const Idx<nDim> idx, const Offset<nDim> off )
  {
      return idx_add_offset( idx, off, Indices{} );
  }

/*
 * helper function to implement addition of Idx and Offset
 *    uses integer sequence to brace-initialise new Idx. This avoids having to construct
 */
   template<int nDim, size_t... Is>
      requires (sizeof...(Is)==nDim)
   Idx<nDim> idx_add_offset( const Idx<nDim> idx, const Offset<nDim> off, std::index_sequence<Is...> )
  {
      return {(idx[Is]+off[Is])...};
  }


/*
 * data type to hold the dimensions of an n-dimensional array
 */
   template<int nDim>
   struct Dims
  {
      std::array<size_t,nDim> dims;
      const size_t& operator[]( const int i ) const { return dims[i]; }
  };

/*
 * check if all elements of Dims are equal
 */
   template<int nDim>
   bool operator==( const Dims<nDim> lhs, const Dims<nDim> rhs )
  {
      bool isSame=true;
      for( int i=0; i<nDim; i++ ){ isSame = isSame && (lhs[i]==rhs[i]); }
      return isSame;
  }


/*
 * helper function to calculate the strides through a flattened array corresponding to elements of a multidimensional array with shape dims
 */
   template<int nDim>
   std::array<size_t,nDim> make_strides( const Dims<nDim>& dims )
  {
      std::array<size_t,nDim> strides;
      strides[nDim-1]=1;
      for( int i=nDim-2; i>=0; i-- )
     {
         strides[i]=strides[i+1]*dims[i+1];
     }
      return strides;
  }

/*
 * return length of flattened array corresponding to multi-dimensional array with dimensions dim
 */
   template<int nDim>
   size_t length( const Dims<nDim>& dims )
  {
      size_t l=1;
      for( int i=0; i<nDim; i++ ){ l*=dims[i]; }
      return l;
  }

/*
 * data type to hold the strides in memory needed to flatten an n-dimensional array with no padding
 *    last index changes fastest (row-major)
 *    eg in 3D, the element in memory corresponding to index i,j,k is: stride[0]*i + stride[1]*j + stride[2]*k
 */
   template<int nDim>
   struct Stride
  {
      std::array<size_t,nDim> strides;

      Stride( const Dims<nDim>& dims ) : strides(make_strides(dims)) {}

      const size_t& operator[]( const int i ) const { return strides[i]; }
  };


/*
 * calculate the index in the flattened array corresponding to the n-dimensional index 'idx', in an n-dimensional array with shape 'stride'
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
 *    Elements are accessed using an Idx type of the same dimension, which can be constructed in-place e.g.:
 *       Idx<2> idx{2,3};
 *       auto x = myArray[ idx ];      // use existing Idx instance
 *       auto y = myArray[ {4,11} ];   // construct temporary in-place
 */
   template<typename ElemT, int nDim>
   struct MDArray
  {
   // shape of multi-dimensional array
      Dims<nDim>   dims;

   // place-values for flattening a multi-dimensional index to a 1D index
      Stride<nDim> stride;

   // the flattened 1D array in memory
      std::vector<ElemT> elems;

   // constructors (no construction without dimensions, no copying, only move)
      MDArray() = delete;
      MDArray( const MDArray&  ) = delete;
      MDArray(       MDArray&& ) = default;

   // assignment (no copying, only move)
      MDArray& operator=( const MDArray&  ) = delete;
      MDArray& operator=(       MDArray&& ) = default;

   // construct with specified dimensions, calculate place-value strides, and initialise flattened array to correct size
      template<typename... Ints>
         requires   (sizeof...(Ints)==nDim)
                 && (is_integer_v<Ints> && ...)
      MDArray( const Ints... Is ) : dims{Is...}, stride({Is...}), elems( (Is*...) ) {}

      MDArray( const Dims<nDim>& d ) : dims(d), stride(d), elems( length(dims) ) {}

   // accessors
      const ElemT& operator[]( const Idx<nDim>& idx ) const { return elems[ stride*idx ]; }
            ElemT& operator[]( const Idx<nDim>& idx )       { return elems[ stride*idx ]; }
  };


