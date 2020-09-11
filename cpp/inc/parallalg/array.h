
# pragma once

# include <utils/type-traits.h>

# include <array>
# include <vector>

# include <cstddef>

namespace par
{

/*
 * ---------------- Atomic types -----------------------------------------
 */

/*
 * An NDIM-dimensional index into an NDIM-dimensional array
 */
   template<int NDIM>
   struct Idx
  {
      std::array<size_t,NDIM> idxs;
      const size_t& operator[]( const unsigned int i ) const { return idxs[i]; }
  };

/*
 * An NDIM-dimensional increment an NDIM-dimensional index. Used for defining stencils
 */
   template<int NDIM>
   struct Offset
  {
      std::array<int,NDIM> offsets;
      const int& operator[]( const unsigned int i ) const { return offsets[i]; }
  };

/*
 * The size in each dimension of an NDIM-dimensional array
 */
   template<int NDIM>
   struct Shape
  {
      std::array<size_t,NDIM> shape;
      const size_t& operator[]( const unsigned int i ) const { return shape[i]; }
  };

/*
 * The strides in memory needed to flatten an NDIM-dimensional array with no padding
 *    last index changes fastest (row-major)
 *
 *    eg in 3D, the element in memory corresponding to index {i,j,k} is: stride[0]*i + stride[1]*j + stride[2]*k
 */
   template<int NDIM>
   struct Stride
  {
      std::array<size_t,NDIM> strides;

      Stride( const Shape<NDIM>& shape )
     {
         strides[NDIM-1]=1;
         for( int i=NDIM-2; i>=0; i-- ){ strides[i]=strides[i+1]*shape[i+1]; }
     }

      const size_t& operator[]( const unsigned int i ) const { return strides[i]; }
  };


/*
 * ---------------- Operator overloads on atomic types -----------------------------------------
 */

/*
 * return index to position 'offset' relative to 'idx'
 */
   template<int NDIM, typename Indices=std::make_index_sequence<NDIM>>
   Idx<NDIM> operator+( const Idx<NDIM> idx, const Offset<NDIM> off )
  {
      Idx<NDIM> result(idx);
      for( unsigned int i=0; i<NDIM; i++ ){ result.idxs[i]+=off[i]; }
      return result;
  }

/*
 * check if two Shapes are equal (sizes in all dimensions are equal)
 */
   template<int NDIM>
   bool operator==( const Shape<NDIM> lhs, const Shape<NDIM> rhs )
  {
      bool isSame=true;
      for( unsigned int i=0; i<NDIM; i++ ){ isSame = isSame && (lhs[i]==rhs[i]); }
      return isSame;
  }

/*
 * calculate the index in the flattened array corresponding to the n-dimensional index 'idx', in an n-dimensional array with shape 'stride'
 */
   template<int NDIM>
   size_t operator*( const Stride<NDIM>& stride, const Idx<NDIM>& idx )
  {
      size_t i=0;
      for( unsigned int m=0; m<NDIM; m++ ){ i+= stride[m]*idx[m]; }
      return i;
  }

   template<int NDIM>
   size_t operator*( const Idx<NDIM>& idx, const Stride<NDIM>& stride )
  {
      return stride*idx;
  }


/*
 * ---------------- Functions on atomic types -----------------------------------------
 */

/*
 * return Shape of node-based array, given Shape of cell-based array
 */
   template<int NDIM>
   Shape<NDIM> nodeDims_from_cellDims( const Shape<NDIM>& cellShape )
  {
      Shape<NDIM> nodeShape(cellShape);
      for( unsigned int i=0; i<NDIM; i++ ){ nodeShape.shape[i]+=1; }
      return nodeShape;
  }

/*
 * return Shape of cell-based array, given Shape of node-based array
 */
   template<int NDIM>
   Shape<NDIM> cellDims_from_nodeDims( const Shape<NDIM>& nodeShape )
  {
      Shape<NDIM> cellShape(nodeShape);
      for( unsigned int i=0; i<NDIM; i++ ){ cellShape.shape[i]-=1; }
      return cellShape;
  }

/*
 * return length of flattened array corresponding to NDIM-dimensional array with given shape
 */
   template<int NDIM>
   size_t length( const Shape<NDIM>& shape )
  {
      size_t l=1;
      for( unsigned int i=0; i<NDIM; i++ ){ l*=shape[i]; }
      return l;
  }


/*
 * ---------------- Memory handling types  -----------------------------------------
 */

/*
 * Flag to set whether an array has fixed size from initialisation, or can be dynamically resized
 */
   enum struct ArraySizing { Static, Dynamic };

/*
 * NDIM-dimensional array type. Last index changes fastest (row-major)
 *    Has 'NDIM' dimensions, and elements of type 'ElemT'
 *    The shape (length in each dimension) must be specified at construction e.g.:
 *       Array<double,2> myArray(5,25);
 *    Elements are accessed using an Idx type of the same dimension, which can be constructed in-place e.g.:
 *       Idx<2> idx{2,3};
 *       auto x = myArray[ idx ];      // use existing Idx instance
 *       auto y = myArray[ {4,11} ];   // construct temporary in-place
 */
   template<typename ElemT, int NDIM, ArraySizing SIZING=ArraySizing::Static>
   struct Array
  {
   private:

   // shape of multi-dimensional array
      Shape<NDIM>   shape_array;

   // place-values for flattening a multi-dimensional index to a 1D index
      Stride<NDIM> stride_array;

   // the flattened 1D array in memory
      std::vector<ElemT> elems_array;

   public:

      using ElemType = ElemT;
      constexpr static int nDim=NDIM;

   // constructors (no construction without dimensions, no copying, only move)
      Array() = delete;
      Array( const Array&  ) = delete;
      Array(       Array&& ) = default;

   // assignment (no copying, only move)
      Array& operator=( const Array&  ) = delete;
      Array& operator=(       Array&& ) = default;

   // construct with specified dimensions, calculate place-value strides, and initialise flattened array to correct size
      template<typename... Ints>
         requires   (sizeof...(Ints)==NDIM)
                 && (is_integer_v<Ints> && ...)
      Array( const Ints... Is ) : shape_array{Is...}, stride_array({Is...}), elems_array( (Is*...) ) {}

      Array( const Shape<NDIM>& s ) : shape_array(s), stride_array(s), elems_array( length(s) ) {}

   // flattened length of array
      size_t flattened_length() const { return elems_array.size(); }

   // accessors
      // elements of array
      const ElemT& operator[]( const Idx<NDIM>& idx ) const { return elems_array[ stride_array*idx ]; }
            ElemT& operator[]( const Idx<NDIM>& idx )       { return elems_array[ stride_array*idx ]; }

      // flattened elements of array
      const ElemT& flatten( const size_t i ) const { return elems_array[i]; }
            ElemT& flatten( const size_t i )       { return elems_array[i]; }

      const std::vector<ElemT>& flatten() const { return elems_array; }

      // array properties
      const size_t& shape(  const unsigned int i ) const { return  shape_array[i]; }
      const size_t& stride( const unsigned int i ) const { return stride_array[i]; }

      const Shape< NDIM>& shape()  const { return  shape_array; }
      const Stride<NDIM>& stride() const { return stride_array; }

   /*
    * Can resize array (same dimensions, different shape) if originally declared as dynamic
    */
      void resize( const Shape<NDIM>& s ) requires (SIZING==ArraySizing::Dynamic)
     {
         shape_array  = Shape< NDIM>(s);
         stride_array = Stride<NDIM>(s);
         elems_array.resize(length(s));
     }
  };

/*
 * returns a vector of arrays with given element type, dimensionality, sizing and shape
 *    utility function using emplace_back because Arrays cannot be copy constructed
 */
   template<typename ElemT, int NDIM, ArraySizing SIZING=ArraySizing::Static>
   std::vector<Array<ElemT,NDIM,SIZING>> vec_of_Arrays( const size_t    narrays,
                                                        const Shape<NDIM> shape )
  {
      std::vector<Array<ElemT,NDIM,SIZING>> vec;
      for( size_t i=0; i<narrays; i++ )
     {
         vec.emplace_back(shape);
     }
      return vec;
  }

}
