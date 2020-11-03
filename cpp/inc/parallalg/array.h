
# pragma once

# include <parallalg/parallalg.h>

# include <utils/type-traits.h>

# include <array>
# include <vector>

namespace par
{

/*
 * ---------------- index manipulation types for Array -----------------------------------------
 */

/*
 * An NDIM-dimensional increment an NDIM-dimensional index. Used for defining stencils
 */
   template<int          NDIM,
            GridType GRIDTYPE= Primal>
   struct Offset
  {
      std::array<int,NDIM> offsets;
      const int& operator[]( const unsigned int i ) const { return offsets[i]; }
  };

/*
 * The size in each dimension of an NDIM-dimensional array
 */
   template<int          NDIM,
            GridType GRIDTYPE= Primal>
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
   template<int          NDIM,
            GridType GRIDTYPE= Primal>
   struct Stride
  {
      std::array<size_t,NDIM> strides;

      Stride( const Shape<NDIM,GRIDTYPE>& shape )
     {
         strides[NDIM-1]=1;
         for( int i=NDIM-2; i>=0; i-- ){ strides[i]=strides[i+1]*shape[i+1]; }
     }

      const size_t& operator[]( const unsigned int i ) const { return strides[i]; }
  };


/*
 * ---------------- Convenience typedefs ---------------------
 */

// Offset types

   // grid type
   template<int NDIM>
   using PrimalOffset = Offset<NDIM,Primal>;

   template<int NDIM>
   using DualOffset = Offset<NDIM,Dual>;

   // dimension
   template<GridType GRIDTYPE>
   using Offset1 = Offset<1,GRIDTYPE>;

   template<GridType GRIDTYPE>
   using Offset2 = Offset<2,GRIDTYPE>;

   template<GridType GRIDTYPE>
   using Offset3 = Offset<3,GRIDTYPE>;

   // grid type and dimension
   using PrimalOffset1 = Offset<1,Primal>;
   using PrimalOffset2 = Offset<2,Primal>;
   using PrimalOffset3 = Offset<3,Primal>;

   using DualOffset1 = Offset<1,Dual>;
   using DualOffset2 = Offset<2,Dual>;
   using DualOffset3 = Offset<3,Dual>;

// Shape types

   // grid type
   template<int NDIM>
   using PrimalShape = Shape<NDIM,Primal>;

   template<int NDIM>
   using DualShape = Shape<NDIM,Dual>;

   // dimension
   template<GridType GRIDTYPE>
   using Shape1 = Shape<1,GRIDTYPE>;

   template<GridType GRIDTYPE>
   using Shape2 = Shape<2,GRIDTYPE>;

   template<GridType GRIDTYPE>
   using Shape3 = Shape<3,GRIDTYPE>;

   // grid type and dimension
   using PrimalShape1 = Shape<1,Primal>;
   using PrimalShape2 = Shape<2,Primal>;
   using PrimalShape3 = Shape<3,Primal>;

   using DualShape1 = Shape<1,Dual>;
   using DualShape2 = Shape<2,Dual>;
   using DualShape3 = Shape<3,Dual>;

// Stride types

   // grid type
   template<int NDIM>
   using PrimalStride = Stride<NDIM,Primal>;

   template<int NDIM>
   using DualStride = Stride<NDIM,Dual>;

   // dimension
   template<GridType GRIDTYPE>
   using Stride1 = Stride<1,GRIDTYPE>;

   template<GridType GRIDTYPE>
   using Stride2 = Stride<2,GRIDTYPE>;

   template<GridType GRIDTYPE>
   using Stride3 = Stride<3,GRIDTYPE>;

   // grid type and dimension
   using PrimalStride1 = Stride<1,Primal>;
   using PrimalStride2 = Stride<2,Primal>;
   using PrimalStride3 = Stride<3,Primal>;

   using DualStride1 = Stride<1,Dual>;
   using DualStride2 = Stride<2,Dual>;
   using DualStride3 = Stride<3,Dual>;

/*
 * ---------------- Operator overloads on atomic types -----------------------------------------
 */

/*
 * return index to position 'offset' relative to 'idx'
 */
   template<int          NDIM,
            GridType GRIDTYPE>
   Idx<NDIM,GRIDTYPE> operator+( const Idx<   NDIM,GRIDTYPE>& idx,
                                 const Offset<NDIM,GRIDTYPE>& off )
  {
      Idx<NDIM,GRIDTYPE> result(idx);
      for( unsigned int i=0; i<NDIM; i++ ){ result.idxs[i]+=off[i]; }
      return result;
  }

/*
 * check if two Shapes are equal (sizes in all dimensions are equal)
 */
   template<int          NDIM,
            GridType GRIDTYPE>
   bool operator==( const Shape<NDIM,GRIDTYPE>& lhs,
                    const Shape<NDIM,GRIDTYPE>& rhs )
  {
      bool isSame=true;
      for( unsigned int i=0; i<NDIM; i++ ){ isSame = isSame && (lhs[i]==rhs[i]); }
      return isSame;
  }

/*
 * calculate the index in the flattened array corresponding to the n-dimensional index 'idx', in an n-dimensional array with shape 'stride'
 */
   template<int          NDIM,
            GridType GRIDTYPE>
   size_t operator*( const Stride<NDIM,GRIDTYPE>& stride,
                     const Idx<   NDIM,GRIDTYPE>&    idx )
  {
      size_t i=0;
      for( unsigned int m=0; m<NDIM; m++ ){ i+= stride[m]*idx[m]; }
      return i;
  }

   template<int          NDIM,
            GridType GRIDTYPE>
   size_t operator*( const Idx<   NDIM,GRIDTYPE>&    idx,
                     const Stride<NDIM,GRIDTYPE>& stride )
  {
      return stride*idx;
  }


/*
 * ---------------- Functions on atomic types -----------------------------------------
 */

/*
 * return Shape of primal array, given Shape of dual array
 */
   template<int NDIM>
   PrimalShape<NDIM> primalShape( const DualShape<NDIM>& dual_shape )
  {
      PrimalShape<NDIM> primal_shape;
      for( unsigned int i=0; i<NDIM; ++i ){ primal_shape.shape[i]=dual_shape[i]+1; }
      return primal_shape;
  }

/*
 * return Shape of dual array, given Shape of primal array
 */
   template<int NDIM>
   DualShape<NDIM> dualShape( const PrimalShape<NDIM>& primal_shape )
  {
      DualShape<NDIM> dual_shape;
      for( unsigned int i=0; i<NDIM; i++ ){ dual_shape.shape[i]=primal_shape[i]-1; }
      return dual_shape;
  }

/*
 * return length of flattened array corresponding to NDIM-dimensional array with given shape
 */
   template<int          NDIM,
            GridType GRIDTYPE>
   size_t length( const Shape<NDIM,GRIDTYPE>& shape )
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
   enum struct ArraySizing { Fixed, Dynamic };
   constexpr ArraySizing   FixedSize = ArraySizing::Fixed;
   constexpr ArraySizing DynamicSize = ArraySizing::Dynamic;

/*
 * NDIM-dimensional array type. Last index changes fastest (row-major)
 *    Has 'NDIM' dimensions, and elements of type 'ElemT'
 *    The shape (length in each dimension) must be specified at construction e.g.:
 *       Array<double,2> myArray(5,25);
 *    Elements are accessed using an Idx type of the same dimension, which can be constructed in-place e.g.:
 *       Idx<2> idx{2,3};
 *       auto x = myArray[ idx ];      // use existing Idx instance
 *       auto y = myArray[ {4,11} ];   // construct temporary in-place
 *    GridType defines the grid layout, either Primal or Dual.
 *       Defaulted to Primal because only matters if arrays of both type are used?
 *    ArraySizing defines whether the array extents can be changed after construction
 *       Fixed by default.
 *       resize method only available if Dynamic
 */
   template<typename       ElemT,
            int             NDIM,
            GridType    GRIDTYPE= Primal,
            ArraySizing   SIZING= FixedSize>
   struct Array
  {
   public: /* typedefs and static members */

      using ElemType = ElemT;
      constexpr static int         nDim        = NDIM;
      constexpr static GridType    gridType    = GRIDTYPE;
      constexpr static ArraySizing arraySizing = SIZING;

      using IdxType    = Idx<   NDIM,GRIDTYPE>;
      using OffsetType = Offset<NDIM,GRIDTYPE>;
      using ShapeType  = Shape< NDIM,GRIDTYPE>;
      using StrideType = Stride<NDIM,GRIDTYPE>;

   private: /* invariant members and underlying memory */

   // shape of multi-dimensional array
      ShapeType   shape_array;

   // place-values for flattening a multi-dimensional index to a 1D index
      StrideType stride_array;

   // the flattened 1D array in memory
      std::vector<ElemT> elems_array;

   public:

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
                 && (is_integer_v<Ints>&&...)
      Array( const Ints... Is ) : shape_array{Is...}, stride_array({Is...}), elems_array( (Is*...) ) {}

   // constructor array with shape s, and no initialisation of array elements
      Array( const ShapeType& s ) : shape_array(s), stride_array(s), elems_array( length(s) ) {}

   // constructor array with shape s, and initialise array elements to val
      Array( const ShapeType& s, const ElemT& val ) : shape_array(s), stride_array(s), elems_array( length(s), val ) {}

   // flattened length of array
      size_t flattened_length() const { return elems_array.size(); }

   // accessors
      // elements of array
      const ElemT& operator()( const IdxType& idx ) const { return elems_array[ stride_array*idx ]; }
            ElemT& operator()( const IdxType& idx )       { return elems_array[ stride_array*idx ]; }

      // flattened elements of array
      const ElemT& flatten( const size_t i ) const { return elems_array[i]; }
            ElemT& flatten( const size_t i )       { return elems_array[i]; }

      const std::vector<ElemT>& flatten() const { return elems_array; }

      // array properties
      const size_t& shape(  const unsigned int i ) const { return  shape_array[i]; }
      const size_t& stride( const unsigned int i ) const { return stride_array[i]; }

      const ShapeType&   shape() const { return  shape_array; }
      const StrideType& stride() const { return stride_array; }

   /*
    * Can resize array (same dimensions, different shape) if originally declared as dynamic
    */
      void resize( const ShapeType& s ) requires (SIZING==DynamicSize)
     {
         shape_array  =  ShapeType(s);
         stride_array = StrideType(s);
         elems_array.resize(length(s));
     }
  };

/*
 * returns a vector of arrays with given element type, dimensionality, sizing and shape
 *    utility function using emplace_back because Arrays cannot be copy constructed
 */
   template<typename       ElemT,
            int           NDIM,
            GridType    GRIDTYPE= Primal,
            ArraySizing   SIZING= FixedSize>
   auto vec_of_Arrays( const size_t             narrays,
                       const Shape<NDIM,GRIDTYPE> shape )
  {
      std::vector<Array<ElemT,NDIM,GRIDTYPE,SIZING>> vec;
      for( size_t i=0; i<narrays; i++ )
     {
         vec.emplace_back(shape);
     }
      return vec;
  }

/*
 * ---------------- Convenience Array typedefs -----------------------------------------
 */

   // grid type
   template<typename       ElemT,
            int             NDIM,
            ArraySizing   SIZING= FixedSize>
   using PrimalArray = Array<ElemT,NDIM,Primal,SIZING>;

   template<typename       ElemT,
            int             NDIM,
            ArraySizing   SIZING= FixedSize>
   using DualArray = Array<ElemT,NDIM,Dual,SIZING>;

   // dimension
   template<typename       ElemT,
            GridType    GRIDTYPE= Primal,
            ArraySizing   SIZING= FixedSize>
   using Array1 = Array<ElemT,1,GRIDTYPE,SIZING>;

   template<typename       ElemT,
            GridType    GRIDTYPE= Primal,
            ArraySizing   SIZING= FixedSize>
   using Array2 = Array<ElemT,2,GRIDTYPE,SIZING>;

   template<typename       ElemT,
            GridType    GRIDTYPE= Primal,
            ArraySizing   SIZING= FixedSize>
   using Array3 = Array<ElemT,3,GRIDTYPE,SIZING>;

   // grid type and dimension
   template<typename       ElemT,
            ArraySizing   SIZING= FixedSize>
   using PrimalArray1 = Array<ElemT,1,Primal,SIZING>;

   template<typename       ElemT,
            ArraySizing   SIZING= FixedSize>
   using PrimalArray2 = Array<ElemT,2,Primal,SIZING>;

   template<typename       ElemT,
            ArraySizing   SIZING= FixedSize>
   using PrimalArray3 = Array<ElemT,3,Primal,SIZING>;

   template<typename       ElemT,
            ArraySizing   SIZING= FixedSize>
   using DualArray1 = Array<ElemT,1,Dual,SIZING>;

   template<typename       ElemT,
            ArraySizing   SIZING= FixedSize>
   using DualArray2 = Array<ElemT,2,Dual,SIZING>;

   template<typename       ElemT,
            ArraySizing   SIZING= FixedSize>
   using DualArray3 = Array<ElemT,3,Dual,SIZING>;

}
