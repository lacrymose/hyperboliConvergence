
# pragma once

# include <parallalg/array.h>

# include <type_traits>

# include <cassert>
# include <cstring>

namespace par
{

/*
 * ------------------------- par::copy ------------------------
 */

/*
 * Copies values from one array to another
 */
   template<typename ElemT, int NDIM, ArraySizing AS1, ArraySizing AS2>
   void copy(       Array<ElemT,NDIM,AS1>& dst,
              const Array<ElemT,NDIM,AS2>& src )
  {
      assert( (dst.shape() == src.shape()) && "par::copy - arrays must be the same shape");

      const size_t len = dst.flattened_length();

      for( size_t i=0; i<len; i++ ){ dst.flatten(i) = src.flatten(i); }

      return;
  }

/*
 * Copies values from one array to another
 *    horrible const cast to use memcpy if ElemT is trivially copyable
 */
   template<typename ElemT, int NDIM, ArraySizing AS1, ArraySizing AS2>
      requires std::is_trivially_copyable_v<ElemT>
   void copy(       Array<ElemT,NDIM,AS1>& dst,
              const Array<ElemT,NDIM,AS2>& src )
  {
      assert( (dst.shape() == src.shape()) && "par::copy - arrays must be the same shape");

      const size_t len = sizeof(ElemT)*dst.flattened_length();

      memcpy( const_cast<ElemT*>(dst.flatten().data()),
                                 src.flatten().data(), len );
      return;
  }


/*
 * Returns a copy of the argument array
 */
   template<typename ElemT, int NDIM, ArraySizing AS>
   Array<ElemT,NDIM,AS> copy( const Array<ElemT,NDIM,AS>& src )
  {
      Array<ElemT,NDIM,AS> dst(src.shape());
      copy( dst, src );
      return dst;
  }


/*
 * ------------------------- par::fill ------------------------
 */

/*
 * Fill an array with a given value
 */
   template<typename ElemT, int NDIM, ArraySizing AS>
      requires std::is_copy_assignable_v<ElemT>
   void fill( Array<ElemT,NDIM,AS>& dst, const ElemT& value )
  {
      const size_t len = dst.flattened_length();

      for( size_t i=0; i<len; i++ ){ dst.flatten(i) = value; }

      return;
  }


/*
 * ------------------------- par::for_each ------------------------
 */

/*
 * Applies a given functor to each element pack in an array pack
 */
   template<typename Functor, int NDIM, typename    ElemT0, ArraySizing    AS0,
                                        typename... ElemTs, ArraySizing... ASs>
   void for_each( const Functor&                   func,
                        Array<ElemT0,NDIM,AS0>&    array0,
                        Array<ElemTs,NDIM,ASs>&... arrays )
  {
      (( assert(   (array0.shape() == arrays.shape())
                && "par::for_each - arrays must be the same shape" ) ),... );

      const size_t len = array0.flattened_length();

      for( size_t i=0; i<len; i++ )
     {
         func( array0.flatten(i),
               arrays.flatten(i)... );
     }
      return;
  }


/*
 * ------------------------- par::transform ------------------------
 */

/*
 * Applies the given functor to each element pack in an array pack, and stores the result in another array
 */

   template<typename Functor, int NDIM, typename    ElemTd, ArraySizing    ASd,
                                        typename    ElemT0, ArraySizing    AS0,
                                        typename... ElemTs, ArraySizing... ASs>
   void transform( const Functor&                   func,
                         Array<ElemTd,NDIM,ASd>&    dst,
                   const Array<ElemT0,NDIM,AS0>&    src0,
                   const Array<ElemTs,NDIM,ASs>&... srcs )
  {
         assert(   (dst.shape() == src0.shape())
                && "par::transform - arrays must be the same shape" );

      (( assert(   (dst.shape() == srcs.shape())
                && "par::transform - arrays must be the same shape" ) ),... );

      const size_t len = dst.flattened_length();
   
      for( size_t i=0; i<len; i++ )
     {
         dst.flatten(i) = func( src0.flatten(i),
                                srcs.flatten(i)... );
     }
      return;
  }

/*
 * Return value overload
 */
   template<typename Functor, int NDIM, typename    ElemT0, ArraySizing    AS0,
                                        typename... ElemTs, ArraySizing... ASs>
   Array<ElemT0,NDIM,AS0> transform( const Functor&                   func,
                                     const Array<ElemT0,NDIM,AS0>&    src0,
                                     const Array<ElemTs,NDIM,ASs>&... srcs )
  {
      Array<ElemT0,NDIM,AS0> dst(src0.shape());
      transform( func, dst, src0, srcs... );
      return dst;
  }


/*
 * ------------------------- par::transform_reduce ------------------------
 */

/*
 * Applies the given functor to each element pack in an array pack, and reduces result with given binary functor
 */

   template<int NDIM, typename ReductionType, typename TransformFunctor, typename    ElemT0, ArraySizing    AS0,
                                              typename    ReduceFunctor, typename... ElemTs, ArraySizing... ASs>
   ReductionType transform_reduce( const TransformFunctor&          tfunc,
                                   const    ReduceFunctor&          rfunc,
                                         ReductionType              init,
                                   const Array<ElemT0,NDIM,AS0>&    src0,
                                   const Array<ElemTs,NDIM,ASs>&... srcs )
  {
      (( assert(   (src0.shape() == srcs.shape())
                && "par::transform_reduce - arrays must be the same shape" ) ),... );

      const size_t len = src0.flattened_length();

      for( size_t i=0; i<len; i++ )
     {
         init = rfunc( std::move(init),
                       tfunc( src0.flatten(i),
                              srcs.flatten(i)... ) );
     }
      return init;
  }
}
