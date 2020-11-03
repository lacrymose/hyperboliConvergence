
# pragma once

# include <parallalg/array.h>
# include <parallalg/parallalg.h>

# include <type_traits>

# include <cassert>
# include <cstring>

# include <omp.h>

namespace par
{

/*
 * ------------------------- par::copy ------------------------
 */

/*
 * if no policy provided, use sequential
 */
   template<typename  ElemT,
            int        NDIM,
            GridType     GT,
            ArraySizing AS0,
            ArraySizing AS1>
   void copy(       Array<ElemT,NDIM,GT,AS0>& dst,
              const Array<ElemT,NDIM,GT,AS1>& src )
  {
      copy( execution::seq, dst, src );
  }

/*
 * Return value overloads with/without policy
 */
   template<typename ElemT,
            int       NDIM,
            GridType    GT,
            ArraySizing AS>
   auto copy( const Array<ElemT,NDIM,GT,AS>& src )
  {
      return copy( execution::seq, src );
  }

   template<execution_policy Policy,
            typename          ElemT,
            int                NDIM,
            GridType             GT,
            ArraySizing          AS>
   auto copy(       Policy policy,
              const Array<ElemT,NDIM,GT,AS>& src )
  {
      Array<ElemT,NDIM,GT,AS> dst(src.shape());
      copy( policy, dst, src );
      return dst;
  }

/*
 * memcopy overload for trivial types
 *    horrible const cast to use memcpy if ElemT is trivially copyable
 */
   template<execution_policy Policy,
            typename          ElemT,
            int                NDIM,
            GridType             GT,
            ArraySizing         AS0,
            ArraySizing         AS1>
      requires std::is_trivially_copyable_v<ElemT>
   void copy(       Policy,
                    Array<ElemT,NDIM,GT,AS0>& dst,
              const Array<ElemT,NDIM,GT,AS1>& src )
  {
      assert( (dst.shape() == src.shape()) && "par::copy - arrays must be the same shape");

      const size_t len = sizeof(ElemT)*dst.flattened_length();

      memcpy( const_cast<ElemT*>(dst.flatten().data()),
                                 src.flatten().data(), len );
      return;
  }


/*
 * Serial execution
 */
   template<typename  ElemT,
            int        NDIM,
            GridType     GT,
            ArraySizing AS0,
            ArraySizing AS1>
      requires !(std::is_trivially_copyable_v<ElemT>)
   void copy(       execution::serial_policy,
                    Array<ElemT,NDIM,GT,AS0>& dst,
              const Array<ElemT,NDIM,GT,AS1>& src )
  {
      assert( (dst.shape() == src.shape()) && "par::copy - arrays must be the same shape");

      const size_t len = dst.flattened_length();

      for( size_t i=0; i<len; ++i ){ dst.flatten(i) = src.flatten(i); }

      return;
  }


/*
 * OpenMP execution
 */
   template<typename  ElemT,
            int        NDIM,
            GridType     GT,
            ArraySizing AS0,
            ArraySizing AS1>
      requires !(std::is_trivially_copyable_v<ElemT>)
   void copy(       execution::openmp_policy,
                    Array<ElemT,NDIM,GT,AS0>& dst,
              const Array<ElemT,NDIM,GT,AS1>& src )
  {
      assert( (dst.shape() == src.shape()) && "par::copy - arrays must be the same shape");

      const size_t len = dst.flattened_length();

# ifdef _OPENMP
   # pragma omp parallel for
# endif
      for( size_t i=0; i<len; ++i ){ dst.flatten(i) = src.flatten(i); }

      return;
  }


/*
 * ------------------------- par::fill ------------------------
 */

/*
 * Fill an array with a given value
 */
   template<typename ElemT,
            int       NDIM,
            GridType    GT,
            ArraySizing AS>
      requires std::is_copy_assignable_v<ElemT>
   void fill(       Array<ElemT,NDIM,GT,AS>& dst,
              const ElemT&                 value )
  {
      fill( execution::seq, dst, value );
      return;
  }

/*
 * Serial execution
 */
   template<typename ElemT,
            int       NDIM,
            GridType    GT,
            ArraySizing AS>
      requires std::is_copy_assignable_v<ElemT>
   void fill(       execution::serial_policy,
                    Array<ElemT,NDIM,GT,AS>& dst,
              const ElemT&                 value )
  {
      const size_t len = dst.flattened_length();

      for( size_t i=0; i<len; ++i ){ dst.flatten(i) = value; }

      return;
  }

/*
 * OpenMP execution
 */
   template<typename ElemT,
            int       NDIM,
            GridType    GT,
            ArraySizing AS>
      requires std::is_copy_assignable_v<ElemT>
   void fill(       execution::openmp_policy,
                    Array<ElemT,NDIM,GT,AS>& dst,
              const ElemT&                 value )
  {
      const size_t len = dst.flattened_length();

# ifdef _OPENMP
   # pragma omp parallel for
# endif
      for( size_t i=0; i<len; ++i ){ dst.flatten(i) = value; }

      return;
  }


/*
 * ------------------------- par::generate ------------------------
 */

/*
 * Fill an array with the result of repeated calls to Generator
 */
   template<typename ElemT,
            int       NDIM,
            GridType    GT,
            ArraySizing AS,
            typename Generator>
   void generate( Array<ElemT,NDIM,GT,AS>& dst,
                  Generator          generator )
  {
      generate( execution::seq, dst, generator );
      return;
  }


/*
 * Serial execution
 */
   template<typename ElemT,
            int       NDIM,
            GridType    GT,
            ArraySizing AS,
            typename Generator>
   void generate( execution::serial_policy,
                  Array<ElemT,NDIM,GT,AS>& dst,
                  Generator          generator )
  {
      const size_t len = dst.flattened_length();

      for( size_t i=0; i<len; ++i ){ dst.flatten(i) = generator(); }

      return;
  }


/*
 * OpenMP execution
 */
   template<typename ElemT,
            int       NDIM,
            GridType    GT,
            ArraySizing AS,
            typename Generator>
   void generate( execution::openmp_policy,
                  Array<ElemT,NDIM,GT,AS>& dst,
                  Generator          generator )
  {
      const size_t len = dst.flattened_length();

# ifdef _OPENMP
   # pragma omp parallel for
# endif
      for( size_t i=0; i<len; ++i ){ dst.flatten(i) = generator(); }

      return;
  }


/*
 * ------------------------- par::generate_idx ------------------------
 */

/*
 * Fill an array with the result of repeated calls to Generator, using the index as an argument
 */
   template<typename     ElemT,
            int           NDIM,
            GridType        GT,
            ArraySizing     AS,
            typename Generator>
   void generate_idx( Array<ElemT,NDIM,GT,AS>& dst,
                      Generator          generator )
  {
      generate_idx( execution::seq, dst, generator );
  }

/*
 * 1D serial execution
 */
   template<typename     ElemT,
            GridType        GT,
            ArraySizing     AS,
            typename Generator>
   void generate_idx( execution::serial_policy,
                      Array1<ElemT,GT,AS>& dst,
                      Generator      generator )
  {
      const size_t n = dst.shape(0);
      for( size_t i=0; i<n; ++i ){ dst({i}) = generator( Idx1<GT>{i} ); }
      return;
  }

/*
 * 1D OpenMP execution
 */
   template<typename     ElemT,
            GridType        GT,
            ArraySizing     AS,
            typename Generator>
   void generate_idx( execution::openmp_policy,
                      Array1<ElemT,GT,AS>& dst,
                      Generator      generator )
  {
      const size_t n = dst.shape(0);

# ifdef _OPENMP
   # pragma omp parallel for
# endif
      for( size_t i=0; i<n; ++i ){ dst({i}) = generator( Idx1<GT>{i} ); }
      return;
  }

/*
 * 2D serial execution
 */
   template<typename     ElemT,
            GridType        GT,
            ArraySizing     AS,
            typename Generator>
   void generate_idx( execution::serial_policy,
                      Array2<ElemT,GT,AS>& dst,
                      Generator      generator )
  {
      const size_t ni = dst.shape(0);
      const size_t nj = dst.shape(1);

      for( size_t i=0; i<ni; ++i )
     {
         for( size_t j=0; j<nj; ++j )
        {
            const Idx2<GT> idx{i,j};
            dst(idx) = generator(idx);
        }
     }
      return;
  }

/*
 * 2D OpenMP execution
 */
   template<typename     ElemT,
            GridType        GT,
            ArraySizing     AS,
            typename Generator>
   void generate_idx( execution::openmp_policy,
                      Array2<ElemT,GT,AS>& dst,
                      Generator      generator )
  {
      const size_t ni = dst.shape(0);
      const size_t nj = dst.shape(1);

# ifdef _OPENMP
   # pragma omp parallel for
# endif
      for( size_t i=0; i<ni; ++i )
     {
         for( size_t j=0; j<nj; ++j )
        {
            const Idx2<GT> idx{i,j};
            dst(idx) = generator(idx);
        }
     }
      return;
  }

/*
 * 3D serial execution
 */
   template<typename     ElemT,
            GridType        GT,
            ArraySizing     AS,
            typename Generator>
   void generate_idx( execution::serial_policy,
                      Array3<ElemT,GT,AS>& dst,
                      Generator      generator )
  {
      const size_t ni = dst.shape(0);
      const size_t nj = dst.shape(1);
      const size_t nk = dst.shape(2);

      for( size_t i=0; i<ni; ++i )
     {
         for( size_t j=0; j<nj; ++j )
        {
            for( size_t k=0; k<nk; ++k )
           {
               const Idx3<GT> idx{i,j,k};
               dst(idx) = generator(idx);
           }
        }
     }
      return;
  }

/*
 * 3D OpenMP execution
 */
   template<typename     ElemT,
            GridType        GT,
            ArraySizing     AS,
            typename Generator>
   void generate_idx( execution::openmp_policy,
                      Array3<ElemT,GT,AS>& dst,
                      Generator      generator )
  {
      const size_t ni = dst.shape(0);
      const size_t nj = dst.shape(1);
      const size_t nk = dst.shape(2);

# ifdef _OPENMP
   # pragma omp parallel for
# endif
      for( size_t i=0; i<ni; ++i )
     {
         for( size_t j=0; j<nj; ++j )
        {
            for( size_t k=0; k<nk; ++k )
           {
               const Idx3<GT> idx{i,j,k};
               dst(idx) = generator(idx);
           }
        }
     }
      return;
  }


/*
 * ------------------------- par::for_each ------------------------
 */

/*
 * Applies a given functor to each element pack in an array pack
 */
   template<typename   FuncObj,
            int           NDIM,
            typename    ElemT0,
            typename... ElemTs,
            GridType        GT,
            ArraySizing    AS0,
            ArraySizing... ASs>
   void for_each( FuncObj                         func,
                  Array<ElemT0,NDIM,GT,AS0>&    array0,
                  Array<ElemTs,NDIM,GT,ASs>&... arrays )
  {
      for_each( execution::seq, func, array0, arrays... );
      return;
  }

/*
 * Serial execution
 */
   template<typename   FuncObj,
            int           NDIM,
            typename    ElemT0,
            typename... ElemTs,
            GridType        GT,
            ArraySizing    AS0,
            ArraySizing... ASs>
   void for_each( execution::serial_policy,
                  FuncObj                         func,
                  Array<ElemT0,NDIM,GT,AS0>&    array0,
                  Array<ElemTs,NDIM,GT,ASs>&... arrays )
  {
      (( assert(   (array0.shape() == arrays.shape())
                && "par::for_each - arrays must be the same shape" ) ),... );

      const size_t len = array0.flattened_length();

      for( size_t i=0; i<len; ++i )
     {
         func( array0.flatten(i),
               arrays.flatten(i)... );
     }
      return;
  }

/*
 * OpenMP execution
 */
   template<typename   FuncObj,
            int           NDIM,
            typename    ElemT0,
            typename... ElemTs,
            GridType        GT,
            ArraySizing    AS0,
            ArraySizing... ASs>
   void for_each( execution::openmp_policy,
                  FuncObj                         func,
                  Array<ElemT0,NDIM,GT,AS0>&    array0,
                  Array<ElemTs,NDIM,GT,ASs>&... arrays )
  {
      (( assert(   (array0.shape() == arrays.shape())
                && "par::for_each - arrays must be the same shape" ) ),... );

      const size_t len = array0.flattened_length();

# ifdef _OPENMP
   # pragma omp parallel for
# endif
      for( size_t i=0; i<len; ++i )
     {
         func( array0.flatten(i),
               arrays.flatten(i)... );
     }
      return;
  }


/*
 * ------------------------- par::for_each_idx ------------------------
 */

/*
 * Applies a given functor to each element pack in an array pack, and also passes the current index
 */
   template<typename   FuncObj,
            int           NDIM,
            typename    ElemT0,
            typename... ElemTs,
            GridType        GT,
            ArraySizing    AS0,
            ArraySizing... ASs>
   void for_each_idx( FuncObj                         func,
                      Array<ElemT0,NDIM,GT,AS0>&    array0,
                      Array<ElemTs,NDIM,GT,ASs>&... arrays )
  {
      for_each_idx( execution::seq, func, array0, arrays... );
  }

   template<typename   FuncObj,
            int           NDIM,
            typename    ElemT0,
            typename... ElemTs,
            GridType        GT,
            ArraySizing    AS0,
            ArraySizing... ASs>
   void for_each_idx(       FuncObj                         func,
                      const Array<ElemT0,NDIM,GT,AS0>&    array0,
                      const Array<ElemTs,NDIM,GT,ASs>&... arrays )
  {
      for_each_idx( execution::seq, func, array0, arrays... );
  }

/*
 * Serial execution overloads
 */
   template<typename   FuncObj,
            typename    ElemT0,
            typename... ElemTs,
            GridType        GT,
            ArraySizing    AS0,
            ArraySizing... ASs>
   void for_each_idx( execution::serial_policy,
                      FuncObj                     func,
                      Array1<ElemT0,GT,AS0>&    array0,
                      Array1<ElemTs,GT,ASs>&... arrays )
  {
      (( assert(   (array0.shape() == arrays.shape())
                && "par::for_each - arrays must be the same shape" ) ),... );

      const size_t n = array0.shape(0);

      for( size_t i=0; i<n; ++i )
     {
         const Idx1<GT> idx{i};
         func( idx, array0(idx),
                    arrays(idx)... );
     }
      return;
  }

   template<typename   FuncObj,
            typename    ElemT0,
            typename... ElemTs,
            GridType        GT,
            ArraySizing    AS0,
            ArraySizing... ASs>
   void for_each_idx( execution::serial_policy,
                      FuncObj                     func,
                      Array2<ElemT0,GT,AS0>&    array0,
                      Array2<ElemTs,GT,ASs>&... arrays )
  {
      (( assert(   (array0.shape() == arrays.shape())
                && "par::for_each - arrays must be the same shape" ) ),... );

      const size_t ni = array0.shape(0);
      const size_t nj = array0.shape(1);

      for( size_t i=0; i<ni; ++i )
     {
         for( size_t j=0; j<nj; ++j )
        {
            const Idx2<GT> idx{i,j};
            func( idx, array0(idx),
                       arrays(idx)... );
        }
     }
      return;
  }

   template<typename   FuncObj,
            typename    ElemT0,
            typename... ElemTs,
            GridType        GT,
            ArraySizing    AS0,
            ArraySizing... ASs>
   void for_each_idx( execution::serial_policy,
                      FuncObj                     func,
                      Array3<ElemT0,GT,AS0>&    array0,
                      Array3<ElemTs,GT,ASs>&... arrays )
  {
      (( assert(   (array0.shape() == arrays.shape())
                && "par::for_each - arrays must be the same shape" ) ),... );

      const size_t ni = array0.shape(0);
      const size_t nj = array0.shape(1);
      const size_t nk = array0.shape(2);

      for( size_t i=0; i<ni; ++i )
     {
         for( size_t j=0; j<nj; ++j )
        {
            for( size_t k=0; k<nk; ++k )
           {
               const Idx3<GT> idx{i,j,k};
               func( idx, array0(idx),
                          arrays(idx)... );
           }
        }
     }
      return;
  }


   template<typename   FuncObj,
            typename    ElemT0,
            typename... ElemTs,
            GridType        GT,
            ArraySizing    AS0,
            ArraySizing... ASs>
   void for_each_idx(       execution::serial_policy,
                            FuncObj                     func,
                      const Array1<ElemT0,GT,AS0>&    array0,
                      const Array1<ElemTs,GT,ASs>&... arrays )
  {
      (( assert(   (array0.shape() == arrays.shape())
                && "par::for_each - arrays must be the same shape" ) ),... );

      const size_t n = array0.shape(0);

      for( size_t i=0; i<n; ++i )
     {
         const Idx1<GT> idx{i};
         func( idx, array0(idx),
                    arrays(idx)... );
     }
      return;
  }

   template<typename   FuncObj,
            typename    ElemT0,
            typename... ElemTs,
            GridType        GT,
            ArraySizing    AS0,
            ArraySizing... ASs>
   void for_each_idx(       execution::serial_policy,
                            FuncObj                     func,
                      const Array2<ElemT0,GT,AS0>&    array0,
                      const Array2<ElemTs,GT,ASs>&... arrays )
  {
      (( assert(   (array0.shape() == arrays.shape())
                && "par::for_each - arrays must be the same shape" ) ),... );

      const size_t ni = array0.shape(0);
      const size_t nj = array0.shape(1);

      for( size_t i=0; i<ni; ++i )
     {
         for( size_t j=0; j<nj; ++j )
        {
            const Idx2<GT> idx{i,j};
            func( idx, array0(idx),
                       arrays(idx)... );
        }
     }
      return;
  }

   template<typename   FuncObj,
            typename    ElemT0,
            typename... ElemTs,
            GridType        GT,
            ArraySizing    AS0,
            ArraySizing... ASs>
   void for_each_idx(       execution::serial_policy,
                            FuncObj                     func,
                      const Array3<ElemT0,GT,AS0>&    array0,
                      const Array3<ElemTs,GT,ASs>&... arrays )
  {
      (( assert(   (array0.shape() == arrays.shape())
                && "par::for_each - arrays must be the same shape" ) ),... );

      const size_t ni = array0.shape(0);
      const size_t nj = array0.shape(1);
      const size_t nk = array0.shape(2);

      for( size_t i=0; i<ni; ++i )
     {
         for( size_t j=0; j<nj; ++j )
        {
            for( size_t k=0; k<nk; ++k )
           {
               const Idx3<GT> idx{i,j,k};
               func( idx, array0(idx),
                          arrays(idx)... );
           }
        }
     }
      return;
  }


/*
 * OpenMP execution overloads
 */
   template<typename   FuncObj,
            typename    ElemT0,
            typename... ElemTs,
            GridType        GT,
            ArraySizing    AS0,
            ArraySizing... ASs>
   void for_each_idx( execution::openmp_policy,
                      FuncObj                     func,
                      Array1<ElemT0,GT,AS0>&    array0,
                      Array1<ElemTs,GT,ASs>&... arrays )
  {
      (( assert(   (array0.shape() == arrays.shape())
                && "par::for_each - arrays must be the same shape" ) ),... );

      const size_t n = array0.shape(0);

# ifdef _OPENMP
   # pragma omp parallel for
# endif
      for( size_t i=0; i<n; ++i )
     {
         const Idx1<GT> idx{i};
         func( idx, array0(idx),
                    arrays(idx)... );
     }
      return;
  }

   template<typename   FuncObj,
            typename    ElemT0,
            typename... ElemTs,
            GridType        GT,
            ArraySizing    AS0,
            ArraySizing... ASs>
   void for_each_idx( execution::openmp_policy,
                      FuncObj                     func,
                      Array2<ElemT0,GT,AS0>&    array0,
                      Array2<ElemTs,GT,ASs>&... arrays )
  {
      (( assert(   (array0.shape() == arrays.shape())
                && "par::for_each - arrays must be the same shape" ) ),... );

      const size_t ni = array0.shape(0);
      const size_t nj = array0.shape(1);

# ifdef _OPENMP
   # pragma omp parallel for
# endif
      for( size_t i=0; i<ni; ++i )
     {
         for( size_t j=0; j<nj; ++j )
        {
            const Idx2<GT> idx{i,j};
            func( idx, array0(idx),
                       arrays(idx)... );
        }
     }
      return;
  }

   template<typename   FuncObj,
            typename    ElemT0,
            typename... ElemTs,
            GridType        GT,
            ArraySizing    AS0,
            ArraySizing... ASs>
   void for_each_idx( execution::openmp_policy,
                      FuncObj                     func,
                      Array3<ElemT0,GT,AS0>&    array0,
                      Array3<ElemTs,GT,ASs>&... arrays )
  {
      (( assert(   (array0.shape() == arrays.shape())
                && "par::for_each - arrays must be the same shape" ) ),... );

      const size_t ni = array0.shape(0);
      const size_t nj = array0.shape(1);
      const size_t nk = array0.shape(2);

# ifdef _OPENMP
   # pragma omp parallel for
# endif
      for( size_t i=0; i<ni; ++i )
     {
         for( size_t j=0; j<nj; ++j )
        {
            for( size_t k=0; k<nk; ++k )
           {
               const Idx3<GT> idx{i,j,k};
               func( idx, array0(idx),
                          arrays(idx)... );
           }
        }
     }
      return;
  }


   template<typename   FuncObj,
            typename    ElemT0,
            typename... ElemTs,
            GridType        GT,
            ArraySizing    AS0,
            ArraySizing... ASs>
   void for_each_idx(       execution::openmp_policy,
                            FuncObj                     func,
                      const Array1<ElemT0,GT,AS0>&    array0,
                      const Array1<ElemTs,GT,ASs>&... arrays )
  {
      (( assert(   (array0.shape() == arrays.shape())
                && "par::for_each - arrays must be the same shape" ) ),... );

      const size_t n = array0.shape(0);

# ifdef _OPENMP
   # pragma omp parallel for
# endif
      for( size_t i=0; i<n; ++i )
     {
         const Idx1<GT> idx{i};
         func( idx, array0(idx),
                    arrays(idx)... );
     }
      return;
  }

   template<typename   FuncObj,
            typename    ElemT0,
            typename... ElemTs,
            GridType        GT,
            ArraySizing    AS0,
            ArraySizing... ASs>
   void for_each_idx(       execution::openmp_policy,
                            FuncObj                     func,
                      const Array2<ElemT0,GT,AS0>&    array0,
                      const Array2<ElemTs,GT,ASs>&... arrays )
  {
      (( assert(   (array0.shape() == arrays.shape())
                && "par::for_each - arrays must be the same shape" ) ),... );

      const size_t ni = array0.shape(0);
      const size_t nj = array0.shape(1);

# ifdef _OPENMP
   # pragma omp parallel for
# endif
      for( size_t i=0; i<ni; ++i )
     {
         for( size_t j=0; j<nj; ++j )
        {
            const Idx2<GT> idx{i,j};
            func( idx, array0(idx),
                       arrays(idx)... );
        }
     }
      return;
  }

   template<typename   FuncObj,
            typename    ElemT0,
            typename... ElemTs,
            GridType        GT,
            ArraySizing    AS0,
            ArraySizing... ASs>
   void for_each_idx(       execution::openmp_policy,
                            FuncObj                     func,
                      const Array3<ElemT0,GT,AS0>&    array0,
                      const Array3<ElemTs,GT,ASs>&... arrays )
  {
      (( assert(   (array0.shape() == arrays.shape())
                && "par::for_each - arrays must be the same shape" ) ),... );

      const size_t ni = array0.shape(0);
      const size_t nj = array0.shape(1);
      const size_t nk = array0.shape(2);

# ifdef _OPENMP
   # pragma omp parallel for
# endif
      for( size_t i=0; i<ni; ++i )
     {
         for( size_t j=0; j<nj; ++j )
        {
            for( size_t k=0; k<nk; ++k )
           {
               const Idx3<GT> idx{i,j,k};
               func( idx, array0(idx),
                          arrays(idx)... );
           }
        }
     }
      return;
  }


/*
 * ------------------------- par::transform ------------------------
 */

/*
 * Applies the given functor to each element pack in an array pack, and stores the result in another array
 */
   template<typename   FuncObj,
            int           NDIM,
            typename    ElemTd,
            typename    ElemT0,
            typename... ElemTs,
            GridType        GT,
            ArraySizing    ASd,
            ArraySizing    AS0,
            ArraySizing... ASs>
   void transform( const FuncObj&                      func,
                         Array<ElemTd,NDIM,GT,ASd>&     dst,
                   const Array<ElemT0,NDIM,GT,AS0>&    src0,
                   const Array<ElemTs,NDIM,GT,ASs>&... srcs )
  {
      transform( execution::seq, func, dst, src0, srcs... );
      return;
  }

/*
 * Return value overload
 */
   template<typename   FuncObj,
            int           NDIM,
            typename    ElemT0,
            typename... ElemTs,
            GridType        GT,
            ArraySizing    AS0,
            ArraySizing... ASs>
   Array<ElemT0,NDIM,GT,AS0> transform( const FuncObj&                      func,
                                        const Array<ElemT0,NDIM,GT,AS0>&    src0,
                                        const Array<ElemTs,NDIM,GT,ASs>&... srcs )
  {
      return transform( execution::seq, func, src0, srcs... );
  }

   template<execution_policy Policy,
            typename        FuncObj,
            int                NDIM,
            typename         ElemT0,
            typename...      ElemTs,
            GridType             GT,
            ArraySizing         AS0,
            ArraySizing...      ASs>
   Array<ElemT0,NDIM,GT,AS0> transform(       Policy                      policy,
                                        const FuncObj&                      func,
                                        const Array<ElemT0,NDIM,GT,AS0>&    src0,
                                        const Array<ElemTs,NDIM,GT,ASs>&... srcs )
  {
      Array<ElemT0,NDIM,GT,AS0> dst(src0.shape());
      transform( policy, func, dst, src0, srcs... );
      return dst;
  }

/*
 * Serial execution
 */
   template<typename   FuncObj,
            int           NDIM,
            typename    ElemTd,
            typename    ElemT0,
            typename... ElemTs,
            GridType        GT,
            ArraySizing    ASd,
            ArraySizing    AS0,
            ArraySizing... ASs>
   void transform(       execution::serial_policy,
                   const FuncObj&                      func,
                         Array<ElemTd,NDIM,GT,ASd>&     dst,
                   const Array<ElemT0,NDIM,GT,AS0>&    src0,
                   const Array<ElemTs,NDIM,GT,ASs>&... srcs )
  {
         assert(   (dst.shape() == src0.shape())
                && "par::transform - arrays must be the same shape" );

      (( assert(   (dst.shape() == srcs.shape())
                && "par::transform - arrays must be the same shape" ) ),... );

      const size_t len = dst.flattened_length();
   
      for( size_t i=0; i<len; ++i )
     {
         dst.flatten(i) = func( src0.flatten(i),
                                srcs.flatten(i)... );
     }
      return;
  }

/*
 * OpenMP execution
 */
   template<typename   FuncObj,
            int           NDIM,
            typename    ElemTd,
            typename    ElemT0,
            typename... ElemTs,
            GridType        GT,
            ArraySizing    ASd,
            ArraySizing    AS0,
            ArraySizing... ASs>
   void transform(       execution::openmp_policy,
                   const FuncObj&                      func,
                         Array<ElemTd,NDIM,GT,ASd>&     dst,
                   const Array<ElemT0,NDIM,GT,AS0>&    src0,
                   const Array<ElemTs,NDIM,GT,ASs>&... srcs )
  {
         assert(   (dst.shape() == src0.shape())
                && "par::transform - arrays must be the same shape" );

      (( assert(   (dst.shape() == srcs.shape())
                && "par::transform - arrays must be the same shape" ) ),... );

      const size_t len = dst.flattened_length();
   
# ifdef _OPENMP
   # pragma omp parallel for
# endif
      for( size_t i=0; i<len; ++i )
     {
         dst.flatten(i) = func( src0.flatten(i),
                                srcs.flatten(i)... );
     }
      return;
  }

/*
 * ------------------------- par::transform_reduce ------------------------
 */

/*
 * Applies the given functor to each element pack in an array pack, and reduces result with given binary functor
 */
   template<int                  NDIM,
            typename TransformFuncObj,
            typename    ReduceFuncObj,
            typename    ReductionType,
            typename           ElemT0,
            typename...        ElemTs,
            GridType               GT,
            ArraySizing           AS0,
            ArraySizing...        ASs>
   ReductionType transform_reduce( const TransformFuncObj&            tfunc,
                                   const    ReduceFuncObj&            rfunc,
                                         ReductionType                 init,
                                   const Array<ElemT0,NDIM,GT,AS0>&    src0,
                                   const Array<ElemTs,NDIM,GT,ASs>&... srcs )
  {
      return transform_reduce( execution::seq, tfunc, rfunc, init, src0, srcs... );
  }

/*
 * Serial execution
 */
   template<int                  NDIM,
            typename TransformFuncObj,
            typename    ReduceFuncObj,
            typename    ReductionType,
            typename           ElemT0,
            typename...        ElemTs,
            GridType               GT,
            ArraySizing           AS0,
            ArraySizing...        ASs>
   ReductionType transform_reduce(       execution::serial_policy,
                                   const TransformFuncObj&            tfunc,
                                   const    ReduceFuncObj&            rfunc,
                                         ReductionType                 init,
                                   const Array<ElemT0,NDIM,GT,AS0>&    src0,
                                   const Array<ElemTs,NDIM,GT,ASs>&... srcs )
  {
      (( assert(   (src0.shape() == srcs.shape())
                && "par::transform_reduce - arrays must be the same shape" ) ),... );

      const size_t len = src0.flattened_length();

      for( size_t i=0; i<len; ++i )
     {
         init = rfunc( std::move(init),
                       tfunc( src0.flatten(i),
                              srcs.flatten(i)... ) );
     }
      return init;
  }

/*
 * OpenMP execution
 */
   template<int                  NDIM,
            typename TransformFuncObj,
            typename    ReduceFuncObj,
            typename    ReductionType,
            typename           ElemT0,
            typename...        ElemTs,
            GridType               GT,
            ArraySizing           AS0,
            ArraySizing...        ASs>
   ReductionType transform_reduce(       execution::openmp_policy,
                                   const TransformFuncObj&            tfunc,
                                   const    ReduceFuncObj&            rfunc,
                                         ReductionType                 init,
                                   const Array<ElemT0,NDIM,GT,AS0>&    src0,
                                   const Array<ElemTs,NDIM,GT,ASs>&... srcs )
  {
      (( assert(   (src0.shape() == srcs.shape())
                && "par::transform_reduce - arrays must be the same shape" ) ),... );

      const size_t len = src0.flattened_length();

# ifdef _OPENMP
      const int nthreads = omp_get_max_threads();
# else
      const int nthreads = 1;
# endif

   // each thread reduces into thread_local_init[thread_num]
      std::vector<ReductionType> thread_local_init(nthreads,init);

   // reduce values for each threads share of elements
# ifdef _OPENMP
   # pragma omp parallel for
# endif
      for( size_t i=0; i<len; ++i )
     {
         const int thread_num = omp_get_thread_num();

         thread_local_init[thread_num] = rfunc( std::move(thread_local_init[thread_num]),
                                                tfunc( src0.flatten(i),
                                                       srcs.flatten(i)... ) );
     }

   // reduce values from each thread
      for( const ReductionType& tli : thread_local_init )
     {
         init = rfunc( std::move(init), tli );
     }

      return init;
  }
}
