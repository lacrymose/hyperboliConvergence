
# pragma once

# include <utils/concepts.h>

# include <array>

# include <type_traits>

/*
 * This header defines a base type corresponding to a member of a vector space
 *
 *    This type must both be inherited from using CRTP. Instances of the derived type will act as members of a vector space, with only the expected arithmetic operations allowed,
 *    i.e. vectors can be added/subtracted to give other vectors,
 *         vectors can be multiplied by scalars to give other vectors,
 *         in-place versions of these operations are also defined
 *    No assumption is made on whether the space is linear or nonlinear.
 *    All instances are assumed to be represented in the standard basis. Transformations between bases must be implemented by the deriving classes.
 *    Individual elements of Vectors can be accessed using the array accessor []
 *
 *    Code example:
 *
 *       template<int N>
 *       struct CartesianVector : VectorSpaceBase<N,CartesianVector<N>> {};
 *
 *       CartesianVector<3> v0,v1;
 *       const double a = ...
 *
 *       // ... initialise ... v0[0]= ...
 *
 *       v1 = v0+v1;
 *       v1 = v0-v1;
 *
 *       v1 = a*v0:
 *       v1 = v0*a:
 *       v1 = v0/a:
 *
 *    Brace initialisation is also possible
 *       CartesianVector<3> p{0.,1.,2.};
 *
 *    Sending Vector to stream sends each element in turn separated by a space
 *    ostream << v0  // -> // ostream << v0[0] << " " << v0[1] ... << v0[n];
 */

// --------------- forward declarations ---------------

   template<int NDIM, typename Derived, floating_point Real> struct VectorSpaceBase;


// --------------- type traits ---------------
/*
 * Note that this type-trait is an implementation detail used to define the necessary operations specifically for the VectorSpaceBase class
 * It is not meant to be used as a concept to determine if any generic type behaves as a vector space
*/

/*
 * returns true if argument is derived from vector space CRTP base
 */
   template<typename T>
   struct has_vectorspace_base : std::false_type {};

   template<typename T>
      requires   requires(){ { T::N } -> int }
              && requires{ typename T::value_type; }
   struct has_vectorspace_base<T> : std::is_base_of<VectorSpaceBase<T::N,T,typename T::value_type>,T> {};


// --------------- concepts ---------------

/*
 * types satisfying this concept behave as a vector space
 * i.e they have the necessary arithmetic operation defined
 */
//   template<typename V>
//   concept bool VectorSpace =


// --------------- vector space element ---------------

/*
 * base CRTP type for an element of a vector space
 */
   template<int NDIM, typename Derived, floating_point Real>
   struct VectorSpaceBase
  {
      using value_type = Real;
      constexpr static int N=NDIM;
      std::array<Real,N> v;

   // accessors
            Real& operator[]( const int i )       { return v[i]; }
      const Real& operator[]( const int i ) const { return v[i]; }

   // in-place arithmetic
      Derived& operator+=( const Derived& other );
      Derived& operator-=( const Derived& other );
      Derived& operator*=( const Real a );
      Derived& operator/=( const Real a );
  };

# include <utils/maths/vectorSpace.ipp>

