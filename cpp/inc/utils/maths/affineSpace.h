
# pragma once

# include <utils/concepts.h>

# include <array>

# include <type_traits>

/*
 * This header defines two base types corresponding to a point and a displacement (delta) in an affine space.
 *
 *    These types must both be inherited from using CRTP. Instances of the two derived types will act as members of an affine space, with only the expected arithmetic operations allowed,
 *    i.e. displacements can be added/subtracted to give other displacements,
 *         displacements can be added/subtracted to/from points to give other points,
 *         points can be subtracted to give displacements,
 *         displacements can be multiplied by scalars to give other displacements,
 *         points cannot be added,nor multiplied by scalars
 *         in-place versions of these operations are also defined
 *    No assumption is made on whether the space is linear or nonlinear.
 *    All instances are assumed to be represented in the standard basis. Transformations between bases must be implemented by the deriving classes.
 *    Individual elements of Points/Deltas can be accessed using the array accessor []
 *
 *    Code example:
 *
 *       template<int N>
 *       struct CartesianPoint : AffinePointBase<N,CartesianPoint<N>,CartesianDelta<N>,double> {};
 *
 *       template<int N>
 *       struct CartesianDelta : AffineDeltaBase<N,CartesianPoint<N>,CartesianDelta<N>,double> {};
 *
 *       CartesianPoint<3> x0,x1;
 *       CartesianDelta<3> d0,d1;
 *       const double a = ...
 *
 *       // ... initialise ... x0[0]= ...
 *
 *       d1 = d0+d1;
 *       d1 = d0-d1;
 *
 *       x1 = x0+d0;
 *       x1 = x0-d0;
 *
 *       d1 = x1-x0:
 *
 *       d1 = a*d0:
 *       d1 = d0*a:
 *       d1 = d0/a:
 *
 *       // x0 = x0+x1;    // will not compile!
 *       // x0 =  a*x1;    // will not compile!
 *
 *    By default, there is no way to convert between a point and a delta, except by manually copying each element over (x0[0]=d0[0]...)
 *    To enable explicit construction of Points/Deltas from each other, the derived types must call "using" on the base constructors:
 *       Explicit conversion is then allowed, but any attempt to implicitly convert will still result in a compiler error
 *       Note that strictly one-way explicit conversion is possible. eg if "using ..." is only in the Point derived type, it will only be possible to construct points from deltas, not the reverse
 *
 *    Code example:
 *
 *       template<int N>
 *       struct CartesianPoint : AffinePointBase<N,CartesianPoint<N>,CartesianDelta<N>>
 *      {
 *          using AffinePointBase<N,CartesianPoint<N>,CartesianDelta<N>>::AffinePointBase;
 *      };
 *
 *       template<int N>
 *       struct CartesianDelta : AffineDeltaBase<N,CartesianPoint<N>,CartesianDelta<N>>
 *      {
 *          using AffineDeltaBase<N,CartesianPoint<N>,CartesianDelta<N>>::AffineDeltaBase;
 *      };
 *
 *       CartesianPoint<3> p = ... initialise ...
 *       // CartesianDelta<3> d = p; // always compiler error
 *       CartesianDelta<3> d = CartesianDelta<3>(p); // allowed only if "using" base class constructor
 *
 *    Brace initialisation is possible, but requires an extra pair of brackets if not "using" the base class constructors, because c++
 *       CartesianPoint<3> p{{0.,1.,2.}};   // allowed only if NOT "using" base class constructor
 *       CartesianPoint<3> p{ 0.,1.,2. };   // allowed only if     "using" base class constructor
 *
 *    Sending either Point/Delta to stream sends each element in turn separated by a space
 *    ostream << x0  // -> // ostream << x0[0] << " " << x0[1] ... << x0[n];
 */

// --------------- forward declarations ---------------

   template<int NDIM, typename Point, typename Delta, floating_point Real> struct AffinePointBase;
   template<int NDIM, typename Point, typename Delta, floating_point Real> struct AffineDeltaBase;


// --------------- type traits ---------------
/*
 * Note that these type-traits are implementation details used to define the necessary operations specifically for Affine(Point/Delta)Base classes
 * They are not meant to be used as concepts to determine if any generic types behave as an affine space
*/

/*
 * returns true if argument is derived from affine point CRTP base
 */
   template<typename T>
   struct has_affinepoint_base : std::false_type {};

   template<typename T>
      requires   requires(){ { T::N } -> int; }
              && requires{ typename T::value_type; }
              && requires{ typename T::point_type; }
              && requires{ typename T::delta_type; }
              && std::is_same_v<T,typename T::point_type>
   struct has_affinepoint_base<T>
            : std::is_base_of<AffinePointBase<T::N,
                                              typename T::point_type,
                                              typename T::delta_type,
                                              typename T::value_type>,T> {};

/*
 * returns true if argument is derived from affine delta CRTP base
 */
   template<typename T>
   struct has_affinedelta_base : std::false_type {};

   template<typename T>
      requires   requires(){ { T::N } -> int; }
              && requires{ typename T::value_type; }
              && requires{ typename T::point_type; }
              && requires{ typename T::delta_type; }
              && std::is_same_v<T,typename T::delta_type>
   struct has_affinedelta_base<T>
            : std::is_base_of<AffineDeltaBase<T::N,
                                              typename T::point_type,
                                              typename T::delta_type,
                                              typename T::value_type>,T> {};

/*
 * returns true if arguments are derived from affine point/delta CRTP bases, have equal dimension, and refer to each other as corresponding point/delta
 */
   template<typename P, typename D>
   struct is_affine_pair : std::false_type {};

   template<typename P, typename D>
      requires   has_affinepoint_base<P>::value
              && has_affinedelta_base<D>::value
              && ( P::N == D::N )
              && std::is_same_v<typename P::value_type,
                                typename D::value_type>
              && std::is_same_v<P,typename D::point_type>
              && std::is_same_v<D,typename P::delta_type>
   struct is_affine_pair<P,D> : std::true_type {};


// --------------- concepts ---------------

/*
 * pairs of types satisfying this concept behave as an affine space point/displacement pairing
 * i.e they have the necessary arithmetic operation defined, but do not define arithmetic operations which are not allowed in an affine space, such as addition of points
 */
//   template<typename P, typename D>
//   concept bool AffineSpace = 


// --------------- affine space point/displacement ---------------

/*
 * base CRTP type for a Point in an affine space
 */
   template<int NDIM, typename Point, typename Delta, floating_point Real>
   struct AffinePointBase
  {
      using value_type = Real;
      using point_type = Point;
      using delta_type = Delta;
      constexpr static int N=NDIM;
      std::array<Real,N> elems;

   // default, copy and move constructors
      AffinePointBase() = default;
      AffinePointBase( const AffinePointBase&  ) = default;
      AffinePointBase(       AffinePointBase&& ) = default;

   // only allow explicit copy/move construction from Delta
      // explicit Delta constructors are only enabled by adding following to derived Point definition:
      // using AffinePointBase<NDIM,Point,Delta>::AffinePointBase
      explicit AffinePointBase( const Delta&  d ) noexcept : elems(d.elems) {}
      explicit AffinePointBase(       Delta&& d ) noexcept : elems(std::move(d.elems)) {}

   // initialiser list constructor, must be of length NDIM list of type Real
      template<typename... Ts>
      requires   (sizeof...(Ts)==NDIM)
              && (std::is_same_v<Real,Ts> && ...)
      AffinePointBase( Ts... t ) noexcept : elems{t...} {}

   // copy/move assignment
      AffinePointBase& operator=( const AffinePointBase&  ) = default;
      AffinePointBase& operator=(       AffinePointBase&& ) = default;

   // accessors
            Real& operator[]( const int i )       { return elems[i]; }
      const Real& operator[]( const int i ) const { return elems[i]; }

   // in-place arithmetic
      Point& operator+=( const Delta& d );
      Point& operator-=( const Delta& d );
  };


/*
 * base CRTP type for a displacement in an affine space
 */
   template<int NDIM, typename Point, typename Delta, floating_point Real>
   struct AffineDeltaBase
  {
      using value_type = Real;
      using point_type = Point;
      using delta_type = Delta;
      constexpr static int N=NDIM;
      std::array<Real,N> elems;

   // default, copy and move constructors
      AffineDeltaBase() = default;
      AffineDeltaBase( const AffineDeltaBase&  ) = default;
      AffineDeltaBase(       AffineDeltaBase&& ) = default;

   // only allow explicit copy/move constructor from Point
      // explicit Point constructors are only enabled by adding following to derived Delta definition:
      // using AffineDeltaBase<NDIM,Point,Delta>::AffineDeltaBase
      explicit AffineDeltaBase( const Point&  p ) noexcept : elems(p.elems) {}
      explicit AffineDeltaBase(       Point&& p ) noexcept : elems(std::move(p.elems)) {}

   // initialiser list constructor, must be of length NDIM list of type Real
      template<typename... Ts>
      requires   (sizeof...(Ts)==NDIM)
              && (std::is_same_v<Real,Ts> && ...)
      AffineDeltaBase( Ts... t ) noexcept : elems{t...} {}

   // copy/move assignment
      AffineDeltaBase& operator=( const AffineDeltaBase&  ) = default;
      AffineDeltaBase& operator=(       AffineDeltaBase&& ) = default;

   // accessors
            Real& operator[]( const int i )       { return elems[i]; }
      const Real& operator[]( const int i ) const { return elems[i]; }

   // in-place arithmetic
      Delta& operator+=( const Delta& d );
      Delta& operator-=( const Delta& d );
      Delta& operator*=( const Real a );
      Delta& operator/=( const Real a );
  };

# include <utils/maths/affineSpace.ipp>

