
# pragma once

# include <tuple>

# include <cassert>
# include <utils/type-traits.h>

namespace utils
{

/*
 * recursively loop through tuple and call DoSomething with tuple element matching thingToMatch and args
      must be called with N=0 first
 */
   template<size_t                   N,
            typename      ThingToMatch,
            typename...     TupleElems,
            typename  GetThingFromElem,
            typename     CompareThings,
            typename       DoSomething,
            typename...           Args>
   void selectFromTuple( const ThingToMatch         thingToMatch,
                         const std::tuple<TupleElems...>     tpl,
                         const GetThingFromElem&       get_thing,
                         const CompareThings&     compare_things,
                         const DoSomething&         do_something,
                               Args&&...                    args )
  {
      if( compare_things( thingToMatch,
                          get_thing( std::get<N>(tpl) ) ) )
     {
         do_something( std::get<N>(tpl), std::forward<Args>(args)... );
         return;
     }
      if constexpr( N+1 < sizeof...(TupleElems) )
     {
         selectFromTuple<N+1>( thingToMatch, tpl,
                               get_thing,
                               compare_things,
                               do_something,
                               std::forward<Args>(args)... );
         return;
     }
      else
     {
         assert( false && "selectFromTuple did not match any member of tuple" );
     }
  }

   template<typename INT>
      requires is_integer_v<INT>
   constexpr INT triangular_number( INT n ){ return n*(n+1)/2; }

}
