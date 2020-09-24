
# pragma once

# include <string>
# include <iostream>
# include <chrono>
# include <tuple>

# include <cmath>
# include <cassert>

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

/*
 * returns string with unit abbreviation of std::chrono
 */
   template<typename   Int,
            typename Ratio>
   constexpr std::string units_string( const std::chrono::duration<Int,Ratio>& )
  {
      using TimeUnits = std::chrono::duration<Int,Ratio>;
      if constexpr( std::is_same_v<TimeUnits,std::chrono::hours> )
     {
         return "h";
     }
      else if constexpr( std::is_same_v<TimeUnits,std::chrono::minutes> )
     {
         return "m";
     }
      else if constexpr( std::is_same_v<TimeUnits,std::chrono::seconds> )
     {
         return "s";
     }
      else if constexpr( std::is_same_v<TimeUnits,std::chrono::milliseconds> )
     {
         return "ms";
     }
      else if constexpr( std::is_same_v<TimeUnits,std::chrono::microseconds> )
     {
         return "us";
     }
      else if constexpr( std::is_same_v<TimeUnits,std::chrono::nanoseconds> )
     {
         return "ns";
     }
      else
     {
      // will always fail
         static_assert( std::is_same_v<TimeUnits,int>, "unit string not defined" );
     }
  }

/*
 * Tests if arguments form a suitable pair for measuring timings with the Clock in units of TimeUnits
 */
   template<typename Clock,
            typename TimeUnits>
   concept bool Clock_and_TimeUnits =
         requires{ typename Clock::time_point; }
      && requires()
                {
                   { Clock::now() } -> typename Clock::time_point;
                }
      && requires( typename Clock::time_point s,
                   typename Clock::time_point e )
                {
                   std::chrono::duration_cast<TimeUnits>(e-s);
                };

/*
 * class which times its own lifetime, and prints it out at destruction
 *    uses Clock to time (eg std::chrono::steady_clock) and measures in TimeUnits (eg std::chrono::seconds)
 *    constructed with an output stream, and a message to print to the stream before printing the lifetime duration
 */
   template<typename Clock     = std::chrono::high_resolution_clock,
            typename TimeUnits = std::chrono::milliseconds>
      requires Clock_and_TimeUnits<Clock,TimeUnits>
   class Timer
  {
   private:

      using time_point = typename Clock::time_point;

      std::ostream&          out;
      const std::string  message;
      const time_point     start;

   public:

      Timer( std::ostream& o, const std::string& m )
             : out(o), message(m), start(Clock::now()) {}

     ~Timer()
     {
         const time_point end = Clock::now();
         const auto timing = std::chrono::duration_cast<TimeUnits>(end-start).count();
         out << message << timing << units_string( TimeUnits{} ) << "\n";
     }
  };
}

