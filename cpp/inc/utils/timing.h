
# pragma once

# include <chrono>
# include <string>

# include <iostream>
# include <iomanip>

# include <cassert>

namespace utils
{
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
   template<typename Clock     = std::chrono::steady_clock,
            typename TimeUnits = std::chrono::milliseconds>
      requires Clock_and_TimeUnits<Clock,TimeUnits>
   class LifetimeTimer
  {
   private:

      using time_point = typename Clock::time_point;

      std::ostream&          out;
      const std::string  message;
      const time_point     start;

   public:

   // record time at construction
      LifetimeTimer( const std::string& m )
                    : out(std::clog), message(m), start(Clock::now()) {}

      LifetimeTimer( const std::ostream& o, const std::string& m )
                    : out(o), message(m), start(Clock::now()) {}

   // calculate and print time between construction and destruction
     ~LifetimeTimer()
     {
         const time_point end = Clock::now();
         const auto timing = std::chrono::duration_cast<TimeUnits>(end-start).count();
         out << message << timing << units_string( TimeUnits{} ) << "\n";
     }
  };

/*
 * class which accumulates the timing of an intermittent process by being started/stopped.
 *    uses Clock to time (eg std::chrono::steady_clock) and measures in TimingDuration (eg std::chrono::seconds)
 *    records timing in RecordingDuration (can be different from TimingDuration)
 *    constructed with an output stream, and a message to print to the stream before printing the lifetime duration
 */
   template<typename Clock             = std::chrono::steady_clock,
            typename TimingDuration    = std::chrono::microseconds,
            typename RecordingDuration = std::chrono::microseconds>
      requires   Clock_and_TimeUnits<Clock,TimingDuration>
              && Clock_and_TimeUnits<Clock,RecordingDuration>
   class StopWatchTimer
  {
   private:

      using time_point = typename Clock::time_point;

      std::ostream&              out;
      const std::string      message;
      time_point        resume_point;
      TimingDuration   clock_reading;
      bool          is_clock_running;

   public:

      StopWatchTimer( const std::string& m )
                    : out(std::clog), message(m) { reset(); }

      StopWatchTimer( const std::ostream& o, const std::string& m )
                    : out(o), message(m) { reset(); }

     ~StopWatchTimer()
     {
         if( is_clock_running ){ pause(); }
         record_time();
     }

   // is clock currently recording?
      bool is_running(){ return is_clock_running; }

   // return current duration reading
      size_t current_timing(){ return std::chrono::duration_cast<RecordingDuration>(clock_reading).count(); }

      void reset()
     {
         clock_reading=TimingDuration{0};
         is_clock_running=false;
     }

   // start/resume timing
      void start()
     {
         assert( !is_clock_running );
         is_clock_running = true;
         resume_point = Clock::now();
     }

   // pause timing
      void pause()
     {
         const time_point pause_point = Clock::now();
         assert( is_clock_running );
         is_clock_running = false;
         clock_reading += std::chrono::duration_cast<TimingDuration>(pause_point-resume_point);
     }

   // send current accumulated time to ostream (reading is not valid if clock is running)
      void record_time()
     {
         assert( !is_clock_running );
         const auto duration_time = current_timing();
         out << message << std::setw(8) << duration_time << units_string( RecordingDuration{} ) << "\n";
     }
  };
}
