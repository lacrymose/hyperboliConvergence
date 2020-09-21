
# include <utils/concepts.h>

# include <array>
//# include <tuple>
# include <iostream>
# include <type_traits>

/*
   template<size_t M, typename... Args>
   struct HoldTuple
  {
      std::tuple<Args...> args;

      HoldTuple( const Args&... as ) : args(as...){};

      template<size_t N>
         requires N<sizeof...(Args)
      auto& arg(){ return std::get<N>(args); }
  };

   template<size_t M, typename... Args>
   auto make_HoldTuple( Args&&... args )
  {
      return HoldTuple<M,Args...>{args...};
  }
*/

   template<floating_point Real, int N>
   using my_array = std::array<Real,N>;

   template<int N, floating_point Real>
   Real sum( my_array<Real,N> arr )
  {
      Real x{0};
      for( const Real a : arr ){ x+=a; }
      return x;
  }

   int main()
  {
      const my_array<float,3> x{1.,2.,3.};
      const float s = sum(x);

      std::cout << sizeof(int) << std::endl;
      std::cout << sizeof(unsigned int) << std::endl;
      std::cout << sizeof(long unsigned int) << std::endl;

      static_assert( std::is_same_v<long unsigned int, size_t> );
      return 0;
  }
