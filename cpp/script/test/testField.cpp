
# include <solutionField/solutionField.h>

# include <vector>
# include <array>

   struct NonCopyable
  {
      NonCopyable() = default;
      NonCopyable( const NonCopyable&  ) = delete;
      NonCopyable(       NonCopyable&& ) = default;

      NonCopyable( int i ){};

      NonCopyable& operator=( const NonCopyable&  ) = delete;
      NonCopyable& operator=(       NonCopyable&& ) = default;

     ~NonCopyable() = default;
  };

   int main()
  {
      std::array<NonCopyable,2> arr{0,1};

//    std::vector<NonCopyable> vec(2,2);
      
      return 0;
  }
