
# include <parallalg/algorithm.h>

# include <parallalg/array.h>

# include <type_traits>

# include <iostream>
# include <iomanip>

/*
 * non-trivially copyable type
 */
   struct NoCopy
  {
      int i;

      NoCopy() = default;
      NoCopy( const int j ){ i=j; }
      NoCopy( const NoCopy& other ){ i=other.i; }
  };

   std::ostream& operator<<( std::ostream& os, const NoCopy& nc )
  {
      os << nc.i;
      return os;
  }

   template<typename ElemT, par::ArraySizing AS>
   void print_array( par::Array<ElemT,2,AS>& array )
  {
      for( size_t i=0; i<array.shape(0); i++ )
     {
         for( size_t j=0; j<array.shape(1); j++ )
        {
            std::cout << std::setw(4) << array[{i,j}] << " ";
        }
         std::cout << std::endl;
     }
      std::cout << std::endl;
  }

   int main()
  {
      par::Shape<2> shape34{3,4};

   // par::fill and par::copy for not trivially copyable type

      std::cout << "par::fill and par::copy for not trivially copyable type" << std::endl;
      std::cout << std::endl;

      std::cout << "true:  " << true  << std::endl;
      std::cout << "false: " << false << std::endl;
      std::cout << "NoCopy is trivially copyable? " << std::is_trivially_copyable_v<NoCopy> << std::endl;
      std::cout << std::endl;

      par::Array<NoCopy,2> arrayNC0(shape34);
      
      print_array( arrayNC0 );
      par::fill( arrayNC0, NoCopy(4) );
      print_array( arrayNC0 );

      par::Array<NoCopy,2> arrayNC1(shape34);

      print_array( arrayNC1 );
      par::copy(   arrayNC1, arrayNC0 );
      print_array( arrayNC1 );

   // par::fill and par::copy for trivially copyable type

      std::cout << "par::fill and par::copy for trivially copyable type" << std::endl;
      std::cout << std::endl;

      par::Array<int,2> array0(shape34);

      print_array( array0 );
      par::fill(   array0, 2 );
      print_array( array0 );

      par::Array<int,2> array1(shape34);

      print_array( array1 );
      par::copy(   array1, array0 );
      print_array( array1 );

   // par::for_each

      std::cout << "par::for_each" << std::endl;
      std::cout << std::endl;


      for( size_t i=0; i<array0.shape(0); i++ )
     {
         for( size_t j=0; j<array0.shape(1); j++ )
        {
            array0[{i,j}] = par::Idx<2>{i,j} * array0.stride();
        }
     }
      print_array( array0 );

      par::for_each( []( const int i ) -> void
                       { std::cout << i << std::endl; },
                     array0 );

      std::cout << std::endl;

      par::for_each( []( const int i, const int j ) -> void
                       { std::cout << i << " " << j << std::endl; },
                     array0,
                     array1 );

      std::cout << std::endl;

   // par::transform

      std::cout << "par::transform" << std::endl;
      std::cout << std::endl;


      par::Array<int,2> array2(shape34);

      print_array( array1 );
      par::transform( []( const int i ) -> int
                        { return i-6; },
                      array1,
                      array0 );
      print_array( array1 );

      par::transform( []( const int i, const int j ) -> int
                        { return i-j; },
                      array2,
                      array0,array1 );
      print_array( array2 );

   // par::transform_reduce

      std::cout << "par::transform_reduce" << std::endl;
      std::cout << std::endl;

   // array2 = 6 for all elements
   // sum3 = flat_length*2 = 12*2 = 24
      const int sum3 = par::transform_reduce(
                                           // transform to a third
                                              []( const int i ) -> int
                                                { return i/3; },
                                           // reduce the sum
                                              []( const int i, const int j ) -> int
                                                { return i+j; },
                                           // start sum from 0
                                              0,
                                              array2
                                            );
      std::cout << sum3 << std::endl;

   // array1 = array0-6
   // array1+array0+6 = 0,1,2,...
   // sum is arithmetic sum, 0+0 + 1+1 = len*(len-1) = 132
      const int sumg = par::transform_reduce(
                                           // transform to a third
                                              []( const int i, const int j ) -> int
                                                { return i+j+6; },
                                           // reduce the sum
                                              []( const int i, const int j ) -> int
                                                { return i+j; },
                                           // start sum from 0
                                              0,
                                              array0,
                                              array1
                                            );

      std::cout << sumg << std::endl;

      return 0;
  }
