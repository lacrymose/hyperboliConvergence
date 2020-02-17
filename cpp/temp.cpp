
# include <iostream>

   template< char C>
   struct CharTemplate
  {
      void cout(){ std::cout << C << std::endl; }

      void two( CharTemplate<C> second ){ std::cout << C << C << std::endl; }
  };

   template<>
   void CharTemplate<'a'>::cout(){ std::cout << "I am an a!" << std::endl; }

   int main()
  {
      CharTemplate<'a'> aTemplate;
      CharTemplate<'b'> bTemplate;

      aTemplate.cout();
      bTemplate.cout();

      aTemplate.two( aTemplate );
//    aTemplate.two( bTemplate );   // won't compile because template parameters do not match

      return 0;
  }
