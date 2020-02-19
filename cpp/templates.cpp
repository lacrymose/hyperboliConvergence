
# include <iostream>

   template< char C >
   struct CharTemplate
  {
      float    f0;
      float    f1;

      CharTemplate<C>();
      CharTemplate<C>( CharTemplate<C>& source );

      template< char D >
      CharTemplate<C>( CharTemplate<D>& source );

      template< char D >
      void add( CharTemplate<D>& source );

      void cout(){ std::cout << C << std::endl; }

      void cout2( CharTemplate<C> second ){ std::cout << C << " " << C << std::endl; }
  };

// simple constructor definitions
   template< char C >
   CharTemplate<C>::CharTemplate(){ f0=0.; f1=1.; }

   template<>
   CharTemplate<'b'>::CharTemplate(){ f0=2.; f1=3.; }

// single template constructor
   template< char C >
   CharTemplate<C>::CharTemplate( CharTemplate<C>& source ){ f0=source.f0; f1=source.f1; }

// double template constructor
   template< char C >
   template< char D >
   CharTemplate<C>::CharTemplate( CharTemplate<D>& source ){ f0=source.f0; f1=source.f1; }

   template<>
   template<>
   CharTemplate<'a'>::CharTemplate( CharTemplate<'b'>& source ){ f0=10.; f1=11; }

// method specializations
   template<>
   void CharTemplate<'a'>::cout(){ std::cout << "I am an a!" << std::endl; }

// method with templated arguments
   template< char C >
   template< char D >
   void CharTemplate<C>::add( CharTemplate<D>& source )
  {
      std::cout << f0+source.f0 << std::endl;
      std::cout << f1+source.f1 << std::endl;
  }

   template<>
   template<>
   void CharTemplate<'b'>::add( CharTemplate<'a'>& source )
  {
      std::cout << -(f0+source.f0) << std::endl;
      std::cout << -(f1+source.f1) << std::endl;
  }

   int main()
  {
      CharTemplate<'a'> aTemplate;
      CharTemplate<'b'> bTemplate;

   // standard constructor for <'a'>
      std::cout << aTemplate.f0 << std::endl;
      std::cout << aTemplate.f1 << std::endl;

   // specialized constructor for <'b'>
      std::cout << bTemplate.f0 << std::endl;
      std::cout << bTemplate.f1 << std::endl;

   // specialized cout for a
      aTemplate.cout();
   // standard cout for b
      bTemplate.cout();

   // single templated cout2
      aTemplate.cout2( aTemplate );
//    aTemplate.cout2( bTemplate );   // won't compile because template parameters do not match

   // single template constructor
      aTemplate.f0=5.;
      aTemplate.f1=6.;
      CharTemplate<'a'> a2Template( aTemplate );
      std::cout << a2Template.f0 << std::endl;
      std::cout << a2Template.f1 << std::endl;

   // double template constructor
      aTemplate.f0=7.;
      aTemplate.f1=8.;
      CharTemplate<'b'> b2Template( aTemplate );
      std::cout << b2Template.f0 << std::endl;
      std::cout << b2Template.f1 << std::endl;

      CharTemplate<'a'> a3Template( bTemplate );
      std::cout << a3Template.f0 << std::endl;
      std::cout << a3Template.f1 << std::endl;

   // double template method
      aTemplate.add( bTemplate );
      bTemplate.add( aTemplate );

      return 0;
  }
