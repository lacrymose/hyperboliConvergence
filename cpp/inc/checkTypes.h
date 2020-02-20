
   template<typename A, typename B>
   struct CheckTypes
  {
      static const bool val=false;
  };

   template<typename A>
   struct CheckTypes<A,A>
  {
      static const bool val=true;
  };


