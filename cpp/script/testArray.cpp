
# include <mdarray/mdarray.h>

# include <iostream>
# include <iomanip>

   template<size_t Nnodes, typename ElemT, typename Functor, typename... ElemTs>
   void for_each_stencil( const std::array<Offset<3>,Nnodes> stencil,
                          const Functor                      functor,
                                MDArray<ElemT, 3>&              array,
                                MDArray<ElemTs,3>&...           arrays );

   int main()
  {
//    Idx<1> idx1{2};
//    Idx<2> idx2{3,4};

//    std::cout << idx1[0] << std::endl;
//    std::cout << idx2[0] << " " << idx2[1] << std::endl;

//    Offset<1> off1{2};
//    Offset<2> off2{21,13};

//    std::cout << (idx1+off1)[0] << std::endl;

//    Idx idx2p(idx2+off2);
//    std::cout << idx2p[0] << " " << idx2p[1] << std::endl;

//    std::cout << std::endl; // --------------------------------

      Dims<3> dims3{4,3,4};

//    Stride<3> stride3(dims3);

//    std::cout << stride3[0] << " "
//              << stride3[1] << " "
//              << stride3[2] << std::endl;

//    Idx<3> idx3{1,2,3};

//    std::cout << stride3*idx3 << std::endl;

//    std::cout << std::endl; // --------------------------------

      MDArray<double,3> array3(dims3);

      for( size_t i=0; i<dims3[0]; i++ )
     {
         for( size_t j=0; j<dims3[1]; j++ )
        {
            for( size_t k=0; k<dims3[2]; k++ )
           {
               array3[{i,j,k}] = 100*i + 10*j + k;
           }
        }
     }

      for( size_t i=0; i<dims3[0]; i++ )
     {
         for( size_t j=0; j<dims3[1]; j++ )
        {
            for( size_t k=0; k<dims3[2]; k++ )
           {
               std::cout << std::setw(4);
               std::cout << array3[{i,j,k}] << " ";
           }
            std::cout << std::endl;
        }
         std::cout << std::endl;
         std::cout << std::endl;
     }

      std::cout << std::endl; // --------------------------------

      Offset<3> here{  0, 0, 0};
      Offset<3> right{ 0, 0, 1};
      Offset<3> left{  0, 0,-1};
      Offset<3> front{ 0, 1, 0};
      Offset<3> back{  0,-1, 0};
      Offset<3> above{ 1, 0, 0};
      Offset<3> below{-1, 0, 0};

//    for( size_t i=0; i<dims3[0]; i++ )
//   {
//       for( size_t k=1; k<dims3[2]-1; k++ )
//      {
//          Idx<3> idx{i,1,k};

//          std::cout << "     " << std::setw(4) << array3[idx+above] << std::endl;
//          std::cout            << std::setw(4) << array3[idx+ left] << " "
//                               << std::setw(4) << array3[idx      ] << " "
//                               << std::setw(4) << array3[idx+right] << std::endl;
//          std::cout << "     " << std::setw(4) << array3[idx+below] << std::endl;
//          std::cout << std::endl;
//      }
//       std::cout << std::endl;
//   }

      const auto print_stencil = []( const std::array<double,7>& elems )->void
     {
         std::cout << std::setw(4) << elems[0] << std::endl;
         std::cout << std::setw(4) << elems[1] << " " << std::setw(4) << elems[2] << std::endl;
         std::cout << std::setw(4) << elems[3] << " " << std::setw(4) << elems[4] << std::endl;
         std::cout << std::setw(4) << elems[5] << " " << std::setw(4) << elems[6] << std::endl;
         std::cout << std::endl;
         return;
     };

      const auto print_stencil2 = [print_stencil]( const std::array<double,7> elems0, const std::array<double,7> elems1 )->void
     {
         print_stencil( elems0 );
         std::cout << std::endl;
         print_stencil( elems1 );
         return;
     };

      for_each_stencil<7>( {here,left,right,
                                 back,front,
                                 below,above},
                           print_stencil2,
                           array3,
                           array3 );

      return 0;
  }

   template<size_t Nnodes, int nDim, typename ElemT,
            typename Stencil_Indices=std::make_index_sequence<Nnodes>>
   std::array<ElemT,Nnodes> get_stencil_elems(        MDArray<ElemT,nDim>&               array,
                                               const Idx<nDim>                         centre,
                                               const std::array<Offset<nDim>,Nnodes>& stencil )
  {
      return get_stencil_elems( array, centre, stencil, Stencil_Indices{} );
  }

   template<size_t Nnodes, int nDim, typename ElemT,
            size_t... Stencil_Is>
   std::array<ElemT,Nnodes> get_stencil_elems(       MDArray<ElemT,nDim>&               array,
                                               const Idx<nDim>                         centre,
                                               const std::array<Offset<nDim>,Nnodes>& stencil,
                                               const std::index_sequence<Stencil_Is...> )
  {
      return { (array[ centre + stencil[Stencil_Is]])... };
  }

   template<size_t Nnodes, typename ElemT, typename Functor, typename... ElemTs>
   void for_each_stencil( const std::array<Offset<3>,Nnodes> stencil,
                          const Functor                      functor,
                                MDArray<ElemT, 3>&             array,
                                MDArray<ElemTs,3>&...          arrays )
  {
      for( size_t i=1; i<array.dims[0]-1; i++ )
     {
         for( size_t j=1; j<array.dims[1]-1; j++ )
        {
            for( size_t k=1; k<array.dims[2]-1; k++ )
           {
               Idx<3> idx{i,j,k};

               functor( get_stencil_elems( array,  idx, stencil ),
                        get_stencil_elems( arrays, idx, stencil )... );
           }
        }
     }
      return;
  }




