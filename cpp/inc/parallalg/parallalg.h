
# pragma once

namespace par
{

/*
 * indicates type of grid that an array is defined over - primal or dual
 *    Primal is standard grid (should be default unless both primal and dual are needed)
 *    Dual grid used for cell based grids
 */
   enum struct GridType { Primal, Dual };
   inline constexpr GridType Primal = GridType::Primal;
   inline constexpr GridType   Dual = GridType::Dual;

/*
 * An NDIM-dimensional index.
 *    Types indexable with Idx are compatible with parallalg algorithms
 */
   template<int          NDIM,
            GridType GRIDTYPE= Primal>
   struct Idx
  {
      std::array<size_t,NDIM> idxs;
      const size_t& operator[]( const unsigned int i ) const { return idxs[i]; }
  };

// convenience typedefs

   // grid type
   template<int NDIM>
   using PrimalIdx = Idx<NDIM,Primal>;

   template<int NDIM>
   using DualIdx = Idx<NDIM,Dual>;

   // dimension
   template<GridType GRIDTYPE>
   using Idx1 = Idx<1,GRIDTYPE>;

   template<GridType GRIDTYPE>
   using Idx2 = Idx<2,GRIDTYPE>;

   template<GridType GRIDTYPE>
   using Idx3 = Idx<3,GRIDTYPE>;

   // grid type and dimension
   using PrimalIdx1 = Idx<1,Primal>;
   using PrimalIdx2 = Idx<2,Primal>;
   using PrimalIdx3 = Idx<3,Primal>;

   using DualIdx1 = Idx<1,Dual>;
   using DualIdx2 = Idx<2,Dual>;
   using DualIdx3 = Idx<3,Dual>;

/*
 * namespace holding classes and functions used to determine how algorithms are parallelised
 */
   namespace execution
  {
   /*
    * unique types indicating how algorithms can be parallelised
    */
      struct serial_policy {};
      inline constexpr serial_policy seq;
   
      struct openmp_policy {};
      inline constexpr openmp_policy omp;
  }

   template<typename T>
   struct is_execution_policy : std::false_type {};

   template<> struct is_execution_policy<execution::serial_policy> : std::true_type {};
   template<> struct is_execution_policy<execution::openmp_policy> : std::true_type {};

   template<typename T>
   inline constexpr bool is_execution_policy_v = is_execution_policy<T>::value;

   template<typename T>
   concept bool execution_policy = is_execution_policy_v<T>;
}
