
# include <conservationLaws/scalarAdvection/scalarAdvection.h>
# include <conservationLaws/euler/euler.h>

# include <geometry/geometry.h>

# include <vector>
# include <iostream>
# include <fstream>

constexpr LawType Law = LawType::ScalarAdvection;
constexpr int nDim = 1;

using BasisT = BasisType<Law>;

constexpr BasisT SolutionBasis = BasisT::Conserved;

using SolutionVarSet   = VariableSet<  Law,nDim,SolutionBasis>;
using SolutionVarDelta = VariableDelta<Law,nDim,SolutionBasis>;

using FluxRes = FluxResult<Law,nDim>;

using Flux = RusanovFlux<Law>;

using Point   = Geometry::Point<nDim>;
using Surface = Geometry::Surface<nDim>;
using Volume  = Geometry::Volume<nDim>;

   auto get_species( ScalarAdvectionBases Basis )
  {
      return Species<LawType::ScalarAdvection>{};
  }

   auto get_species( EulerBases Basis )
  {
      return get_air_species();
  }

   auto initialLeft( ScalarAdvectionBases Basis )
  {
      using VarSet = VariableSet<LawType::ScalarAdvection,
                                 1,
                                 ScalarAdvectionBases::Conserved>;
      constexpr Types::Real u=1.;
      constexpr Types::Real q=2.;
      return VarSet{u,q};
  }

   auto initialRight( ScalarAdvectionBases Basis )
  {
      using VarSet = VariableSet<LawType::ScalarAdvection,
                                 1,
                                 ScalarAdvectionBases::Conserved>;
      constexpr Types::Real u=1.;
      constexpr Types::Real q=1.;
      return VarSet{u,q};
  }

   auto initialLeft( EulerBases Basis )
  {
      using VarSet = VariableSet<LawType::Euler,
                                 1,
                                 EulerBases::Conserved>;
      constexpr Types::Real r =1.0;
      constexpr Types::Real ru=0.0;
      constexpr Types::Real re=2.5;
      return VarSet{ru,r,re};
  }
   
   auto initialRight( EulerBases Basis )
  {
      using VarSet = VariableSet<LawType::Euler,
                                 1,
                                 EulerBases::Conserved>;
      constexpr Types::Real r =0.125;
      constexpr Types::Real ru=0.0;
      constexpr Types::Real re=0.25;
      return VarSet{ru,r,re};
  }

   void writeState( std::ofstream& os, const State<LawType::ScalarAdvection,1>& s )
  {
      os << s.velocity(0) << " "
         << s.scalar()    << std::endl;
      return;
  }
   
   void writeState( std::ofstream& os, const State<LawType::Euler,1>& s )
  {
      os << s.density()   << " "
         << s.velocity(0) << " "
         << s.pressure()  << std::endl;
      return;
  }
   

   int main()
  {
   // parameters

      // number of cells and timesteps
      const int nx=64;
      const int nt=64;

      // courant-friedrichs-lewy number
      const Types::Real cfl=0.5;

   // setup
      const Species<Law> species = get_species( SolutionBasis );
      const Flux         flux{};

   // solutions arrays
      std::vector<SolutionVarSet>  q(nx);
      std::vector<SolutionVarSet> q1(nx);
      std::vector<FluxRes>       res(nx);

   // initialise mesh
      const std::vector<Point> nodes = [&]()
     {
         std::vector<Point> nds(nx+1);
         for( int i=0; i<nx+1; i++ ){ nds[i][0]=i; }
         return nds;
     }();
      const std::vector<Volume> cells = Geometry::dual( nodes );

   // initialise solution to left/right states
      const SolutionVarSet ql0 = initialLeft(  SolutionBasis );
      const SolutionVarSet qr0 = initialRight( SolutionBasis );
      std::cout << "left  state: " << ql0 << std::endl;
      std::cout << "right state: " << qr0 << std::endl;
      for( int i=0;    i<nx/2; i++ ){ q[i]=ql0; }
      for( int i=nx/2; i<nx;   i++ ){ q[i]=qr0; }

   // integrate forward in time
      for( int t=0; t<nt; t++ )
     {
      // reset residual
         for( FluxRes& q0 : res ){ q0 = FluxRes{}; }

      // accumulate cell residual contributions from the flux across each face
         for( int i=0; i<nx-1; i++ )
        {
            const SolutionVarSet ql=q[i];
            const SolutionVarSet qr=q[i+1];

            const Surface face = surface( nodes[i] );

            const FluxRes fr = flux( species, face,  ql, qr );

            res[i]  -=fr;
            res[i+1]+=fr;
        }

      // left boundary : dirichlet condition at initial state
        {
            const int i=0;
            const SolutionVarSet qr=q[i];
            const Surface face = surface( nodes[i] );
            const FluxRes fr = flux( species, face, ql0, qr );
            res[i]+=fr;
        }

      // right boundary : outflow, upwind flux
        {
            const int i=nx-1;
            const SolutionVarSet ql=q[i];
            const Surface face = surface( nodes[i] );
            const FluxRes fr = flux( species, face, ql, ql );
            res[i]-=fr;
        }

      // calculate maximum stable timestep
         const Types::Real dt = [&]()
        {
            Types::Real dmin = std::numeric_limits<Types::Real>::max();
            for( int i=0; i<nx; i++ )
           {
               const Types::Real d = cells[i].volume
                                    /res[i].lambda;
               dmin = std::min( dmin, d );
           }
            return cfl*dmin;
        }();

      // integrate cell residuals forward by dt and average over cell volume
         for( int i=0; i<nx; i++ )
        {
            const Types::Real d = dt/cells[i].volume;

            const SolutionVarDelta dq = d*res[i].flux;

            q[i]+=dq;
        }
     }

      std::ofstream solutionFile( "data/testSA.dat" );
      if( solutionFile.is_open() )
     {
         for( const SolutionVarSet& q0 : q )
        {
            writeState( solutionFile, set2State( species, q0 ) );
        }
     }
      else
     {
         std::cout << "cannot open \"testSA.dat for writing\"\n" << std::endl;
         return 1;
     }
      solutionFile.close();

      return 0;
  }

/*
         grad = gradients( cells, q );

         res = fluxResidual(        nodes, flux,          species,       q );
         res = fluxResidual( cells, nodes, flux, limiter, species, grad, q );

         const Types::Real dt = timestep( cells, res, cfl );

         q1 = eulerForwardUpdate( cells, species, dt, q, res );
         std::swap( q,q1 );
*/


