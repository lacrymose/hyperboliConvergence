# ifndef ODE_H
# define ODE_H

# include <types.h>

namespace ODE
{
   namespace Implicit
  {
      struct MultiStep
     {
   // du/dt = R(u)

   // ( b_0*u_{n+1} + b_1*u_{n} + b_2*u_{n-1} )/dt
   //       = g_0*R( u_{n+1} ) + g_1*R( u_{n} )

      // order of accuracy
         int order;

      // number of nsteps
         int nsteps;

      // number of residual evaluations
         int nresid;

      // coefficients
         Types::Real  beta[3]={0,0,0};
         Types::Real gamma[2]={0,0};

         inline void eulerBackward1()
        {
            order=1;
            nsteps=2;
            nresid=1;

            beta[0]= 1.;
            beta[1]=-1.;

            gamma[0]=1.;
        }

         inline void backwardDifference2()
        {
            order=2;
            nsteps=3;
            nresid=1;

            beta[0]= 1.5;
            beta[1]=-2.0;
            beta[2]= 0.5;

            gamma[0]=1.;
        }

         inline void trapeziumRule2()
        {
            order=2;
            nsteps=2;
            nresid=2;

            beta[0]= 1.;
            beta[1]=-1.;

            gamma[0]=0.5;
            gamma[1]=0.5;
        }
     };
  }

   namespace Explicit
  {
      struct RungeKutta
     {
      // order of accuracy
         int order;

      // number of stages
         int nstages;

      // maximum allowable cfl
         Types::Real maxCFL;

      // coefficients
      // alpha[s][t] is coefficient of residual t at stage s
         Types::Real alpha[6][6]= {{0,0,0,0,0,0},
                             {0,0,0,0,0,0},
                             {0,0,0,0,0,0},
                             {0,0,0,0,0,0},
                             {0,0,0,0,0,0},
                             {0,0,0,0,0,0}};

      // beta[s] is the size of the timestep of stage s
         Types::Real beta[6]={0,0,0,0,0,0};

         inline void ssp11()
        {
            order=1;
            nstages=1;
            maxCFL=1.;

            beta[0]=1.0;

            alpha[0][0]=1.0;
        }

         inline void ssp22()
        {
            order=2;
            nstages=2;
            maxCFL=1.;

            beta[0]=1.0;
            beta[1]=1.0;

            alpha[0][0]=1.0;

            alpha[1][0]=0.5;
            alpha[1][1]=0.5;
        }

         inline void ssp33()
        {
            order=3;
            nstages=3;
            maxCFL=1.;

            beta[0]=1.0;
            beta[1]=0.5;
            beta[2]=1.0;

            alpha[0][0]=1.;

            alpha[1][0]=0.5;
            alpha[1][1]=0.5;

            alpha[2][0]=1./6.;
            alpha[2][1]=1./6.;
            alpha[2][2]=2./3.;
        }

         inline void ssp34()
        {
            order=3;
            nstages=4;
            maxCFL=2.;

            beta[0]=0.5;
            beta[1]=1.0;
            beta[2]=0.5;
            beta[3]=1.0;

            alpha[0][0]=1.0;

            alpha[1][0]=0.5;
            alpha[1][1]=0.5;

            alpha[2][0]=1./3.;
            alpha[2][1]=1./3.;
            alpha[2][2]=1./3.;

            alpha[3][0]=1./6.;
            alpha[3][1]=1./6.;
            alpha[3][2]=1./6.;
            alpha[3][3]=1./2.;
        }
     };
  }
}

# endif
