
# pragma once

# include <utils/concepts.h>

# include <array>

namespace ODE
{
   namespace Implicit
  {
      template<floating_point Real>
      struct MultiStep
     {
   // du/dt = R(u)

   // ( b_0*u_{n+1} + b_1*u_{n} + b_2*u_{n-1} )/dt
   //       = g_0*R( u_{n+1} ) + g_1*R( u_{n} )

      // order of accuracy
         unsigned int order;

      // number of nsteps
         unsigned int nsteps;

      // number of residual evaluations
         unsigned int nresid;

      // coefficients
         std::array<Real,3>  beta;
         std::array<Real,3> gamma;
     };

      template<floating_point Real>
      MultiStep<Real> backwardDifference1()
     {
         MultiStep<Real> ms{};
         ms.order  = 1;
         ms.nsteps = 2;
         ms.nresid = 1;

         ms.beta[0] = 1.0;
         ms.beta[1] =-1.0;

         ms.gamma[0]=1.0;

         return ms;
     }

      template<floating_point Real>
      MultiStep<Real> backwardDifference2()
     {
         MultiStep<Real> ms{};
         ms.order  = 2;
         ms.nsteps = 3;
         ms.nresid = 1;

         ms.beta[0] = 1.5;
         ms.beta[1] =-2.0;
         ms.beta[2] = 0.5;

         ms.gamma[0]=1.0;

         return ms;
     }

      template<floating_point Real>
      MultiStep<Real> trapeziumRule2()
     {
         MultiStep<Real> ms{};
         ms.order  = 2;
         ms.nsteps = 2;
         ms.nresid = 2;

         ms.beta[0] = 1.0;
         ms.beta[1] =-1.0;

         ms.gamma[0]=0.5;
         ms.gamma[1]=0.5;

         return ms;
     }
  }

   namespace Explicit
  {
      template<floating_point Real>
      struct RungeKutta
     {
      // order of accuracy
         unsigned int order;

      // number of stages
         unsigned int nstages;

      // maximum allowable cfl
         Real maxCFL;

      // coefficients
      // alpha[s][t] is coefficient of residual t at stage s
         std::array<std::array<Real,6>,6> alpha;

      // beta[s] is the size of the timestep of stage s
         std::array<Real,6> beta;
     };

      template<floating_point Real>
      RungeKutta<Real> ssp11()
     {
         RungeKutta<Real> rk{};

         rk.order  =1;
         rk.nstages=1;
         rk.maxCFL =1.;

         rk.beta[0]=1.0;

         rk.alpha[0][0]=1.0;

         return rk;
     }

      template<floating_point Real>
      RungeKutta<Real> ssp22()
     {
         RungeKutta<Real> rk{};

         rk.order=2;
         rk.nstages=2;
         rk.maxCFL=1.;

         rk.beta[0]=1.0;
         rk.beta[1]=1.0;

         rk.alpha[0][0]=1.0;

         rk.alpha[1][0]=0.5;
         rk.alpha[1][1]=0.5;

         return rk;
     }

      template<floating_point Real>
      RungeKutta<Real> ssp33()
     {
         RungeKutta<Real> rk{};

         rk.order=3;
         rk.nstages=3;
         rk.maxCFL=1.;

         rk.beta[0]=1.0;
         rk.beta[1]=0.5;
         rk.beta[2]=1.0;

         rk.alpha[0][0]=1.;

         rk.alpha[1][0]=0.5;
         rk.alpha[1][1]=0.5;

         rk.alpha[2][0]=1./6.;
         rk.alpha[2][1]=1./6.;
         rk.alpha[2][2]=2./3.;

         return rk;
     }

      template<floating_point Real>
      RungeKutta<Real> ssp34()
     {
         RungeKutta<Real> rk{};

         rk.order=3;
         rk.nstages=4;
         rk.maxCFL=2.;

         rk.beta[0]=0.5;
         rk.beta[1]=1.0;
         rk.beta[2]=0.5;
         rk.beta[3]=1.0;

         rk.alpha[0][0]=1.0;

         rk.alpha[1][0]=0.5;
         rk.alpha[1][1]=0.5;

         rk.alpha[2][0]=1./3.;
         rk.alpha[2][1]=1./3.;
         rk.alpha[2][2]=1./3.;

         rk.alpha[3][0]=1./6.;
         rk.alpha[3][1]=1./6.;
         rk.alpha[3][2]=1./6.;
         rk.alpha[3][3]=1./2.;

         return rk;
     }
  }
}

