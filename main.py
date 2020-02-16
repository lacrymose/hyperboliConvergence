
import numpy as np
import matplotlib.pyplot as plt
import timeseries_plotting as tsp

def difference2( x, u, isPeriodic ):
    """
    calculates second order difference operator of field u over periodic domain x
    x: point coordinate array
    u: solution value array
    """
    du=np.zeros_like( u )

    du[1:-1] = u[2:] - u[:-2]

    if isPeriodic:
        du[-1] = u[ 0] - u[-2]
        du[ 0] = u[ 1] - u[-1]
    else:
        du[-1] = u[-1] - u[-2]
        du[-1]*= 2.
        du[ 0] = 0.

    return 0.5*du

def difference4( x, u, isPeriodic ):
    """
    calculates second order fourth difference operator of field u over periodic domain x
    x: point coordinate array
    u: solution value array
    """
    du=np.zeros_like( u )

    du[2:-2] = u[:-4] - 4*u[1:-3] + 6*u[2:-2] - 4*u[3:-1] + u[4:]

    if isPeriodic:
        du[-2] = u[-4] - 4*u[-3] + 6*u[-2] - 4*u[-1] + u[0]
        du[-1] = u[-3] - 4*u[-2] + 6*u[-1] - 4*u[ 0] + u[1]
        du[ 0] = u[-2] - 4*u[-1] + 6*u[ 0] - 4*u[ 1] + u[2]
        du[ 1] = u[-1] - 4*u[ 0] + 6*u[ 1] - 4*u[ 2] + u[3]
    else:
        du[-2] = du[-3]
        du[-1] = du[-3]
        du[ 0] = 0.
        du[ 1] = du[ 2]

    return du

def cirs( mdt, r0, isPeriodic ):
    """
    implicit residual smoothing of r to increase timestep by factor of mdt.
    smoothing coefficient according to Enander1993
    timestep adjustment for fourier symbol overshoot of 1.2
    uses jacobi iterations to solve implicit system
    """
    beta = 0.25*( mdt*mdt - 1)

    r2=r0
    for i in range(100):
        r1 =r2
        r2 = mdt*r0/1.2

        r2[1:-1]+= beta*( r1[:-2] + r1[2:] )
        r2[1:-1]/= 1+2*beta

        if isPeriodic:
            r2[ 0]+= beta*( r1[ -1] + r1[1 ] )
            r2[-1]+= beta*( r1[ -2] + r1[0 ] )
            r2[ 0]/= 1+2*beta
            r2[-1]/= 1+2*beta
        else:
            r2[-1]+= beta*( r1[-2] )
            r2[-1]/= 1+beta

    return r2

def cers( mdt, r0, isPeriodic ):
    """
    implicit residual smoothing of r to increase timestep by factor of mdt.
    uses jacobi iterations to solve implicit system
    """
    beta  = 0.25*( mdt*mdt - 1 )
    gamma = 0.25*( ( ( 1+4*beta )/mdt ) - 1 )

#   r1=( 1+2*gamma )*r0
#   r1[1:-1]+= gamma*( -r0[:-2] - r0[2:] )

    r1=r0
    r1[1:-1]+= gamma*( r0[ :-2] - r0[1:-1] )
    r1[1:-1]+= gamma*( r0[2:  ] - r0[1:-1] )
    if isPeriodic:
        r1[ 0]+= gamma*( r0[-1] - r0[ 0] )
        r1[ 0]+= gamma*( r0[ 1] - r0[ 0] )

        r1[-1]+= gamma*( r0[-2] - r0[-1] )
        r1[-1]+= gamma*( r0[ 0] - r0[-1] )
    else:
        r1[-1]+= gamma*( r0[-2] - r0[-1] )

    return r1

# parameters
L   = 2*np.pi
n   = 128
h   = L/n
T   = 1.0*L
cfl = 1.5
omega=1

isPeriodic=True

# fourth order smoothing coefficient
l4  = 1./32

# residual smoothing timestep increase
mdt = 2

# runge-kutta coefficients
alpha=np.array( [0.6, 0.6, 1.0] )

# initialise
x = np.linspace( 0, L-h, n )
if isPeriodic:
    u0= np.cos( omega*x )
else:
    u0= 1. - x/(2.*np.pi)

dt = cfl*h
nt = int(T/dt)

r = np.zeros_like( u0 )
u = u0

residuals=np.zeros( nt )
history = np.zeros( ( nt+1, len(u0) ) )
history[0,:] = u0

# time stepping
for i in range(nt):
    u0 = u
    for k in range( len(alpha) ):
        r0 =    difference2( x, u, isPeriodic )
        r0+= l4*difference4( x, u, isPeriodic )
        r0/= h

        r1= cers( mdt, r0, isPeriodic )
        r = cirs( mdt, r1, isPeriodic )
#       r = cirs( mdt, r0, isPeriodic )

        u = u0 - alpha[k]*dt*r

    residuals[i] = np.sum( np.abs( u-u0 ) )/dt
    history[i+1,:]=u

tsp.view_timeseries1D( x, history )

#plt.semilogy( residuals )
#plt.show()








