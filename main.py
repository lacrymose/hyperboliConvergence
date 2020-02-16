
import numpy as np
import matplotlib.pyplot as plt
import timeseries_plotting as tsp

import residualSmoothing as rs
import finiteDifferences as fd

# parameters
L   = 2*np.pi
n   = 64
h   = L/n
T   = 0.1*L
cfl = 1.0
omega=1

isPeriodic=False
smoothed=False

# fourth order smoothing coefficient
l4  = 1./32

# residual smoothing timestep increase
mdt = 10

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
    # runge kutta iterations
    for k in range( len(alpha) ):
        # flux residual calculation
        r0 =    fd.difference2( x, u, isPeriodic )
        r0+= l4*fd.difference4( x, u, isPeriodic )
        r0/= h

        if smoothed:
            r0[ 0]=0.
#           r0[-1]=0.
            r1= rs.cers( mdt, r0, isPeriodic )
            r1[ 0]=0.
#           r1[-1]=0.
            r = rs.cirs( mdt, r1, isPeriodic )

#           r = rs.cers( mdt, r0, isPeriodic )

#           r = rs.cirs( mdt, r0, isPeriodic )
        else:
            r = r0

        # evolve
        u = u0 - alpha[k]*dt*r
        # boundary conditions
        u[0]=1.
#       u[-1] = u[-2] + ( u[-2] - u[-3] )

    residuals[i] = np.sum( np.abs( u-u0 ) )/dt
    history[i+1,:]=u

tsp.view_timeseries1D( x, history )

#plt.semilogy( residuals )
#plt.show()








