import numpy as np

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
        r1[ 0]+= gamma*( r0[-1] - r0[ 0] )
        r1[-1]+= gamma*( r0[-2] - r0[-1] )

    return r1


