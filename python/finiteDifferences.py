import numpy as np

def difference2( x, u, isPeriodic ):
    """
    calculates second order difference operator of field u over periodic domain x
    x: point coordinate array
    u: solution value array
    """
    du=np.zeros_like( u )

    du[1:-1] = 0.5*( u[2:] - u[:-2] )

    if isPeriodic:
        du[-1] = 0.5*( u[ 0] - u[-2] )
        du[ 0] = 0.5*( u[ 1] - u[-1] )
    else:
        du[-1] = u[-1] - u[-2]
        du[ 0] = u[ 1] - u[ 0]

    return du

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
        du[ 0] = du[ 2]
        du[ 1] = du[ 2]

    return du


