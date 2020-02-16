import numpy              as np
import matplotlib.pyplot  as plt

def previous_slice( ax ):
    ax.idx = ( ax.idx -1 ) % ax.field.shape[0]
    ax.lines[0].set_ydata( ax.field[ax.idx,:] )

def next_slice( ax ):
    ax.idx = ( ax.idx +1 ) % ax.field.shape[0]
    ax.lines[0].set_ydata( ax.field[ax.idx,:] )

def process_key( event ):
    fig = event.canvas.figure
    ax = fig.axes[0]
    if event.key == 'up':
        previous_slice( ax )
    elif event.key == 'down':
        next_slice( ax )
    fig.canvas.draw()

def process_scroll( event ):
    fig = event.canvas.figure
    ax = fig.axes[0]
    if event.button == 'up':
        previous_slice( ax )
    elif event.button == 'down':
        next_slice( ax )
    fig.canvas.draw()

def view_timeseries1D( domain, field ):
    fig, ax = plt.subplots()
    ax.idx = 0
    ax.field = field
    ax.domain = domain

    lo=np.min( field )
    hi=np.max( field )
    if( lo > 0 ):
        lo = 0.95*lo
    else:
        lo = 1.05*lo

    if( hi > 0 ):
        hi = 1.05*hi
    else:
        hi = 0.95*hi

    ax.set_ylim( [ lo, hi ] )
    ax.plot( domain, field[0,:] )
    fig.show()

    fig.canvas.mpl_connect(    'scroll_event', process_scroll )
    fig.canvas.mpl_connect( 'key_press_event', process_key    )

