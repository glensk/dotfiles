#!/usr/bin/env python
 # -*- coding: utf-8 -*-
from __future__ import print_function
import numpy as np
import pandas as pd
import os,sys,argparse

#def help(p = None):
#    string = ''' helptext '''
#    p = argparse.ArgumentParser(description=string,
#            formatter_class=argparse.RawTextHelpFormatter)
#    p.add_argument('-i','--inputfile', required=True, type=str,default=False, help="name of the inputfile inputfile")
#    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)
#    return p
#
#def todo(args):
#    if not os.path.isfile(args.inputfile):
#        sys.exit("inputfile "+args.inputfile+" does not exist")
#    return

##########################################################
# prapare data for plots
##########################################################
def xyz_line_to_dataframe(x=False,y=None,z=None,data=None,xlabel='x',ylabel='y',zlabel='z'):
    if data is not None:
        x = data[:,0]
        y = data[:,1]
        z = data[:,2]
    df = np.array([x,y,z]).T
    df = pd.DataFrame(df)
    df.columns = [xlabel,ylabel,zlabel]
    return df

def xyz_line_to_xyz_unique(x=False,y=None,z=None,data=None):
    '''
    input:
    -------
    x, y,z : array-like, 1d

    output:
    -------
    x, y: array-like, 1d
    z   : array-like, 2d
    '''
    if data is not None:
        x = data[:,0]
        y = data[:,1]
        z = data[:,2]
    else:
        data = np.zeros((len(x),3))
        data[:,0] = x
        data[:,1] = y
        data[:,2] = z

    xunique = np.unique(x)
    yunique = np.unique(y)
    X, Y = np.meshgrid(xunique, yunique, copy=False)
    zz = np.zeros((len(xunique),len(yunique)))
    #xx = np.zeros((len(xunique),len(yunique)))
    for idx,i in enumerate(data):
        xidx=np.where(data[idx,0]==xunique)[0][0]
        yidx=np.where(data[idx,1]==yunique)[0][0]
        zz[xidx,yidx]=data[idx,2]
        #xx[xidx,yidx]=data[idx,2]
        print('-->',data[idx],xidx,yidx)
    print('X')
    print(X)
    print('Y')
    print(Y)
    print('zz')
    print(zz)
    sys.exit('tt')
    return xunique,yunique,zz

##########################################################
# make plots
##########################################################
def polyfit2d(x, y, z, kx=3, ky=3, order=None):
    '''
    Two dimensional polynomial fitting by least squares.
    Fits the functional form f(x,y) = z.

    Notes
    -----
    Resultant fit can be plotted with:
    np.polynomial.polynomial.polygrid2d(x, y, soln.reshape((kx+1, ky+1)))

    Parameters
    ----------
    x, y: array-like, 1d
        x and y coordinates.
    z: np.ndarray, 2d
        Surface to fit.
    kx, ky: int, default is 3
        Polynomial order in x and y, respectively.
    order: int or None, default is None
        If None, all coefficients up to maxiumum kx, ky, ie. up to and including x^kx*y^ky, are considered.
        If int, coefficients up to a maximum of kx+ky <= order are considered.

    Returns
    -------
    Return paramters from np.linalg.lstsq.

    soln: np.ndarray
        Array of polynomial coefficients.
    residuals: np.ndarray
    rank: int
    s: np.ndarray

    '''

    # grid coords
    x, y = np.meshgrid(x, y)
    # coefficient array, up to x^kx, y^ky
    coeffs = np.ones((kx+1, ky+1))

    # solve array
    a = np.zeros((coeffs.size, x.size))

    # for each coefficient produce array x^i, y^j
    for index, (j, i) in enumerate(np.ndindex(coeffs.shape)):
        # do not include powers greater than order
        if order is not None and i + j > order:
            arr = np.zeros_like(x)
        else:
            arr = coeffs[i, j] * x**i * y**j
        a[index] = arr.flatten()

    # do leastsq fitting and return leastsq result
    return np.linalg.lstsq(a.T, np.ravel(z), rcond=None)

def make_nice_2D_scatterplot(df,tags=None,x="a",y="b",color=None,symbols=False):
    ''' this is a 2D scatter plot

        make_nice_scatterplot(dataframe,tags=None,x="a",y="b",color=range(ntot))
        df is a pandas dataframe e.g.
        df = pd.DataFrame(projs) # where projs is a numpy array
        df.columns = ["a","b"]
        make_nice_scatterplot(df,tags=tags,x="a",y="b",color=range(ntot))

        or:
        x = np.arange(len(rd))
        y = rd  # error in KMC structures
        color = np.log(fd)
        from myutils import make_nice_scatterplot as sp
        dfa = np.array([x,y]).T
        import pandas as pd
        df = pd.DataFrame(dfa)
        df.columns = ["struct","error"]
        sp(df,x="struct",y="error",color=color)

        example in: fps_correlation_scatterplot.ipynb


    '''
    import plotly.graph_objects as go
    fig = go.Figure(data=go.Scatter(x=df[x],
                                    y=df[y],
                                    mode='markers',
                                    text=tags,
                                    marker=dict(
                                        size=5,
                                        color=color, #set color equal to a variable
                                        colorscale='Viridis', # one of plotly colorscales
                                        showscale=True
                                    )
                                    )) # hover text goes here
    if symbols != False:
    	fig['data'][0]['marker']['symbol'] = 'triangle-left'
    size = 400
    fig.update_layout(
        autosize=False,
        width=size,
        height=size*0.5,
        margin=go.layout.Margin(
            l=0,
            r=0,
            b=0,
            t=0,
            pad=0
       ),
        #paper_bgcolor="LightSteelBlue",
    )
    fig.show()
    return

if __name__ == '__main__':
    p = help()
    args = p.parse_args()
    todo(args)


