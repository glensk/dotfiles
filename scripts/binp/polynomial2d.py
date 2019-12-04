#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

import warnings
import numpy as np
import sys
import myplotutils as myp
import myutils as my
np.set_printoptions(suppress=True)   # display arrays withou 000000
np.set_printoptions(precision=6)    # print only 6 digist after .

#__all__ = ['polyval2d', 'polygrid2d', 'polyfit2d']
#
#
#class RankWarning(UserWarning):
#    pass


def _polyvander2d(x, y, deg):
    """Return the pseudo-Vandermonde matrix for a given degree.
    [[1, x[0], y[0], ... , x[0]*y[0]^(n-1), y[0]^n]
     [1, x[1], y[1], ... , x[1]*y[1]^(n-1), y[1]^n]
     ...                                       ...
     [1, x[M], y[M], ... , x[M]*y[M]^(n-1), y[M]^n]]
    where `n` is `deg`.
    Parameters
    ----------
    x, y : array_like, shape (M,)
        Array of points. `x` and `y` must have the same shape.
    deg : int
        Degree of the resulting matrix.
    Returns
    -------
    v : ndarray
        The Vandermonde matrix. The shape of the matrix is ``x.shape +
        ((deg+1)*(deg+2) // 2, ))``.
    """
    x = np.array(x, copy=False, ndmin=1) + 0.0
    y = np.array(y, copy=False, ndmin=1) + 0.0
    if x.ndim != 1:
        raise ValueError("x must be 1-dim.")
    if y.ndim != 1:
        raise ValueError("y must be 1-dim.")

    dims = ((deg+1)*(deg+2) // 2, ) + x.shape
    v = np.empty(dims, dtype=x.dtype)
    v[0] = x * 0 + 1.0
    i = 1
    for j in range(1, deg+1):
        v[i:i+j] = x * v[i-j:i]
        v[i+j] = y * v[i-1]
        i += j + 1
    return np.rollaxis(v, 0, v.ndim)


def polyfit2d(x, y, z, deg=1, rcond=None, full_output=False):
    """Return the coefficients of a polynomial of degree `deg`.
    The coefficients are determined by the least square fit to given data values
    `z` at given points ``(x, y)``.  The fitting assumes the polynomial in a
    form::
    .. math:: p(x,y) = \\sum_{i,j} c_{i,j} * x^i * y^j
    with a constraint of ``i + j <= n`` where `n` is `deg`.
    Parameters
    ----------
    x, y : array_like, shape (M,)
        x- and y-oordinates of the M data points ``(x[i], y[i])``.
    z : array_like, shape (M,)
        z-coordinates of the M data points.
    deg : int, optional
        Degree of the polynomial to be fit.
    rcond : float, optional
        Relative condition of the fit.  Singular values smaller than
        `rcond`, relative to the largest singular value, will be
        ignored.  The default value is ``len(x)*eps``, where `eps` is
        the relative precision of the platform's float type, about 2e-16
        in most cases.
    full_output : {True, False}, optional
        Just the coefficients are returned if False, and diagnostic
        information from the SVD is also returned if True.
    Returns
    -------
    coef : ndarray
        Array of coefficients.
    [residuals, rank, singular_values, rcond] : if `full_output` = True
        Sum of the squared residuals of the least-squares fit; the
        effective rank of the scaled pseudo-Vandermonde matrix; its
        singular values, and the specified value of `rcond`. For more
        details, see `numpy.linalg.lstsq`.
    Warns
    -----
    RankWarning
        The rank of the coefficient matrix in the least-squares fit is
        deficient.  The warning is only raised if `full_output` = False.
    See Also
    --------
    polyval2d, polygrid2d
    """
    x = np.asarray(x) + 0.0
    y = np.asarray(y) + 0.0
    z = np.asarray(z) + 0.0

    deg = int(deg)
    if deg < 1:
        raise ValueError("deg must be larger than 1.")

    # Check inputs.
    if x.ndim != 1:
        raise ValueError("x must be 1-dim.")
    if y.ndim != 1:
        raise ValueError("y must be 1-dim.")
    if z.ndim != 1:
        raise ValueError("z must be 1-dim.")
    if x.size != y.size or x.size != z.size:
        raise ValueError("x, y, and z must have the same size.")

    # Set up the matrices for the problem in transposed form.
    lhs = _polyvander2d(x, y, deg).T
    rhs = z.T

    # Set rcond.
    if rcond is None:
        rcond = x.size * np.finfo(x.dtype).eps

    # Determine the norms of the design maxtirx columns.
    if issubclass(lhs.dtype.type, np.complexfloating):
        scl = np.sqrt((np.square(lhs.real) + np.square(lhs.imag)).sum(1))
    else:
        scl = np.sqrt(np.square(lhs).sum(1))
    scl[scl == 0] = 1

    # Solve the least squares problem.
    c1, resids, rank, s = np.linalg.lstsq(lhs.T / scl, rhs.T, rcond)
    c1 = (c1.T / scl).T

    # Warn on rank reduction.
    if rank != lhs.shape[0] and not full_output:
        msg = "The fit may be poorly conditioned."
        warnings.warn(msg, RankWarning)

    # Allocate the coefficients in a 2-dim array.
    inds = []
    for m in range(deg + 1):
        for j in range(m + 1):
            for i in range(m + 1):
                if i + j != m:
                    continue
                inds.append((i, j))
    c2 = np.zeros((deg + 1, deg + 1))
    c2[zip(*inds)] = c1

    if full_output:
        return c2, [resids, rank, s, rcond]
    else:
        return c2


def polygrid2d(x, y, c):
    """
    Evaluate a 2-D polynomial on the Cartesian product of x and y.
    This function returns the values:
    .. math:: p(a,b) = \sum_{i,j} c_{i,j} * a^i * b^j
    where the points `(a, b)` consist of all pairs formed by taking
    `a` from `x` and `b` from `y`. The resulting points form a grid with
    `x` in the first dimension and `y` in the second.
    The parameters `x` and `y` are converted to arrays only if they are
    tuples or a lists, otherwise they are treated as a scalars. In either
    case, either `x` and `y` or their elements must support multiplication
    and addition both with themselves and with the elements of `c`.
    If `c` has fewer than two dimensions, ones are implicitly appended to
    its shape to make it 2-D. The shape of the result will be c.shape[2:] +
    x.shape + y.shape.
    Parameters
    ----------
    x, y : array_like, compatible objects
        The two dimensional series is evaluated at the points in the
        Cartesian product of `x` and `y`.  If `x` or `y` is a list or
        tuple, it is first converted to an ndarray, otherwise it is left
        unchanged and, if it isn't an ndarray, it is treated as a scalar.
    c : array_like
        Array of coefficients ordered so that the coefficients for terms of
        degree i,j are contained in ``c[i,j]``. If `c` has dimension
        greater than two the remaining indices enumerate multiple sets of
        coefficients.
    Returns
    -------
    values : ndarray, compatible object
        The values of the two dimensional polynomial at points in the Cartesian
        product of `x` and `y`.
    See Also
    --------
    polyval2d, polyfit2d
    numpy.polynomial.polynomial.polygrid2d
    """
    from numpy.polynomial.polynomial import polygrid2d

    c = np.asarray(c)
    if c.ndim != 2 or c.shape[0] != c.shape[1]:
        raise ValueError("c must be a squard 2-dim array.")
    return polygrid2d(x, y, c)


def polyval2d(x, y, c):
    """Evaluate a 2-D polynomial at points (x, y).
    This function returns the value
    .. math:: p(x,y) = \\sum_{i,j} c_{i,j} * x^i * y^j
    Parameters
    ----------
    x, y : array_like, compatible objects
        The two dimensional series is evaluated at the points (x, y), where x
        and y must have the same shape. If x or y is a list or tuple, it is
        first converted to an ndarray, otherwise it is left unchanged and, if it
        is not an ndarray, it is treated as a scalar.
    c : array_like
        Array of coefficients ordered so that the coefficient of the term of
        multi-degree i,j is contained in `c[i,j]`.
    Returnes
    --------
    values : ndarray, compatible object
        The values of the two dimensional polynomial at points formed with pairs
        of corresponding values from x and y.
    See Also
    --------
    polygrid2d, polyfit2d
    numpy.polynomial.polynomial.polyval2d
    """
    from numpy.polynomial.polynomial import polyval2d

    c = np.asarray(c)
    if c.ndim != 2 or c.shape[0] != c.shape[1]:
        raise ValueError("c must be a squared 2-dim array.")
    return polyval2d(x, y, c)


def xyz_line_to_XYZ(x,y,z):
    #np.set_printoptions(precision=4)
    #np.set_printoptions(suppress=False)
    my.fakefun()
    sys.exit()
    print('x',x)
    print('y',y)
    print('z',z)
    if x.ndim != 1:
        raise ValueError("x must be 1-dim.")
    if y.ndim != 1:
        raise ValueError("y must be 1-dim.")
    if z.ndim != 1:
        raise ValueError("z must be 1-dim.")
    if x.size != y.size or x.size != z.size:
        print('x.size',x.size)
        print('y.size',y.size)
        print('z.size',z.size)
        raise ValueError("x, y, and z must have the same size.")
    xx,yy,zz = myp.xyz_line_to_xyz_unique(x=x,y=y,z=z,data=None)
    sys.exit('k777')
    print('xx')
    print(xx)
    print('yy')
    print(yy)
    print('zz')
    print(zz)
    # diese x, y, und z sind wowas wie die unique elements!
    if False:
        x = np.array([100,100,200])
        y = np.array([11.1,11.2,11.3])
        z = np.array([0.1,0.2,0.3])

    data = np.zeros((len(x),3))
    data[:,0] = x
    data[:,1] = y
    data[:,2] = z
    print('data')
    print(data)
    X, Y = np.meshgrid(xx, yy, copy=False)
    Z = X**2 + Y**2
    #Z = np.zero(np.shape(X))
    print('X')
    print(X)
    print('Y')
    print(Y)
    print('Z')
    print(Z)

    X = X.flatten()
    Y = Y.flatten()
    B = Z.flatten()
    print('X',len(X))
    print('x',len(x))
    print(X)
    print(x)
    print('Y',len(Y))
    print(Y)
    print('B',len(Z))
    print(B)
    sys.exit('k888')
    A = np.array([X*0+1, X, Y, X**2, X**2*Y, X**2*Y**2, Y**2, X*Y**2, X*Y]).T
    coeff, r, rank, s = np.linalg.lstsq(A, B,rcond=None)

    def poly2Dreco(X, Y, c):
        return (c[0] + X*c[1] + Y*c[2] + X**2*c[3] + X**2*Y*c[4] + X**2*Y**2*c[5] +
           Y**2*c[6] + X*Y**2*c[7] + X*Y*c[8])

    print('done')
    print('coeff')
    print(coeff)
    out = poly2Dreco(X, Y, coeff)
    print('out')
    print(out)

    # make 2d plots
    if False:
        for idv,v in enumerate(np.unique(y)):
            allt = np.arange(0,x.max()+1,1)
            allz = polynomial2d.polyval2d(allt,v,c)
            #for idt,t in enumerate(np.arange(0,xgrid.max()+1,1)):
            #    z = polynomial2d.polyval2d(t,v,c)
            #    print(v,t,z)
            np.savetxt("out_new_"+str(v)+'fit.dat',np.array([allt,allz]).T)
            ka = data[np.where(data[:,1]==v)[0]][:,[0,2]]
            arr = ka[ka[:,0].argsort()]
            np.savetxt("out_new_"+str(v)+'dat.dat',arr)

    return

if __name__ == "__main__":
    my.fakefun()
    sys.exit()
    np.set_printoptions(precision=4)
    data = G = np.loadtxt("/Users/glensk/Downloads/cu/fah/thermo/Fah_surface")[:,[0,1,2]]
    x = xline = t = data[:,0];
    y = yline = v = data[:,1];
    z = zline = e = data[:,2];
    xyz_line_to_XYZ(x=x,y=y,z=z)
    df1 = myp.xyz_line_to_dataframe(x=x,y=y,z=z+0,xlabel="temp",ylabel='vol',zlabel='fah')

    # Build figure
    import plotly
    import plotly.graph_objs as go
    fig = go.Figure()
    fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))

    # Add first scatter trace with medium sized markers
    fig.add_trace(go.Scatter3d(mode='markers',x=x,y=y,z=z,
            opacity=0.9,marker=dict(color='red',size=5,line=dict(color='MediumPurple',width=2)),name='data'))
    fig.show()

    # try normal plot
    if False: # works
        print('nn')
        import matplotlib.pyplot as plt
        plt.plot([1, 2, 3, 4])
        plt.ylabel('some numbers')
        plt.show()

    if False: # does not work
        import plotly.express as px

        iris = px.data.iris()
        fig = px.scatter(iris, x="sepal_width", y="sepal_length", color="species")

        fig.update_traces(marker=dict(size=12,
                                      line=dict(width=2,
                                                color='DarkSlateGrey')),
                          selector=dict(mode='markers'))
        fig.show()

    if False:
        import plotly.graph_objects as go
        fig = go.Figure(
            data=[go.Bar(y=[2, 1, 3])],
            layout_title_text="A Figure Displayed with fig.show()"
        )
        fig.show()

