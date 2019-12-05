#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

import warnings
import numpy as np
import sys,os
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

def polyfit2dbetter_DOES_NOT_WORK(x, y, z, kx=3, ky=3, order=None):
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
    #x, y = np.meshgrid(x, y)
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

def lmfit2d(x,y,z,verbose=False):
    if verbose:
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
    data = np.zeros((len(x),3))
    data[:,0] = x
    data[:,1] = y
    data[:,2] = z

    if True:
        all1v = np.unique(data[:,1])
        for v in all1v:
            data = np.append(data, [[0, v, 0]], axis=0)
    #print('data')
    #print(data)
    #sys.exit()

    def exp1(x,y):
        return x*0+1, x,       y,       x**2,       x**2*y,       x**2*y**2,       y**2,       x*y**2,       x*y

    def exp1sw(x,y):
        return x*0+1, y,       x**2,       y**2,       x**2*y,       x**2*y**2,       x,       x*y**2,       x*y

    def exp2(x,y):
        return x,  y,  x**2,       x**2*y,       x**2*y**2,       y**2,       x*y**2,       x*y

    def exp2ca(xy,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17):
        x = xy[:,0]
        y = xy[:,1]
        return c1*x**0*y**1 + \
               c2*x**0*y**2 + \
               c3*x**1*y**0 + \
               c4*x**1*y**1 + \
               c5*x**1*y**2 + \
               c6*x**2*y**0 + \
               c7*x**2*y**1 + \
               c8*x**2*y**2 + \
               c9 *x**3*y**0 + \
               c10*x**3*y**1 + \
               c11*x**3*y**2 + \
               c12*x**3*y**3 + \
               c13*x**4*y**0 + \
               c14*x**4*y**1 + \
               c15*x**4*y**2 + \
               c16*x**4*y**3 + \
               c17*x**4*y**4 + \
               0

    #def exp2c(xy,c):
    #    x = xy[:,0]
    #    y = xy[:,1]
    #    return x*c[0] + y*c[1] + c[2]*x**2 + c[3]*x**2*y+c[4]*x**2*y**2+c[5]*y**2+c[6]*x*y**2+c[7]*x*y

    #def exp3(x,y):
    #    return x,y,x**2,x**3,x**2*y,x**2*y**2,y**2, y**3,x*y**2,x*y

    if True:
        from lmfit import Model,minimize
        func = exp2ca
        model = Model(func)
        parameter_names = model.param_names
        independent_variable = model.independent_vars
        #print('parameter_names',parameter_names)
        #print('independent_variable',independent_variable)
        #print('parameter_names',type(parameter_names))
        for i in parameter_names: model.set_param_hint(i ,value=1)
        result = model.fit(data[:, 2], xy=data[:, 0:2])
        print(result.fit_report())
        #print('vgl (obtained by fit)',result.best_values)
        coef_lmfit = []
        for i in parameter_names: coef_lmfit.append(result.best_values.get(i))
        print('coef_lmfit',coef_lmfit)
        #print()
        #print(result.best_fit)
        diffmax = np.abs(result.best_fit-data[:,2]).max()
        print('diffmax',diffmax)
        print('888',data[0,0:2])
        denset = np.arange(0,x.max()+1,1)
        all1v = np.unique(data[:,1])
        for v in all1v:
            densev = np.repeat(v,len(denset))
            dd = np.array([denset,densev]).T
            out2dense = func(dd,*coef_lmfit)
            #print(out2dense)
            np.savetxt("out_OUTLMdense_"+str(v)+'dat.dat',out2dense)
        #print(exp2c(np.array([data[0,0:2]]),*coef_lmfit)[0])
    return

def oldfit_lstsq():
    expuse = exp2
    def expression(x,y,exprr,coefs=None):
        if len(x) != len(y):
            print('len(x)',len(x))
            print('len(y)',len(y))
            sys.exit('lenx not leny')
        expr = exprr(x,y)
        #print('expr[0] :(corr)',expr[0].shape)
        if coefs is None:
            Amatrix = np.array(expr).T
            return Amatrix
        else:
            Aresult = np.sum(np.array(expr).T*coefs,axis=1)
            return Aresult

    coeff2, r2, rank2, s2 = np.linalg.lstsq(expression(x,y,expuse), z,rcond=None)
    print('----------'*3)
    print('len(coeff2)',len(coeff2))
    print('r2         ',r2)
    #print('rank2      ',rank2)
    #print('s2         ',s2)

    denset = np.arange(0,x.max()+1,1)
    all1v = np.unique(data[:,1])
    for v in all1v:
        out2dense = expression(denset, np.repeat(v,len(denset)), expuse, coeff2)
        np.savetxt("out_OUT2dense_"+str(v)+'dat.dat',out2dense)

    # get diffs
    def get_diffmax(data,expuse,coeff2):
        diffmax = 0
        for i in data:
            #print('i[0]',i[0],i[1],expuse(i[0],i[1]))
            fah_fit = expression(np.array([i[0]]), np.array([i[1]]), expuse, coeff2)[0]
            diff = np.abs(fah_fit-i[2])
            if diff > diffmax:
                diffmax = diff
            #print(i[0],i[1],'fah_fit',fah_fit,i[2],fah_fit-i[2])
        return diffmax

    diffmax = get_diffmax(data,expuse,coeff2)
    print('diffmax',diffmax)
    print('----------'*3)


    #def expression(x,y,exprr,coefs=None):

    def save_v_lines_from2darray(usearray=None,string=None):
        if usearray is None:
            sys.exit('provide usearray1 and string')
        if string is None:
            sys.exit('provide usearray2 and string')
        if type(string) != str:
            sys.exit('has tobe a sring')
        all1t = np.unique(usearray[:,0])
        all1v = np.unique(usearray[:,1])
        for v in all1v:
            ka = usearray[np.where(usearray[:,1]==v)[0]][:,[0,2]]
            arr = ka[ka[:,0].argsort()]
            np.savetxt("out_"+string+"_"+str(v)+'dat.dat',arr)
        return

    save_v_lines_from2darray(usearray=data,string='data')
    #save_v_lines_from2darray(usearray=OUT1,string='OUT1')
    #save_v_lines_from2darray(usearray=OUT2,string='OUT2')
    #X,Y,Z = myp.xyz_line_to_xyz_unique(x=x,y=y,z=z,data=None)
    #print('X')
    #print(X)
    #print(X[0])
    #print('Y')
    #print(Y)
    #print(Y[:,0])
    #kx=3
    #ky=3
    #soln,ka,kb,kc = polyfit2dbetter(X, Y, Z, kx=kx, ky=ky, order=None)
    #fitted_surf = np.polynomial.polynomial.polygrid2d(X, Y, soln.reshape((kx+1,ky+1)))
    #print('xxxcoefs',soln)
    #print('xxxcoefs',fitted_surf)
    return

if __name__ == "__main__":
    np.set_printoptions(precision=4)
    filepath = "/Users/glensk/Downloads/cu/fah/thermo/Fah_surface"
    data = np.loadtxt(filepath)[:,[0,1,2]]
    print(filepath)
    print(os.path.basename(filepath))
    print(os.path.dirname(filepath))
    sys.exit()
    x = data[:,0];
    y = data[:,1];
    z = data[:,2];
    lmfit2d(x=x,y=y,z=z)
    sys.exit('6789')
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

