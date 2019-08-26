#!/usr/bin/env python
 # -*- coding: utf-8 -*-
from __future__ import print_function
import numpy as np
import os,sys,argparse

def _gen_strainmatrix(e1=0.,e2=0.,e3=0.,e4=0.,e5=0.,e6=0.):
    e_M = np.array([
        [e1,     e6/2.,  e5/2.],
        [e6/2.,  e2,     e4/2.],
        [e5/2.,  e4/2.,  e3   ],
    ])
    return e_M

def _apply_strain(x_M, e_M):
    I_M = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    r_M = np.matmul(x_M, I_M+e_M)
    return r_M

def _find_angle(u, v):
    """
    Returns the angle between vectors (numpy arrays) u and v
    """
    n_u = np.linalg.norm(u)
    n_v = np.linalg.norm(v)
    ang_uv = np.arccos( np.clip(np.dot(u, v)/(n_u*n_v), -1, 1))
    return ang_uv

def _2lammpslattice(b_M):
    """
    Converts a basis array, b_M to a lammps-friendly lattice definition
    i.e. one that uses xlo xhi ylo yhi zlo zhi xy xz yz
    """
    n_x = np.linalg.norm(b_M[0])
    n_y = np.linalg.norm(b_M[1])
    n_z = np.linalg.norm(b_M[2])

    ang_yz = _find_angle(b_M[1], b_M[2]) #alpha
    ang_xz = _find_angle(b_M[0], b_M[2]) #beta
    ang_xy = _find_angle(b_M[0], b_M[1]) #gamma

    xhi = n_x
    xy = n_y * np.cos(ang_xy)
    xz = n_z * np.cos(ang_xz)
    yhi = np.sqrt(n_y ** 2. - xy ** 2.)
    yz = (n_y*n_z*np.cos(ang_yz)-xy*xz)/yhi
    zhi = np.sqrt(n_z ** 2. - xz ** 2. - yz ** 2.)

    l_M = np.array([
        [xhi, 0., 0.],
        [xy, yhi, 0.],
        [xz, yz, zhi],
    ])

    #l_M = np.round(l_M, decimals=15)
    return l_M
