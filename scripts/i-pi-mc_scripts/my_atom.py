#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
import numpy as np
import sys
import math


np.set_printoptions(suppress=True)   # display arrays withou 000000
np.set_printoptions(precision=6)    # print only 6 digist after .

atomic_symbols = [
    '',  'H',  'He', 'Li', 'Be',
    'B',  'C',  'N',  'O',  'F',
    'Ne', 'Na', 'Mg', 'Al', 'Si',
    'P',  'S',  'Cl', 'Ar', 'K',
    'Ca', 'Sc', 'Ti', 'V',  'Cr',
    'Mn', 'Fe', 'Co', 'Ni', 'Cu',
    'Zn', 'Ga', 'Ge', 'As', 'Se',
    'Br', 'Kr', 'Rb', 'Sr', 'Y',
    'Zr', 'Nb', 'Mo', 'Tc', 'Ru',
    'Rh', 'Pd', 'Ag', 'Cd', 'In',
    'Sn', 'Sb', 'Te', 'I',  'Xe',
    'Cs', 'Ba', 'La', 'Ce', 'Pr',
    'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
    'Tb', 'Dy', 'Ho', 'Er', 'Tm',
    'Yb', 'Lu', 'Hf', 'Ta', 'W',
    'Re', 'Os', 'Ir', 'Pt', 'Au',
    'Hg', 'Tl', 'Pb', 'Bi', 'Po',
    'At', 'Rn', 'Fr', 'Ra', 'Ac',
    'Th', 'Pa', 'U',  'Np', 'Pu',
    'Am', 'Cm', 'Bk', 'Cf', 'Es',
    'Fm', 'Md', 'No', 'Lr']

atomic_z = np.arange(len(atomic_symbols))

atomic_names = [
    '', 'Hydrogen', 'Helium', 'Lithium', 'Beryllium', 'Boron',
    'Carbon', 'Nitrogen', 'Oxygen', 'Fluorine', 'Neon', 'Sodium',
    'Magnesium', 'Aluminium', 'Silicon', 'Phosphorus', 'Sulfur',
    'Chlorine', 'Argon', 'Potassium', 'Calcium', 'Scandium',
    'Titanium', 'Vanadium', 'Chromium', 'Manganese', 'Iron',
    'Cobalt', 'Nickel', 'Copper', 'Zinc', 'Gallium', 'Germanium',
    'Arsenic', 'Selenium', 'Bromine', 'Krypton', 'Rubidium',
    'Strontium', 'Yttrium', 'Zirconium', 'Niobium', 'Molybdenum',
    'Technetium', 'Ruthenium', 'Rhodium', 'Palladium', 'Silver',
    'Cadmium', 'Indium', 'Tin', 'Antimony', 'Tellurium',
    'Iodine', 'Xenon', 'Caesium', 'Barium', 'Lanthanum',
    'Cerium', 'Praseodymium', 'Neodymium', 'Promethium',
    'Samarium', 'Europium', 'Gadolinium', 'Terbium',
    'Dysprosium', 'Holmium', 'Erbium', 'Thulium', 'Ytterbium',
    'Lutetium', 'Hafnium', 'Tantalum', 'Tungsten', 'Rhenium',
    'Osmium', 'Iridium', 'Platinum', 'Gold', 'Mercury',
    'Thallium', 'Lead', 'Bismuth', 'Polonium', 'Astatine',
    'Radon', 'Francium', 'Radium', 'Actinium', 'Thorium',
    'Protactinium', 'Uranium', 'Neptunium', 'Plutonium',
    'Americium', 'Curium', 'Berkelium', 'Californium',
    'Einsteinium', 'Fermium', 'Mendelevium', 'Nobelium',
    'Lawrencium']

atomic_masses_ase = np.array([
       0.00000, # X
       1.00794, # H
       4.00260, # He
       6.94100, # Li
       9.01218, # Be
      10.81100, # B
      12.01100, # C
      14.00670, # N
      15.99940, # O
      18.99840, # F
      20.17970, # Ne
      22.98977, # Na
      24.30500, # Mg
      26.98154, # Al
      28.08550, # Si
      30.97376, # P
      32.06600, # S
      35.45270, # Cl
      39.94800, # Ar
      39.09830, # K
      40.07800, # Ca
      44.95590, # Sc
      47.88000, # Ti
      50.94150, # V
      51.99600, # Cr
      54.93800, # Mn
      55.84700, # Fe
      58.93320, # Co
      58.69340, # Ni
      63.54600, # Cu
      65.39000, # Zn
      69.72300, # Ga
      72.61000, # Ge
      74.92160, # As
      78.96000, # Se
      79.90400, # Br
      83.80000, # Kr
      85.46780, # Rb
      87.62000, # Sr
      88.90590, # Y
      91.22400, # Zr
      92.90640, # Nb
      95.94000, # Mo
        np.nan, # Tc
     101.07000, # Ru
     102.90550, # Rh
     106.42000, # Pd
     107.86800, # Ag
     112.41000, # Cd
     114.82000, # In
     118.71000, # Sn
     121.75700, # Sb
     127.60000, # Te
     126.90450, # I
     131.29000, # Xe
     132.90540, # Cs
     137.33000, # Ba
     138.90550, # La
     140.12000, # Ce
     140.90770, # Pr
     144.24000, # Nd
        np.nan, # Pm
     150.36000, # Sm
     151.96500, # Eu
     157.25000, # Gd
     158.92530, # Tb
     162.50000, # Dy
     164.93030, # Ho
     167.26000, # Er
     168.93420, # Tm
     173.04000, # Yb
     174.96700, # Lu
     178.49000, # Hf
     180.94790, # Ta
     183.85000, # W
     186.20700, # Re
     190.20000, # Os
     192.22000, # Ir
     195.08000, # Pt
     196.96650, # Au
     200.59000, # Hg
     204.38300, # Tl
     207.20000, # Pb
     208.98040, # Bi
     np.nan, # Po
     np.nan, # At
     np.nan, # Rn
     np.nan, # Fr
     226.02540, # Ra
     np.nan, # Ac
     232.03810, # Th
     231.03590, # Pa
     238.02900, # U
     237.04820, # Np
     np.nan, # Pu
     np.nan, # Am
     np.nan, # Cm
     np.nan, # Bk
     np.nan, # Cf
     np.nan, # Es
     np.nan, # Fm
     np.nan, # Md
     np.nan, # No
     np.nan])# Lw

# from http://en.wikipedia.org/wiki/List_of_elements (see there for corresponding REFS)
atomic_masses = np.array([
    np.nan     ,    1.008     ,    4.002602  ,    6.94      ,
    9.012182   ,   10.81      ,   12.011     ,   14.007     ,
    15.999     ,   18.9984032 ,   20.1797    ,   22.98976928,
    24.305     ,   26.9815386 ,   28.085     ,   30.973762  ,
    32.06      ,   35.45      ,   39.948     ,   39.0983    ,
    40.078     ,   44.955912  ,   47.867     ,   50.9415    ,
    51.9961    ,   54.938045  ,   55.845     ,   58.933195  ,
    58.6934    ,   63.546     ,   65.38      ,   69.723     ,
    72.63      ,   74.9216    ,   78.96      ,   79.904     ,
    83.798     ,   85.4678    ,   87.62      ,   88.90585   ,
    91.224     ,   92.90638   ,   95.96      ,        np.nan,
    101.07     ,  102.9055    ,  106.42      ,  107.8682    ,
    112.411    ,  114.818     ,  118.71      ,  121.76      ,
    127.6      ,  126.90447   ,  131.293     ,  132.9054519 ,
    137.327    ,  138.90547   ,  140.116     ,  140.90765   ,
    144.242    ,        np.nan,  150.36      ,  151.964     ,
    157.25     ,  158.92535   ,  162.5       ,  164.93032   ,
    167.259    ,  168.93421   ,  173.054     ,  174.9668    ,
    178.49     ,  180.94788   ,  183.84      ,  186.207     ,
    190.23     ,  192.217     ,  195.084     ,  196.966569  ,
    200.592    ,  204.38      ,  207.2       ,  208.9804    ,
    np.nan,             np.nan,        np.nan,        np.nan,
    np.nan,             np.nan,  232.03806   ,  231.03588   ,
    238.02891  ,        np.nan,        np.nan,        np.nan,
    np.nan,             np.nan,        np.nan,        np.nan,
    np.nan,             np.nan,        np.nan,        np.nan])

atomic_masses_potcar = np.array([
         np.nan,    1.   ,    4.   ,    7.01 ,    9.013,   10.811,
         12.011,   14.001,   16.   ,   18.998,   20.18 ,   22.99 ,
         24.305,   26.982,   28.085,   30.974,   32.066,   35.453,
         39.949,      np.nan,      np.nan,      np.nan,   47.88 ,   50.941,
         51.996,   54.938,   55.847,   58.933,   58.69 ,   63.546,
         65.39 ,   69.723,   72.61 ,   74.922,   78.96 ,   79.904,
         83.8  ,      np.nan,      np.nan,      np.nan,   91.224,      np.nan,
         95.94 ,   98.906,  101.07 ,  102.906,  106.42 ,  107.868,
        112.411,  114.82 ,  118.71 ,  121.75 ,  127.6  ,  126.904,
        131.294,      np.nan,      np.nan,  138.9  ,  140.115,  140.907,
        144.24 ,  146.915,  150.36 ,  151.965,  157.25 ,      np.nan,
            np.nan,      np.nan,      np.nan,  168.93 ,  173.04 ,  174.967,
        178.49 ,  180.948,  183.85 ,  186.207,  190.2  ,  192.22 ,
        195.08 ,  196.966,  200.59 ,  204.38 ,  207.2  ,  208.98 ,
            np.nan,      np.nan,      np.nan,      np.nan,      np.nan,  227.028,
        232.039,  231.036,  238.029,  237.048,  244.064,      np.nan,
            np.nan,      np.nan,      np.nan,      np.nan,      np.nan,      np.nan,
            np.nan,      np.nan])
# Covalent radii from:
#
#  Covalent radii revisited,
#  Beatriz Cordero, Verónica Gómez, Ana E. Platero-Prats, Marc Revés,
#  Jorge Echeverría, Eduard Cremades, Flavia Barragán and Santiago Alvarez,
#  Dalton Trans., 2008, 2832-2838 DOI:10.1039/B801115J
_missing = 0.2
atomic_covalent_radii = np.array([
    _missing,  # X
    0.31,  # H
    0.28,  # He
    1.28,  # Li
    0.96,  # Be
    0.84,  # B
    0.76,  # C
    0.71,  # N
    0.66,  # O
    0.57,  # F
    0.58,  # Ne
    1.66,  # Na
    1.41,  # Mg
    1.21,  # Al
    1.11,  # Si
    1.07,  # P
    1.05,  # S
    1.02,  # Cl
    1.06,  # Ar
    2.03,  # K
    1.76,  # Ca
    1.70,  # Sc
    1.60,  # Ti
    1.53,  # V
    1.39,  # Cr
    1.39,  # Mn
    1.32,  # Fe
    1.26,  # Co
    1.24,  # Ni
    1.32,  # Cu
    1.22,  # Zn
    1.22,  # Ga
    1.20,  # Ge
    1.19,  # As
    1.20,  # Se
    1.20,  # Br
    1.16,  # Kr
    2.20,  # Rb
    1.95,  # Sr
    1.90,  # Y
    1.75,  # Zr
    1.64,  # Nb
    1.54,  # Mo
    1.47,  # Tc
    1.46,  # Ru
    1.42,  # Rh
    1.39,  # Pd
    1.45,  # Ag
    1.44,  # Cd
    1.42,  # In
    1.39,  # Sn
    1.39,  # Sb
    1.38,  # Te
    1.39,  # I
    1.40,  # Xe
    2.44,  # Cs
    2.15,  # Ba
    2.07,  # La
    2.04,  # Ce
    2.03,  # Pr
    2.01,  # Nd
    1.99,  # Pm
    1.98,  # Sm
    1.98,  # Eu
    1.96,  # Gd
    1.94,  # Tb
    1.92,  # Dy
    1.92,  # Ho
    1.89,  # Er
    1.90,  # Tm
    1.87,  # Yb
    1.87,  # Lu
    1.75,  # Hf
    1.70,  # Ta
    1.62,  # W
    1.51,  # Re
    1.44,  # Os
    1.41,  # Ir
    1.36,  # Pt
    1.36,  # Au
    1.32,  # Hg
    1.45,  # Tl
    1.46,  # Pb
    1.48,  # Bi
    1.40,  # Po
    1.50,  # At
    1.50,  # Rn
    2.60,  # Fr
    2.21,  # Ra
    2.15,  # Ac
    2.06,  # Th
    2.00,  # Pa
    1.96,  # U
    1.90,  # Np
    1.87,  # Pu
    1.80,  # Am
    1.69,  # Cm
    _missing,  # Bk
    _missing,  # Cf
    _missing,  # Es
    _missing,  # Fm
    _missing,  # Md
    _missing,  # No
    _missing,  # Lr
    ])

# This data is from Ashcroft and Mermin.
atomic_reference_states = [\
    None, #X
    {'symmetry': 'diatom', 'd': 0.74}, #H
    {'symmetry': 'atom'}, #He
    {'symmetry': 'bcc', 'a': 3.49}, #Li
    {'symmetry': 'hcp', 'c/a': 1.567, 'a': 2.29}, #Be
    {'symmetry': 'tetragonal', 'c/a': 0.576, 'a': 8.73}, #B
    {'symmetry': 'diamond', 'a': 3.57},#C
    {'symmetry': 'diatom', 'd': 1.10},#N
    {'symmetry': 'diatom', 'd': 1.21},#O
    {'symmetry': 'diatom', 'd': 1.42},#F
    {'symmetry': 'fcc', 'a': 4.43},#Ne
    {'symmetry': 'bcc', 'a': 4.23},#Na
    {'symmetry': 'hcp', 'c/a': 1.624, 'a': 3.21},#Mg
    {'symmetry': 'fcc', 'a': 4.05},#Al
    {'symmetry': 'diamond', 'a': 5.43},#Si
    {'symmetry': 'cubic', 'a': 7.17},#P
    {'symmetry': 'orthorhombic', 'c/a': 2.339, 'a': 10.47,'b/a': 1.229},#S
    {'symmetry': 'orthorhombic', 'c/a': 1.324, 'a': 6.24, 'b/a': 0.718},#Cl
    {'symmetry': 'fcc', 'a': 5.26},#Ar
    {'symmetry': 'bcc', 'a': 5.23},#K
    {'symmetry': 'fcc', 'a': 5.58},#Ca
    {'symmetry': 'hcp', 'c/a': 1.594, 'a': 3.31},#Sc
    {'symmetry': 'hcp', 'c/a': 1.588, 'a': 2.95},#Ti
    {'symmetry': 'bcc', 'a': 3.02},#V
    {'symmetry': 'bcc', 'a': 2.88},#Cr
    {'symmetry': 'cubic', 'a': 8.89},#Mn
    {'symmetry': 'bcc', 'a': 2.87},#Fe
    {'symmetry': 'hcp', 'c/a': 1.622, 'a': 2.51},#Co
    {'symmetry': 'fcc', 'a': 3.52},#Ni
    {'symmetry': 'fcc', 'a': 3.61},#Cu
    {'symmetry': 'hcp', 'c/a': 1.856, 'a': 2.66},#Zn
    {'symmetry': 'orthorhombic', 'c/a': 1.695, 'a': 4.51, 'b/a': 1.001},#Ga
    {'symmetry': 'diamond', 'a': 5.66},#Ge
    {'symmetry': 'rhombohedral', 'a': 4.13, 'alpha': 54.10},#As
    {'symmetry': 'hcp', 'c/a': 1.136, 'a': 4.36},#Se
    {'symmetry': 'orthorhombic', 'c/a': 1.307, 'a': 6.67, 'b/a': 0.672},#Br
    {'symmetry': 'fcc', 'a': 5.72},#Kr
    {'symmetry': 'bcc', 'a': 5.59},#Rb
    {'symmetry': 'fcc', 'a': 6.08},#Sr
    {'symmetry': 'hcp', 'c/a': 1.571, 'a': 3.65},#Y
    {'symmetry': 'hcp', 'c/a': 1.593, 'a': 3.23},#Zr
    {'symmetry': 'bcc', 'a': 3.30},#Nb
    {'symmetry': 'bcc', 'a': 3.15},#Mo
    {'symmetry': 'hcp', 'c/a': 1.604, 'a': 2.74},#Tc
    {'symmetry': 'hcp', 'c/a': 1.584, 'a': 2.70},#Ru
    {'symmetry': 'fcc', 'a': 3.80},#Rh
    {'symmetry': 'fcc', 'a': 3.89},#Pd
    {'symmetry': 'fcc', 'a': 4.09},#Ag
    {'symmetry': 'hcp', 'c/a': 1.886, 'a': 2.98},#Cd
    {'symmetry': 'tetragonal', 'c/a': 1.076, 'a': 4.59},#In
    {'symmetry': 'tetragonal', 'c/a': 0.546, 'a': 5.82},#Sn
    {'symmetry': 'rhombohedral', 'a': 4.51, 'alpha': 57.60},#Sb
    {'symmetry': 'hcp', 'c/a': 1.330, 'a': 4.45},#Te
    {'symmetry': 'orthorhombic', 'c/a': 1.347, 'a': 7.27, 'b/a': 0.659},#I
    {'symmetry': 'fcc', 'a': 6.20},#Xe
    {'symmetry': 'bcc', 'a': 6.05},#Cs
    {'symmetry': 'bcc', 'a': 5.02},#Ba
    {'symmetry': 'hcp', 'c/a': 1.619, 'a': 3.75},#La
    {'symmetry': 'fcc', 'a': 5.16},#Ce
    {'symmetry': 'hcp', 'c/a': 1.614, 'a': 3.67},#Pr
    {'symmetry': 'hcp', 'c/a': 1.614, 'a': 3.66},#Nd
    None,#Pm
    {'symmetry': 'rhombohedral', 'a': 9.00, 'alpha': 23.13},#Sm
    {'symmetry': 'bcc', 'a': 4.61},#Eu
    {'symmetry': 'hcp', 'c/a': 1.588, 'a': 3.64},#Gd
    {'symmetry': 'hcp', 'c/a': 1.581, 'a': 3.60},#Th
    {'symmetry': 'hcp', 'c/a': 1.573, 'a': 3.59},#Dy
    {'symmetry': 'hcp', 'c/a': 1.570, 'a': 3.58},#Ho
    {'symmetry': 'hcp', 'c/a': 1.570, 'a': 3.56},#Er
    {'symmetry': 'hcp', 'c/a': 1.570, 'a': 3.54},#Tm
    {'symmetry': 'fcc', 'a': 5.49},#Yb
    {'symmetry': 'hcp', 'c/a': 1.585, 'a': 3.51},#Lu
    {'symmetry': 'hcp', 'c/a': 1.582, 'a': 3.20},#Hf
    {'symmetry': 'bcc', 'a': 3.31},#Ta
    {'symmetry': 'bcc', 'a': 3.16},#W
    {'symmetry': 'hcp', 'c/a': 1.615, 'a': 2.76},#Re
    {'symmetry': 'hcp', 'c/a': 1.579, 'a': 2.74},#Os
    {'symmetry': 'fcc', 'a': 3.84},#Ir
    {'symmetry': 'fcc', 'a': 3.92},#Pt
    {'symmetry': 'fcc', 'a': 4.08},#Au
    {'symmetry': 'rhombohedral', 'a': 2.99, 'alpha': 70.45},#Hg
    {'symmetry': 'hcp', 'c/a': 1.599, 'a': 3.46},#Tl
    {'symmetry': 'fcc', 'a': 4.95},#Pb
    {'symmetry': 'rhombohedral', 'a': 4.75, 'alpha': 57.14},#Bi
    {'symmetry': 'sc', 'a': 3.35},#Po
    None,#At
    None,#Rn
    None,#Fr
    None,#Ra
    {'symmetry': 'fcc', 'a': 5.31},#Ac
    {'symmetry': 'fcc', 'a': 5.08},#Th
    {'symmetry': 'tetragonal', 'c/a': 0.825, 'a': 3.92},#Pa
    {'symmetry': 'orthorhombic', 'c/a': 2.056, 'a': 2.85, 'b/a': 1.736},#U
    {'symmetry': 'orthorhombic', 'c/a': 1.411, 'a': 4.72, 'b/a': 1.035},#Np
    {'symmetry': 'monoclinic'},#Pu
    None,#Am
    None,#Cm
    None,#Bk
    None,#Cf
    None,#Es
    None,#Fm
    None,#Md
    None,#No
    None]#Lw

# http://www.webelements.com
atomic_ground_state_magnetic_moments = np.array([
    0.0, # X
    1.0, # H
    0.0, # He
    1.0, # Li
    0.0, # Be
    1.0, # B
    2.0, # C
    3.0, # N
    2.0, # O
    1.0, # F
    0.0, # Ne
    1.0, # Na
    0.0, # Mg
    1.0, # Al
    2.0, # Si
    3.0, # P
    2.0, # S
    1.0, # Cl
    0.0, # Ar
    1.0, # K
    0.0, # Ca
    1.0, # Sc
    2.0, # Ti
    3.0, # V
    6.0, # Cr
    5.0, # Mn
    4.0, # Fe
    3.0, # Co
    2.0, # Ni
    1.0, # Cu
    0.0, # Zn
    1.0, # Ga
    2.0, # Ge
    3.0, # As
    2.0, # Se
    1.0, # Br
    0.0, # Kr
    1.0, # Rb
    0.0, # Sr
    1.0, # Y
    2.0, # Zr
    5.0, # Nb
    6.0, # Mo
    5.0, # Tc
    4.0, # Ru
    3.0, # Rh
    0.0, # Pd
    1.0, # Ag
    0.0, # Cd
    1.0, # In
    2.0, # Sn
    3.0, # Sb
    2.0, # Te
    1.0, # I
    0.0, # Xe
    1.0, # Cs
    0.0, # Ba
    1.0, # La
    1.0, # Ce
    3.0, # Pr
    4.0, # Nd
    5.0, # Pm
    6.0, # Sm
    7.0, # Eu
    8.0, # Gd
    5.0, # Tb
    4.0, # Dy
    3.0, # Ho
    2.0, # Er
    1.0, # Tm
    0.0, # Yb
    1.0, # Lu
    2.0, # Hf
    3.0, # Ta
    4.0, # W
    5.0, # Re
    4.0, # Os
    3.0, # Ir
    2.0, # Pt
    1.0, # Au
    0.0, # Hg
    1.0, # Tl
    2.0, # Pb
    3.0, # Bi
    2.0, # Po
    1.0, # At
    0.0, # Rn
    1.0, # Fr
    0.0, # Ra
    1.0, # Ac
    2.0, # Th
    3.0, # Pa
    4.0, # U
    5.0, # Np
    6.0, # Pu
    7.0, # Am
    8.0, # Cm
    5.0, # Bk
    4.0, # Cf
    4.0, # Es
    2.0, # Fm
    1.0, # Md
    0.0, # No
    np.nan])# Lw

atomic_melting_webelements = np.array([
        np.nan,    14.01,      0.95,   453.69,  1560.  ,  2349.  ,
        3800.,    63.05,    54.8 ,    53.53,      24.56,   370.87,
        923.  ,   933.47,  1687.  ,      883.,      388.36,   171.6 ,
        83.8,   336.53,  1115.  ,  1814.  ,  1941.  ,  2183.  ,
        2180.  ,  1519.  ,  1811.  ,  1768.  ,  1728.  ,  1357.77,
        692.68,   302.91,  1211.4 ,  1090.  ,   494.  ,   265.8 ,
        115.79,   312.46,  1050.  ,  1799.  ,  2128.  ,  2750.  ,
        2896.  ,  2430.  ,  2607.  ,  2237.  ,  1828.05,  1234.93,
        594.22,   429.75,   505.08,   903.78,   722.66,   386.85,
        161.4 ,   301.59,  1000.  ,  1193.  ,  1068.  ,  1208.  ,
        1297.  ,  1373.  ,  1345.  ,  1099.  ,  1585.  ,  1629.  , 1680.
        ,  1734.  ,  1770.  ,  1818.  ,  1097.  ,  1925.  , 2506.
        ,  3290.  ,  3695.  ,  3459.  ,  3306.  ,  2739.  , 2041.4 ,
        1337.33,   234.32,   577.  ,   600.61,   544.4 ,
         527.  ,   575.  ,   202.  ,      300.,   973.  ,  1323.  ,
        2020.  ,  1843.  ,  1407.  ,   913.  ,   913.  ,  1449.  , 1613.
        ,      1323.,  1173.  ,  1133.  ,      1800.,      1100.,
            1100.,      1900.])

atomic_melting = atomic_melting_webelements
#['al, '1', 'cu', '3' ] --> ['al', 'cu', 'cu', 'cu']
def get_elements(listin, verbose = False, get_elements_count = False):
    '''
    ['al, '1', 'cu', '3' ] --> ['Al', 'Cu']
    '''
    import sys

    def is_int(x):
        try:
            a = float(x)
            b = int(a)
        except ValueError:
            return False
        else:
            return a == b


    elements = []
    elements_count = []
    elements_count_str = []
    for i in listin:
        if is_int(i) == True:
            if verbose:
                print(i,"yes")
            elements_count.append(int(i))
            elements_count_str.append(i)
        else:
            if verbose:
                print(i,"no")
            elements.append(i)
    if verbose:
        print("elements      :",elements,len(elements))
        print("elements_count:",elements_count,len(elements_count))

    if len(elements_count) == 0:
        elements_count = [1]*len(elements)

    if len(elements) != len(elements_count):
        m1 = "Length of elements and numbers has to be equal. Examples: \n"
        m2 = "GOOD: si ag 1 31;    GOOD :  si 1 al 31\n"
        m3 = "NOT GOOD: al si 31  ; NOT GOOD: al 7 si 6 cu\n"
        m4 = "You have "+str(len(elements))+" elements ("+" ".join(elements)+") "
        m5 = "but only "+str(len(elements_count))+" numbers. ("
        m6 = " ".join(elements_count_str)+")\n"
        sys.exit(m4+m5+m6+m1+m2+m3)

    # convert elements to spropper uppercase elements and check if available
    all_lower = [a.lower() for a in elements]
    symbol = [word[0].upper() + word[1:] for word in all_lower]

    # check if element is available
    number = [atomic_symbols.index(a) for a in symbol]  # better, this gives traceback in error in ipython
    #try:
    #    number = [atomic_symbols.index(a) for a in symbol]
    #except ValueError as e:
    #    #print "I/o:",e
    #    sys.exit("Atomic species:"+str(e)+" of valid elements:"+str(atomic_symbols))

    if get_elements_count == True:
        return symbol, elements_count
    else:
        return symbol
def get_element_from_path(path):
    ''' gets a path and returns, if not ambigious, the element symbol
    the element is expected to be the symbol somewhere in the path like .../Al/...'''
    elements_out = []
    for i in path.split("/"):
        #print i,len(i)
        if len(i) <= 2 and len(i) >= 1: # element symbols are 1 or 2 chars long
            elements_out.append(i)
            #print "ele:",i
    #print elements_out,len(elements_out)
    if len(elements_out) == 1:
        return elements_out[0]
    else:
        print("element(s):",elements_out)
        return False





def get_elements_list(listin, verbose = False):
    ''' ['al, '1', 'cu', '3' ] --> ['al', 'cu', 'cu', 'cu'] '''
    symbol, elements_count = get_elements(
            listin = listin,
            verbose = verbose,
            get_elements_count = True)
    elements_list = []
    for i in range(len(elements_count)):
        elements_list.extend([symbol[i]]*elements_count[i])

    if verbose:
        print("||elements    ",elements)
        print("||elements_list",elements_list)

    return elements_list

def data(elementlistin = None):
    ''' returnes all data of elements given '''
    if elementlistin == None:
            sys.exit("please specify at least one element")

    elementlist = get_elements_list(elementlistin)

    if len(elementlist) == 0:
        sys.exit("please specify at least one element")

    length = [len(a) for a in elementlist]
    if max(length) > 2:
        sys.exit("one or more from elementlist: "+str(_elementlist)+" is longer than 2 cars; unknown")

    all_lower = [a.lower() for a in elementlist]


    symbol = [word[0].upper() + word[1:] for word in all_lower]
    number = np.array([atomic_symbols.index(a) for a in symbol])
    mass = atomic_masses[number]
    melting = atomic_melting[number]

    getVar = lambda searchList, ind: [searchList[i] for i in ind]
    name = getVar(atomic_names, number)
    reference_state = getVar(atomic_reference_states, number)
    return number,symbol,name,mass,melting,reference_state

def number(listin):
    return data(listin)[0]

def symbol(listin):
    return data(listin)[1]

def name(listin):
    return data(listin)[2]

def mass(listin):
    return data(listin)[3]

def melting(listin):
    return data(listin)[4]

def reference_state(listin):
    return data(listin)[5]

def calphad_free_energy(listin):
    pass

def mass_potcar(listin):
    pass

class plots(object):
    def __init__(self):
        self.hallo = 1

    def __call__(self):
        #atom = atom() # we need atom.atomic_z, ....
        #print "az:",atom.atomic_z
        self.plot_mass_difference()
        return

    def plot_mass_difference(self):
        import utils_plot
        #utils_plot.plot_string_at_point(atom.atomic_z,atom.atomic_masses_ase - atom.atomic_masses_potcar,atom.atomic_symbols)
        #utils_plot.plot_string_at_point(atom.atomic_z,atom.atomic_masses - atom.atomic_masses_potcar,atom.atomic_symbols)

        f1 = utils_plot.plot_string_at_point(atomic_z,atomic_masses_ase - atomic_masses_potcar,atomic_symbols)
        f2 = utils_plot.plot_string_at_point(atomic_z,atomic_masses - atomic_masses_potcar,atomic_symbols)
        f3 = utils_plot.plot_string_at_point(atomic_z,atomic_masses_ase - atomic_masses_potcar,atomic_symbols)


class atom(object):     # advantage: once this is executed we know that everything is defined
                        # im __init__ sollte nie ein systemabbrch erfolgen (keine checks)
    def __init__(self, listin = None):  # elements: [aL cU 1 3], [SI], [Ag 7]
        # load in other script with:
            # import atom
            # ka = atom.atom([element])
            # print ka.melting_rounded[0]
            #
            # or
            # import my_atom
            # atom = my_atom.atom([args.element])
            # print atom.melting_rounded[0]
            #
            #
        self._z         = False
        self._n         = False
        self._s         = False
        self._m         = False
        self._r         = False
        self._t         = False
        self._tr        = False
        self._a         = False
        self._str       = False
        self._verbose   = False
        self._listin    = listin    # self._listin: [aL cU 1 3], [SI], [Ag 7]

        # define self_listin
        if self._listin == None:
            #help (may generally only change __init__ optons when calling from system shell)
            self._help() # defines self._listin

        if self._verbose:
            print("self._listin:",self._listin)

        elementlist = get_elements_list(self._listin, verbose = self._verbose)
        self.number = number(elementlist)
        self.symbol= symbol(elementlist)
        self.name = name(elementlist)
        self.mass = mass(elementlist)
        self.melting = melting(elementlist)
        self.melting_rounded = [ int(math.ceil(b)) for b in self.melting ]
        self.reference_state = reference_state(elementlist)

        #print "self._str:",self._str
        def out_as_str(a):
            #a=data(listin)[3]
            s=""
            for i in a:
                s+=str(i)+"  "
            return s

        if self._str == False:
            if self._z  :print(self.number)
            if self._n  :print(self.name)
            if self._s  :print(self.symbol)
            if self._m  :print(self.mass)
            if self._r  :print(self.reference_state)
            if self._t  :print(self.melting)
            if self._tr :print(self.melting_rounded)
            if self._a  :
                print(self.number)
                print(self.name)
                print(self.symbol)
                print(self.mass)
                print(self.reference_state)
                print(self.melting)
                print(self.melting_rounded)
        if self._str == True:
            if self._z  :print(out_as_str(self.number))
            if self._n  :print(out_as_str(self.name))
            if self._s  :print(out_as_str(self.symbol))
            if self._m  :print(out_as_str(self.mass))
            if self._r  :print(out_as_str(self.reference_state))
            if self._t  :print(out_as_str(self.melting))
            if self._tr :print(out_as_str(self.melting_rounded))
            if self._a  :
                print(out_as_str(self.number))
                print(out_as_str(self.name))
                print(out_as_str(self.symbol))
                print(out_as_str(self.mass))
                print(out_as_str(self.reference_state))
                print(out_as_str(self.melting))
                print(out_as_str(self.melting_rounded))
        return


    def __call__(self):
        return self.symbol

    def _help(self):
        import argparse
        #from argparse import ArgumentDefaultsHelpFormatter
        from argparse import RawTextHelpFormatter
        string = ''' element database;

        References (masses):
            Royal Society of Chemistry – Visual Element Periodic Table
            data taken from http://en.wikipedia.org/wiki/List_of_elements

        References (melting points):
            Reference1: http://www.ptable.com (which is from wiki)
            Reference2: J. Phys. Chem. Ref. Data, Vol. 39, No. 4, 2010 for actinides

        References (reference state):
            Ashcroft and Mermin

        Available elements:
        Z : 'Element'
        1 : 'H'  'He' 'Li' 'Be' 'B'  'C'  'N'  'O'  'F'  'Ne' 'Na' 'Mg' 'Al'
        14: 'Si' 'P'  'S'  'Cl' 'Ar' 'K'  'Ca' 'Sc' 'Ti' 'V'  'Cr' 'Mn' 'Fe'
        27: 'Co' 'Ni' 'Cu' 'Zn' 'Ga' 'Ge' 'As' 'Se' 'Br' 'Kr' 'Rb' 'Sr' 'Y'
        40: 'Zr' 'Nb' 'Mo' 'Tc' 'Ru' 'Rh' 'Pd' 'Ag' 'Cd' 'In' 'Sn' 'Sb' 'Te'
        53: 'I'  'Xe' 'Cs' 'Ba' 'La' 'Ce' 'Pr' 'Nd' 'Pm' 'Sm' 'Eu' 'Gd' 'Tb'
        66: 'Dy' 'Ho' 'Er' 'Tm' 'Yb' 'Lu' 'Hf' 'Ta' 'W'  'Re' 'Os' 'Ir' 'Pt'
        79: 'Au' 'Hg' 'Tl' 'Pb' 'Bi' 'Po' 'At' 'Rn' 'Fr' 'Ra' 'Ac' 'Th' 'Pa'
        92: 'U'  'Np' 'Pu' 'Am' 'Cm' 'Bk' 'Cf' 'Es'
        '''
        parser = argparse.ArgumentParser(description=string,
                formatter_class=RawTextHelpFormatter)
            #formatter_class=ArgumentDefaultsHelpFormatter)
        parser.add_argument('elements', metavar='elements', nargs='+',
            help='list of elements, e.g. Mg; AL cU; si 2 u 7; u o 8 6')
        parser.add_argument('-z', action='store_true', default=False,
            help='get atom number (proton number)')
        parser.add_argument('-n', action='store_true', default=False,
            help='get atom name')
        parser.add_argument('-s', action='store_true', default=False,
            help='get atom symbol')
        parser.add_argument('-m', action='store_true', default=False,
            help='get atom mass (no dimensions)')
        parser.add_argument('-r', action='store_true', default=False,
            help='get atomic reference state')
        parser.add_argument('-t', action='store_true', default=False,
            help='get atomic melting temperature (K)')
        parser.add_argument('-tr', action='store_true', default=False,
            help='get atomic melting temperature rounded up (K)')
        parser.add_argument('-a', action='store_true', default=False,
            help='get all data')
        parser.add_argument('-str', action='store_true', default=False,
            help='output as string and not as python object')
        parser.add_argument('-v','--verbose',
                help='verbose', action='store_true', default=self._verbose)
        args = parser.parse_args()

        if args.verbose:
            print("args:",args)

        self._listin = args.elements
        self._z = args.z
        self._n = args.n
        self._s = args.s
        self._m = args.m
        self._r = args.r
        self._t = args.t
        self._tr = args.tr
        self._a = args.a
        self._str = args.str
        self._verbose = args.verbose


if __name__ == '__main__':
    atom = atom()
    #print "at:",atom
    #print "at:",atom()
    #atom = atom(['u','o',1,2])
    #print "KK",atom.name,atom.mass

    d = np.array([])
    if False:
        for i,s in enumerate(atomic_names):
            print(i,atomic_masses[i]-atomic_masses_wiki[i])
            d = np.append(d,atomic_masses[i]-atomic_masses_wiki[i])
    if False:
        import my
        a=np.array([])
        for ind,element in enumerate(atomic_symbols):
            if element == '':
                x=''
            else:
                #f_element=`grep -i "^[0-9]* $element " $file_we | awk '{print $3}'`
                #x = my.run('grep '+str(element)+' datafile.txt | grep -o "{{sort|[0-9]*|[0-9.]*" | sed "s|.*sort ||" | sed "s|^[0-9]* ||" | sed "s|\|| |g" | xargs -n 1 | tail -1')
                #x = my.run('grep '+str(element)+' datafile.txt | sed "s|\(.*\)\|\|.*|\\1|" | sed "s|\(.*\)\|\|.*|\\1|" | sed "s|\(.*\)\|\|.*|\\1|" | sed "s|\(.*\)\|\|.*|\\1|" | sed "s|\(.*\)\|\|.*|\\1|" | sed "s|\(.*\)\|\|.*|\\1|" | sed "s|.*\|\|||" | sed "s|{{ref.*||" | sed "s|{{sort||" | sed "s|}}||" | sed "s|\|| |g" | xargs -n 1 | tail -1')
                #xx = my.run2('echo '+x+' | sed -s "s|(| |"')
                x = my.run('grep mass /Users/glensk/Thermodynamics/vasp_potentials/PAW-GGA-PBE/'+element+'/POTCAR | awk \'{print $3}\' | sed "s|;||"')
                #x = x.split('(')[0]

                print(element,"x::",x)
                #mp=my.run('getMeltingPoint.sh '+str(element))
            if my.is_float(x):
                x = float(x)
            else:
                x=np.nan
            d = np.append(d,atomic_masses[i]-atomic_masses_wiki[i])
            print(element,", x:",x)
            print("")
            print("")
            a=np.append(a,x)
        np.savetxt("x.dat",a,fmt='%.12f')

