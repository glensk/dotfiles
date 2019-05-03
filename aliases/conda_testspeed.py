#!/usr/bin/env python
import numpy as np
# from https://www.infoworld.com/article/3187484/software/how-does-a-20x-speed-up-in-python-grab-you.html
# to test:
# python $dotfiles/bin/conda_testspeed.py
# source activate intelpy
# python $dotfiles/bin/conda_testspeed.py   # yields about 10x speedup
import time

N = 102400

x = np.linspace(0.0123, 4567.89, N)

def mine(x,Z,func,name):

  print(name)

  start = time.time()

  for z in range ( 0, Z ) :

    y = func(x);

  end = time.time()

  print(N, Z, end - start)

  return

mine(x,10000,np.sin,'np.sin')

mine(x,10000,np.cos,'np.cos')

mine(x,10000,np.tan,'np.tan')

