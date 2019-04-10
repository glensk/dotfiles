# to see which sidepackages python is using
import site; site.getsitepackages()

#to find path to module
print a_module.__file__


import sys
print sys.path

# !!!!! Great Stuff like in mathematika!!!!!  from: http://stackoverflow.com/questions/3676805/python-multidimensional-list-how-to-grab-one-dimension
        >>> import numpy as np
        >>> foo = np.array([[0,1,2],[3,4,5],[6,7,8]])

In [21]: foo
Out[21]:
array([[0, 1, 2],
       [3, 4, 5],
       [6, 7, 8]])

In [22]: foo[0]
Out[22]: array([0, 1, 2])

In [36]: foo[0,2]
Out[36]: 2


In [23]: foo[1]
Out[23]: array([3, 4, 5])

In [34]: foo[[0,2]]
Out[34]:
array([[0, 1, 2],
       [6, 7, 8]])


In [14]: foo[:2]
Out[14]:
array([[0, 1, 2],
       [3, 4, 5]])

In [15]: foo[2:]
Out[15]: array([[6, 7, 8]])

n [16]: foo[:,0]
Out[16]: array([0, 3, 6])

In [18]: foo[:,0:1]
Out[18]:
array([[0],
       [3],
       [6]])

In [19]: foo[:,0:2]
Out[19]:
array([[0, 1],
       [3, 4],
       [6, 7]])

In [20]: foo[:,0:3]
Out[20]:
array([[0, 1, 2],
       [3, 4, 5],
       [6, 7, 8]])


----------------------Mathematica style---------------------------
In [29]: foo[[0,2]]
Out[29]:
array([[0, 1, 2],
       [6, 7, 8]])

In [30]: foo[[0,2],1]
Out[30]: array([1, 7])


shutil.rmtree(nonemptydir)
os.remove(path)   ## i think this is for files
os.chdir(path)      ## cd path
os.getcwd()         ## pwd
os.rename("oldname","newname")  
os.makedirs(dir)    ## mkdir

### cp file1 file2
import shutil
shutil.copyfile(file1,file2)
shutil.copy(file1,tofolder)


    filew = open("KPOINTS", "w")
    filew.write("K-Points\n")
    filew.close()

with open("x.txt") as f:
    data = f.read()
    do something with data
f.closed


with open(NBANDS_file,'w') as f:
     for line in atoms_nbands:
         f.write(str(line)+"\n")

###################################################
###################################################
python qtconsole --pylab=inline

defines everything in a range
def xx(n):
    if 3 < n < 6:
        print "in"
    else:
        print "out"

############# to measure time .... example in pot_energy_forces.py
python -m cProfile -o prof /Users/glensk/Thermodynamics/python_thermodynamics/pot_energy_forces.py -dudlref -s 20 -> writes fiel 'prof'
python -m cProfile -o prof ~/Thermodynamics/python_thermodynamics/lammps_pos_to_sum.py -N 10 -a 4.14 -dt 40 -q 8 8 20 -ps -calclast   -> schreibt file 'prof'
python -m cProfile -s cumtime ~/Thermodynamics/python_thermodynamics/lammps_pos_to_sum.py -N 10 -a 4.14 -dt 40 -q 8 8 20 -ps -calclast  -> writes timings to screen

to evaluate file 'prof'
pstats.Stats('prof').strip_dirs().sort_stats("cumulative").print_stats(30)'
