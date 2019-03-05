see https://linux.die.net/man/1/module

MODULESHOME:
    The location of the master Modules package file directory containing module command initialization scripts, the executable program modulecmd, and a directory containing a collection of master modulefiles.

MODULEPATH:
    The path that the module command searches when looking for modulefiles. Typically, it is set to a default value by the bootstrap procedure. MODULEPATH can be set using 'module use' or by the module initialization script to search group or personal modulefile directories before or after the master modulefile directory.



source $MODULESHOME/init/bash   # to load module command
source $MODULESHOME/init/zsh 
echo $MODULEPATH

