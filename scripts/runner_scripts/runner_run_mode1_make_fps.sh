#!/bin/sh

###################################################
echo "make fps to get structures to calculate"
echo "cursel.distances.dat shows the distances"
###################################################
# using runner_scratch_new_data_all function_data (3.2GB) and logfile 
#CurSel.py -t 1e-3 --landmarks 20 ../function.data ../logfile_mode1.1    -> 569 sec on 20 struct 
#CurSel.py -t 1e-3 --landmarks 40 ../function.data ../logfile_mode1.1    -> 616 sec on 40 struct 
#CurSel.py -t 1e-3 --landmarks 80 ../function.data ../logfile_mode1.1    -> 553 sec on 80 struct 
#CurSel.py -t 1e-3 --landmarks 800 ../function.data ../logfile_mode1.1   -> 578 sec on 800 struct 
#CurSel.py -t 1e-3 --landmarks 2000 ../function.data ../logfile_mode1.1  -> 625 sec on 2000 struct 
#CurSel.py -t 1e-3 --landmarks `grep -c begin ../input.data` ../function.data ../logfile_mode1.1
#                                                                        -> 577 sec on 2500 struct


#frame_selector.py input.data.all precomp cursel.landmarks  # crates input.dat_selected.data with 40 strut

# This makes apparently the full list out of which a subset can be chosen.
CurSel.py -t 1e-3 --landmarks `grep -c begin input.data` function.data logfile_mode1.1
head -20 cursel.landmarks > cursel.landmarks.20
head -40 cursel.landmarks > cursel.landmarks.40

frame_selector.py input.data precomp cursel.landmarks.20
mv input_selected.data input_selected.data.20

frame_selector.py input.data precomp cursel.landmarks.40
mv input_selected.data input_selected.data.40

# this will create a file cursel.landmarks which contains indexes;
# x cursel.distances.dat (out.dat) will then give the distances
