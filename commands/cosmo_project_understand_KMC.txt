$HOME/google_drive/Glensk/KMC-Tutorial/test_20180809

$HOME/google_drive/scripts_epfl/lammps_cosmo/lammps_cosmopc_20180821/lammps/src/lmp_serial_cosmopc_runner-lammps < in.nn 
 ---> resulted in 
ERROR: neighbor list array too small, increase MAXNEIGH and recompile.
--> change in pair_runner.h the #define MAXNEIGH from 192 to 500 in line 26 in the scr folder of lammps


------------------------------------------------------------------------------------------


cd /home/glensk/google_drive/Glensk/KMC-Tutorial/test
../../../scripts_epfl/i-pi-mc/bin/i-pi input-runner.xml&> log &
cat log # -->(need to install numpy with .... probably change to env the i-pi script)
    lmp_serial < in.nn  # this uses the wrong lammps (not the one compiled with runner)
    

    study ~/google_drive/scripts_epfl/i-pi-mc/ipi/engine/motion/dynamics.py
    study ~/google_drive/scripts_epfl/i-pi-mc/ipi/engine/motion/al6xxx_kmc.py ????
    
    STUDY: /Users/glensk/Dropbox/Albert/google_drive/scripts_epfl/i-pi-mc/ipi/engine/motion/al6xxx_kmc.py seems the most relevant to me.



    cd /home/glensk/google_drive/Glensk/KMC-Tutorial/test_3
    i-pi input-runner.xml&> log &
    /home/glensk/google_drive/scripts_epfl/lammps_cosmo/lammps_cosmopc_20180821/lammps/src/lmp_serial_cosmopc_runner-lammps2 < in.nn


------------------------------------------------------------------------------------------

GENERAL: n_{rates} = (n_{types}+1)^{n_{site}} 

## UNARY CASE (214 sites)            = 2^214 ~ 3 10^64
## UNARY CASE only NN (fcc)          = 2^12  ~ 1 10^4
## BINARY CASE only NN (fcc)         = 3^12  ~ 5 10^5
## TERNARY CASE only NN (fcc)        = 4^12  ~ 2 10^7
## QUARTERNARY CASE only NN (fcc)    = 5^12  ~ 2 10^8  (Al-Mg-Si-vac)
## QUARTERNARY CASE only NN+2NN (fcc)= 5^18  ~ 3 10^12 (Al-Mg-Si-vac)

If only 1 chemical species -> 214 possibilities.
If already 2 species -> 
------------------------------------------------------------------------------------------
files in /Users/glensk/Dropbox/Albert/google_drive/Glensk/KMC-Tutorial/ref:
was a run of 200000 steps
init-al202mg6si6.xyz    --> contains the initial starting positions having 214 atoms 
                            (and 2 vac) (maybe this is relaxed? seems different from data.nn)
data.nn (structure.lmp) --> dummy starting configuration for lammps; is never written to simulation.{out,pos_0.xyz} 
simulation.out          --> contains the step and corresponding total energy
                            1column: The current simulation time step.
                            2column: potential energy
KMC_AL6XXX              --> contains for every step one line: the corresponding string of 
                            positions (216 which includes 2*V)
                            and has other ordering then init-al202mg6si6.xyz
                            1column: step but in a different unit
                            2column: ??
                            3column: potential energy
                            4column: string of positions
simulation.pos_0.xyz    --> contains the exact positions for every step (here 200001) 
KMC_QCACHE              --> ??
KMC_ECACHE              --> ?? (here one can apparently quickly check if for a given structure the energy has already been calculated)
in.nn                   --> lammps input file;
------------------------------------------------------------------------------------------
Timing:
-------
run some test calculations and  (~12h single core -> 3223 out of 200003 (==1.6% of long run)
        --> ~on 16 cores, this would take about 2 days.
------------------------------------------------------------------------------------------

Bill's TODO:
- partitioning the domain: so some part would be with NN pot and some with some linear theory.
- not relaxing after every event: I currently dont see how this is possible since i think that this might heavily change the dE's and therefore the rates.
- doing a faster ML for the transition energy barriers vs. local composition:
        - get "real" energy barriers (maybe from NEB) vs. local composition (maybe distinguish 1 and 2 vacancy case)
        - get delta Energy vs. local composition

Done: 
-----
- rerun ref calc
- tried to make appointment with d.g. but no luck, so far no input
- started looking at the details of the code 


Questions general: 
------------------
- are results in ref folder the one's presented in the thesis for the 6x6x6 sc? (is 8x8x8 supercell also there?)
- are there any post processing scripts to analyze the KMC (e.g. precipitate formation).
- which results from the KMC can be compared/validated to experiment?
- concentration of vacancies in thesis is done at 873K, KMC however at 300K. Is this supposed to mimic the behaviour after quenching (as used in a realworld example?)
++- we want to study / do kmc at particular temperature? study several temp? or do some ramping?
- are there any quantities out of the KMC you are particularly interested in? 
- what are first aims: - checking physical sanity of results? checking input assumptions (set diffusion barriers, actual diffusion barriers, temperatures, time evaluations, supercell size)
- should I just start some calcs with this settings and obtain longer KMC runs, larger cells, several seeds?
- should I also take some KMC structures and recalculate with DFT?
- when keeping SF, how long does it take to fit NN pot?


Questions for i-pi KMC:
-----------------------
++- in input: dffusion Barrier for Al defined; diffusion barriers for Mg or Si not defined explicitly; when many neighbors are Mg or Si, explicit barriers for Mg and Si can become important?
- is the given nn potential the first (worse) or the second including KMC structures recalculated by DFT? 
ok- total time... in 

Checks:
-------
- check how results are changed when migration barrier is changed by from 0.52eV to 0.42eV or 0.62eV
- check for few structures the currently approximated barrier vs. actual barrier.
- check how results are change if actual barriers would be calculated
- check alat from 4.06 to 4.045.... ok wobei at 300K this mighte be around 4.06

Checks today: 
- check the Time spent in different positions
- check timing 1 core vs 2 cores
- check positions before and after the geo opt
- check energies for different jumps and list those

------------------------------------------------------------------------------------------
-- Meeting with Michele on 2018_09_26
------------------------------------------------------------------------------------------
- push recalc 
- conc of vac is open quest. 
- 1% si 0.5% mg  rerun + rerun in larger cell. 
- box which can hold one si6mg5 should be enough. (should hold 1FU at experimental concentrations)
TODO:
    - 8x8x8 (512 atoms) cell and put 6si (1.12%) and 5mg (1%) and one vac. (this can be run in DFT)
    - diff random seeds. build average. start 4 seeds.
    - make 1000 atoms cell 10x10x10 
    - average the energies as function of time (get slope of curves) for several seeds. is the 1000 atom cell half as slow? (one could also put 2 vacs in 1000atom cell and ask weather as fast as ...)
    - constellium is interested : do clusters act as vacancy traps. 

    - set up fidis: 
        - how well parallel? try 6 cores or 12 cores. threads. 
        - 4 cores per lammps and 6 lammps nthreads = 6
        - tottime set 3 days
        - KMC 300K? 350? lets ask bill 
        - barriers of diffusion by paper of gulio (should be tried but michele does not think that this influences the KMC significantly)

------------------------------------------------------------------------------------------
gnuplot

  G N U P L O T
  Version 5.2 patchlevel 4    last modified 2018-06-01

  Copyright (C) 1986-1993, 1998, 2004, 2007-2018
  Thomas Williams, Colin Kelley and many others

  gnuplot home:     http://www.gnuplot.info
  faq, bugs, etc:   type "help FAQ"
  immediate help:   type "help"  (plot window: hit 'h')

Terminal type is now 'aqua'
gnuplot> p 'KMC_AL6XXX'
Icon^M                KMC_QCACHE            in.nn                 log                   simulation.out
KMC_AL6XXX            RESTART               init-al202mg6si6.xyz  log.lammps            simulation.pos_0.xyz
KMC_ECACHE            data.nn               input-runner.xml      nohup.out             tmp.xyz
gnuplot> p 'KMC_AL6XXX' u 1:2 w l
gnuplot> p 'KMC_AL6XXX' u 1:3 w l
gnuplot> p 'KMC_AL6XXX' u ($1*2.4e-17):3 w l
gnuplot> p 'KMC_AL6XXX' u ($1*2.4e-17):3 w l, 'simulation.out' u ($1/50000):2 w l
gnuplot> p 'KMC_AL6XXX' u ($1*2.4e-17):3 w l, 'simulation.out' u ($1/200000):2 w l
gnuplot>

------------------------------------------------------------------------------------------
xmgrace
awk '{print $1*2.4e-17,$3}' ../ref/KMC_AL6XXX | x -
------------------------------------------------------------------------------------------

Michele: 
 - Why did previous parallization not worK: definig srun with -c did not work in normal que but did work in debug que. defining srun with -c which used only one core in the end. 
 - I would use all the migration barriers from wolverton paper.
 - 8x8x8 vs 10x10x10 5Mg6Si using for both 5Mg and 6Si  -> keep as was  

 - need --exclusive 

Bill: 
 - what concentrations of Mg and Si exactly?

------------------------------------------------------------------------------------------
Meeting Michele 16. October
------------------------------------------------------------------------------------------
1% Si 1% Mg    8x8x8 1vac
1% Si 1% Mg 10x10x10 1vac
1% Si 1% Mg    8x8x8 2vac
1% Si 1% Mg 10x10x10 2vac
4x random seeds each

------------------------------------------------------------------------------------------
Meeting Michele 19. October
------------------------------------------------------------------------------------------
are vacancies bound to Mg/Si? -> check expecially the big runs
Questions: 
- can i reduce the amout writte to log.i-pi? at a known position, all the rates are known, should go directly to choosing from possible events list. 
- why is not the entire state saved with all the rates and cdf? 
- would it make sense to check (for a few structures) to calculate the real barriers (e.g. by neb) and compare to the estimated barrieres?
- al6xxx_kmc.py line 230 incicates that fist si and then mg. is this mandatory? (probably not since everything is shuffled and we have same si/mg concentration, if concentrations differ, i do not know yet)
- from point of view of comparing to experiment we should rather use 4Mg4Si of 1000 (alloy H) and 6Mg8Si of 1000 (alloy F).
