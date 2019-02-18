10000000000 is a high prior

qstat -j <jobID>		#information of job
qstat		look at que status
qstat -u "*"    zeigt alle user
qstat -j "*"    zeigt alle in der que befindlichen jobs von allen usern
qstat -u '*' | grep " r " | grep cmdft              shows all running cmdft jobs


qhost  oder qhost -q -h zeitgt alle infos der im netzwerk befindlichen rechner
qhost -q -h cmmd100   zeigt die DERZEITIGE belegung an (auch seriell/parallel jobs)

qdel jobid	delete job

qalter -js 1 63651  (where 63651 is the job id) makes this job being calculated first
qalter 80584 -q parallel  ## definiert im nachhinein irgendeine que (hier die que namens parallel)
qalter 88888 -pe mpie24 24 -l cmmd   >> chagees job from cmdft to cmmd
qalter 88888 -pe mpie8 8 -l cmdft    >> chagees job from cmdft to cmmd


qls -r | grep glensk  --> wie lange laeuft dre job schono?
qls -u 		show all users qued
qls -q		show all ques and usage
qls -w glensk	status meiner auftraege
qls -n  list all nodes and status


specify system (cmmd cmdft)

qsub -l cmdft pe mpi 24 skript
qsub -l cmmd

qsub -pe mpi $cores /home/glensk/scripts/vasp_/runvasp/vasp.cmmd.par  ## submitted nur auf cmmd
qsub -pe mpi* $cores /home/glensk/scripts/vasp_/runvasp/vasp.cmmd.par ## submitted auf cmmd und cmdft


qconf -spl                  show all parallel environments
qconf -sp mpi               einstellungen fuer parallel environmet mpi

in der sge (sungridengine) gibt es 
  -hardresources
  -softresources (nach prioritaet)


qsub -pe "serial" 16 -l cmdft --> FUNKTIONIERT NICHT!!! serial kann nur fure einen knoten angesetzt werden
qsub -pe "mpi*" 24 -l cmdft sollte aber gehen laut alexej


Dear cluster users,

please note, that (as explained in our wiki) `mpi` parallel environment is deprecated. Therefore, all waiting jobs requesting `mpi` will never be spooled. As a result, despite there are many waiting jobs, the cluster is not fully loaded.

One should use either `mpi24` (mpi8) for cmmd (cmdft), or `mpi*` instead. You can change the requested parallel environment for already submitted jobs with `qalter -pe mpi* nSlots job_id`.
qstat -pri

if job is in Eqw (eqw) error status:
    qmod -cj <jobID>          to clear the error status
