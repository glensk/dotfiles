verdi computer list
verdi computer show daint
verdi computer setup
    Info: enter "?" for help
    Computer label: fidis
    Hostname: fidis.epfl.ch
    Description []: fidis
    Enable the computer? [True]:
    Transport plugin: ssh
    Scheduler plugin: slurm
    Shebang line (first line of each script, starting with #!) [#!/bin/bash]:
    Work directory on the computer [/scratch/{username}/aiida/]:
    Mpirun command [mpirun -np {tot_num_mpiprocs}]: srun -n {tot_num_mpiprocs}
    Default number of CPUs per machine: 28

verdi computer list
verdi computer configure ssh fidis
    Info: enter "?" for help
    User name [glensk]:
    port Nr [22]:
    Look for keys [False]:
    SSH key file []:
    Error: Path "" does not exist.
    Info: Manually pass a key file if default path is not set in ssh config
    SSH key file []: ?
    Info: Manually pass a key file if default path is not set in ssh config
    SSH key file []: /home/glensk/.ssh/id_rsa
    Connection timeout in s [60]:
    Allow ssh agent [False]: True
    SSH proxy command []:
    Compress file transfers [True]:
    GSS auth [False]:
    GSS kex [False]:
    GSS deleg_creds [False]:
    GSS host [fidis.epfl.ch]:
    Load system host keys [True]:
    Key policy [RejectPolicy]:
    Connection cooldown time (sec) [5]: 30
    Info: Configuring computer fidis for user albert.glensk@gmail.com.
    Success: fidis successfully configured for albert.glensk@gmail.com
    %(aiida) 14:46:43 glensk@lammm@cosmopc18 /home/glensk/aiida

verdi computer test fidis
verdi group show 8
verdi code setup
    Info: enter "?" for help
    Label: pw-v6.3_fidis
    Description []: version 6.3 on fidis
    Default calculation input plugin: quantumespresso.pw
    Installed on target computer? [True]:
    Computer: fidis
    Remote absolute path: /home/glensk/sources/q-e-29c951ebbfc3b2248d0308102d05c4c081db891e/PW/src/pw.x
    
    
    1 #====================================================#
    2 #=               Pre execution script               =#
    3 #====================================================#
    4
    5 #SBATCH --constraint=E5v4
    6
    7 set +e
    8 source $MODULESHOME/init/bash    # necessary in the case of zsh or other init shells
    9 module load intel intel-mpi intel-mkl #quantum-espresso
    0 export OMP_NUM_THREADS=1p

verdi code list
verdi group list
verdi data upf listfamilies
verdi deamon status


 3036  ls /local/scratch/glensk/Dropbox/Albert/scripts/dotfiles/scripts/runner_scripts/
 3037  vi /local/scratch/glensk/Dropbox/Albert/scripts/dotfiles/scripts/runner_scripts/runner_run_mode1_make_fps.py
 3038  exit
 3039  ssh fidis
 3041  cd
 3042  cd aiida
 3044  verdi computer setup
 3045* verdi computer show daint
 3046  verdi computer test
 3048  verdi computer configure fidis
 3049  verdi computer list
 3050  verdi computer configure ssh fidis
 3051  verdi computer test fidis
 3054  verdi code setup
 3057  verdi group show 8
 3058  ~/Dropbox/Albert/git/aiida-alloy/launch_workflow_alalloy_scf.py -h
 3059  ~/Dropbox/Albert/git/aiida-alloy/launch_workflow_alalloy_scf.py --help
 3060  cat ~/Dropbox/Albert/git/aiida-alloy/project_scripts/launch_Al6xxxDB.sh
 3061* /Dropbox/Albert/git/aiida-alloy/project_scripts/
 3062* ~/Dropbox/Albert/git/aiida-alloy/project_scripts/
 3065  dot
 3066  cd ../../
 3067  cd gi
 3068  cd git
 3069  cd aiida-alloy
 3070  verdi code list
 3071  verdi group list
 3072  ls
 3078* vi test_daniel.sh
 3079  verdi data upf list families
 3080  verdi data upf listfamilies
 3081* ./test_daniel.sh
 3083* verdi deamon status
 3084* verdi daemon status
 3088* verdi process list
 3089* verdi calculation list
 3090* verdi calculation list -a -p1
 3091* history

