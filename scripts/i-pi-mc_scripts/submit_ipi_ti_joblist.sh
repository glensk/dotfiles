#!/bin/bash

#SBATCH --job-name=NNP-mpi
#SBATCH --get-user-env
#SBATCH --output=_scheduler-stdout.txt
#SBATCH --error=_scheduler-stderr.txt
#SBATCH --nodes=1
#SBATCH --ntasks 36
#SBATCH --time=00-60:29:00

set +e
#export LD_LIBRARY_PATH=
#source $MODULESHOME/init/bash    # necessary for zsh or other init shells
#module load intel intel-mpi intel-mkl fftw python/2.7.16
#export LD_LIBRARY_PATH=/home/glensk/sources/n2p2//lib:${LD_LIBRARY_PATH} # necessary for n2p2
#export OMP_NUM_THREADS=2

#myutils.py -ef ipi_sart_job
#fah.py -ef go_through_all_fah_jobs

hier=`pwd`
echo "hier:$hier"
[ ! -e "joblist_tmp.dat" ] &&  [ ! -e "joblist.dat" ] && echo "one of joblist{,_tmp}.dat has to exist" && exit
[ ! -e "joblist_tmp.dat" ] && cp joblist.dat joblist_tmp.dat
todo=`wc -l $hier/joblist_tmp.dat | awk '{print $1}'`
echo "XX todo1:$todo:"
while [ "$todo" -ge "1" ];do
    echo ""
    cd $hier
    ############################# 
    # check if we are done
    ############################# 
    todo=`wc -l $hier/joblist_tmp.dat | awk '{print $1}'`
    [ "$todo" = "0" ] && echo "XX FINISHED from todo" && break
    

    ############################# 
    # get node status 
    ############################# 
    nextjob=`head -1 $hier/joblist_tmp.dat`
    running_jobs_count=`ps aux | grep glensk | grep i-pi | wc -l`
    cpu2=`top -bn1 | grep "Cpu(s)" | sed "s/.*, *\([0-9.]*\)%* id.*/\1/" | awk '{print (100 - $1)*100}'`
    echo "XX jobs todo $todo, running_jobs_count $running_jobs_count, cpu2:$cpu2, nextjob: $nextjob"

    ############################# 
    # get next job 
    ############################# 
    cd $hier
    if [ "$cpu2" -le "8300" ];then   # 8000 = 80% or 9500 = 95%; 9900 = 99%
        # actually the cpu2 never gets beyond 84% even if too many jobs submitted;
        # the last jobs just "wait" until others are done; seems to be ok strategy
        cd $nextjob
        myutils.py -ef ipi_sart_job &
        sed '1d' $hier/joblist_tmp.dat > tmpfile; mv tmpfile $hier/joblist_tmp.dat # POSIX
        echo "XX started $nextjob"
        sleep 1.0
        cd $hier
    else
        sleep 50 
        echo 'XX sleeping 30 sec'
    fi
done
echo "XX FINISHED THE LOOP"
wait
echo "now all jobs done"
cd $hier
fah.py -ef fah_go_through_all_angK_folder_and_exec_function get_Fah_and_avg_dudl_in_one_ang_K_folder_from_ipi_job
fah.py -ef fah_get_Fah_surface
fah.py -ef fah_get_thermo
