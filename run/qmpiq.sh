#PBS -l nodes=1:ppn=4
#PBS -q dque
#PBS -N NOB
#PBS -j eo
echo
echo
echo ======================================================
echo Working directory is $PBS_O_WORKDIR
echo
echo Running on host `hostname`
echo
echo Start time is `date`
echo
echo Running Directory is `pwd`
echo
echo This jobs runs on the following nodes :
echo `cat $PBS_NODEFILE`

NOP=$(wc -l $PBS_NODEFILE | awk '{print $1}')

echo
echo This job runs on the following process : $NOP
echo
echo ========================================================

## Loop through nodes before executing my job to be sure there aren't any MPI shm semaphores left
#for i in `cat $PBS_NODEFILE | sort -u` ; do
#      echo "removing IPC shm segments on $i"
#      rsh $i /home/prog/mpich-1.2.7p1/sbin/cleanipcs
#done

rm ./out.1

export P4_GLOBMEMSIZE=1073741824

cd $PBS_O_WORKDIR

cat $PBS_NODEFILE > ./hostfile
#cat $PBS_NODEFILE > ./jobnode

mpirun -machinefile $PBS_NODEFILE -np $NOP ./a.out >> MONITOR.out
