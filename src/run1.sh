# Request nodes/processors
#PBS -l nodes=1:ppn=1
# Request wall time
#PBS -l walltime=30:00:00
# ask for memory
#PBS -l pmem=6gb,pvmem=12gb
# Name of job.
#PBS -N n2s
# Join STDERR TO STDOUT.  (omit this if you want separate STDOUT AND STDERR)
#PBS -j oe
# Send me mail on job start, job end and if job aborts
####    #PBS -M phines@uvm.edu
####    #PBS -m bea

cd ~/dcsimsep
sh ./sim_all_n2s $RUN_ID $N_RUNS > ml_$RUN_ID.out

 
