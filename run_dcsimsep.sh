# This job needs 1 compute node with 4 processor per node.
#PBS -l nodes=1:ppn=4
# It should be allowed to run for up to 30 hour.
#PBS -l walltime=30:00:00
# ask for memory
#PBS -l pmem=6gb,pvmem=8gb
# Name of job.
#PBS -N myjob
# Join STDERR TO STDOUT.  (omit this if you want separate STDOUT AND STDERR)
#PBS -j oe   
# Send me mail on job start, job end and if job aborts
#PBS -M phines@uvm.edu
#PBS -m bea
# ask for a MATLAB license
#PBS -W GRES:matlab

cd $HOME/matlab/PowerSystems/dcsimsep
echo "This is myjob running on " `hostname`
matlab -nodesktop -r 'sim_all_k2345, exit'
echo 'done'
