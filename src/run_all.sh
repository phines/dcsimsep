#!/bin/bash
qsub_script="run1.sh"

nruns=50

for run_id in $(eval echo {1..$nruns})
do
	cmd="qsub -v LD_LIBRARY_PATH,XAPPLRESDIR,RUN_ID=$run_id,N_RUNS=$nruns $qsub_script"
	echo "$cmd"
	$cmd
done
