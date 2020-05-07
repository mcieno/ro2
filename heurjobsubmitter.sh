#!/bin/bash
# See: https://www.dei.unipd.it/bladecluster
#
# Submit 3 heur-benchmarking jobs, 2 random seeds each
#
# Use qstat to check what jobs are running
# Use qdel JOB_ID [JOB_ID...] to un-submit

qsub -cwd -m ea heurbenchmark.job 1728 4612
sleep 10

qsub -cwd -m ea heurbenchmark.job 7849 4988
sleep 10

qsub -cwd -m ea heurbenchmark.job 7302 3890

echo "[+] All jobs submitted"
