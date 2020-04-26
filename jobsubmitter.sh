#!/bin/bash
# See: https://www.dei.unipd.it/bladecluster
#
# Submit 3 benchmarking jobs, 2 random seeds each
#
# Use qstat to check what jobs are running
# Use qdel JOB_ID [JOB_ID...] to un-submit

qsub -cwd -m ea benchmark.job 21357 78986
sleep 10

qsub -cwd -m ea benchmark.job 86874 10856
sleep 10

qsub -cwd -m ea benchmark.job 83214 13476
echo "[+] All jobs submitted"
