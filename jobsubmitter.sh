#!/bin/bash

qsub -cwd -m ea benchmark.job 21357 78986
sleep 1
qsub -cwd -m ea benchmark.job 86874 10856
sleep 1
qsub -cwd -m ea benchmark.job 83214 13476

echo "[+] All jobs submitted"
