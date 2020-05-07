#!/usr/bin/env sh

echo "[*] Building with make all"
make -f Makefile.blade all > /dev/null || exit 1

timelimit=600       # 10 minutes for finding the first solution
heurtime=600        # 10 minutes for optimizing it

bmdir="benchmarks"
bmsig="$bmdir/bm_$(date +%F_%T)"

bmfile="$bmsig.heur.csv"
bmfile_png="$bmsig.heur.png"

mkdir -p $bmdir || exit 1

echo "[*] Saving benchmark to $bmfile"

models=(
    HeurHardfix
    HeurLocalBranching
)

testbed=(
    data/burma14.tsp
    data/ulysses16.tsp
    data/ulysses22.tsp
    data/att48.tsp
    data/eil51.tsp
    data/berlin52.tsp
    data/st70.tsp
    data/eil76.tsp
    data/pr76.tsp
    data/gr96.tsp
    data/rat99.tsp
    data/kroA100.tsp
    data/kroB100.tsp
    data/kroC100.tsp
    data/kroD100.tsp
    data/kroE100.tsp
    data/rd100.tsp
    data/eil101.tsp
    data/lin105.tsp
    data/pr107.tsp
    data/pr124.tsp
    data/bier127.tsp
    data/ch130.tsp
    data/pr136.tsp
    data/gr137.tsp
    data/pr144.tsp
    data/ch150.tsp
    data/kroA150.tsp
    data/kroB150.tsp
    data/pr152.tsp
    data/u159.tsp
    data/rat195.tsp
    data/d198.tsp
    data/kroA200.tsp
    data/kroB200.tsp
    data/gr202.tsp
    #data/ts225.tsp
    #data/tsp225.tsp
    #data/pr226.tsp
    #data/gr229.tsp
    #data/gil262.tsp
    data/pr264.tsp
    data/a280.tsp
    data/pr299.tsp
    data/lin318.tsp
    data/linhp318.tsp
    data/rd400.tsp
    data/fl417.tsp
    data/gr431.tsp
    data/pr439.tsp
    data/pcb442.tsp
    data/d493.tsp
    data/att532.tsp
    data/ali535.tsp
    data/u574.tsp
    data/rat575.tsp
    data/p654.tsp
    data/d657.tsp
    data/gr666.tsp
    data/u724.tsp
    data/rat783.tsp
    data/dsj1000.tsp
    ###data/pr1002.tsp
    ###data/u1060.tsp
    ###data/vm1084.tsp
    ###data/pcb1173.tsp
    ###data/d1291.tsp
    ###data/rl1304.tsp
    ###data/rl1323.tsp
    ###data/nrw1379.tsp
    ###data/fl1400.tsp
    ###data/u1432.tsp
    ###data/fl1577.tsp
    ###data/d1655.tsp
    ###data/vm1748.tsp
    ###data/u1817.tsp
    ###data/rl1889.tsp
    ###data/d2103.tsp
    ###data/u2152.tsp
    ###data/u2319.tsp
    ###data/pr2392.tsp
    ###data/pcb3038.tsp
    ###data/fl3795.tsp
    ###data/fnl4461.tsp
    ###data/rl5915.tsp
    ###data/rl5934.tsp
    ###data/pla7397.tsp
    ###data/rl11849.tsp
    ###data/usa13509.tsp
    ###data/brd14051.tsp
    ###data/d15112.tsp
    ###data/d18512.tsp
    ###data/pla33810.tsp
    ###data/pla85900.tsp
)

# Check if seeds are provided as command line args
if [ "$#" -eq 0 ]; then
    seeds=(
        21357
        78986
        86874
    )
else
    seeds=( $@ )
fi

echo "$((${#models[@]}+1)) Optimal,${models[@]}" | tr -s ' ' ',' > $bmfile

echo "[*] Testing ${#models[@]} models (${models[@]}) on ${#testbed[@]} files (${#seeds} seeds each)"


for file in "${testbed[@]}"; do
    testresult=( $(timeout 7200 ./bin/tsp $file --model=GenericConcordeRand -j16 --noplot --quiet) )

    if [ $? -ne 0 ]; then
        # very high cost on timeout
        echo -e "\n[!] Timeout while computing optimal solution for $file. Skipping..."

    else

        trueopt="${testresult[2]}"
        echo -e "\n[O] Optimal solution for $file found: ${seeds[@]}. Starting heuristics..."

        for seed in "${seeds[@]}"; do
            echo -e "\n[*] ================ $file : $seed ================"

            echo -n "$file:$seed,$trueopt" >> $bmfile

            for model in "${models[@]}"; do
                testresult=( $(./bin/tsp $file --model=$model --seed $seed --timelimit $timelimit --heurtime $heurtime -j4 --noplot --quiet) )

                if [ $? -ne 0 ]; then
                    # very high cost on timeout
                    testresult=( $timelimit 0 10000000 )
                    echo "[-] $model timed out ($timelimit)"
                else
                    echo "[>] $model finished in ${testresult[0]}s (cost: ${testresult[2]})"
                fi

                echo -n ",${testresult[2]}" >> $bmfile

                sleep 1
            done

            echo "" >> $bmfile
        done
    fi
done

sleep 1  # let file streams flush

python2 ./perfprof.py        \
    -D ','                   \
    -S 0                     \
    -M 3                     \
    $bmfile                  \
    $bmfile_png              \
    -X "Cost ratio"          \
    -P "Heuristc solutions ($(($heurtime/60)) min)"  > /dev/null

echo -e "\n\n[+] All done"
