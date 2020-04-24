#!/usr/bin/env sh
# See: https://www.dei.unipd.it/bladecluster

echo "[*] Building with make all"
make -f Makefile.blade all > /dev/null || exit 1

timelimit=3600      # 60 minutes
nodelimit=10000000  # 10 million nodes

bmdir="benchmarks"
bmsig="bm_$(date +%F_%T)"

bmfile_nodes="$bmdir/$bmsig.nodes.csv"
bmfile_times="$bmdir/$bmsig.times.csv"

bmfile_nodes_png="$bmdir/$bmsig.nodes.png"
bmfile_times_png="$bmdir/$bmsig.times.png"

mkdir -p $bmdir || exit 1

echo "[*] Saving benchmark to $bmdir/$bmsig.[nodes|times].csv"

models=(
    #Dummy
    #MTZ
    #Flow1
    #LazyMTZ
    #LazyFlow1
    Loop
    Legacy
    Generic
    LegacyConcorde
    GenericConcorde
    LegacyConcordeShallow
    GenericConcordeShallow
    LegacyConcordeRand
    GenericConcordeRand
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
    data/ts225.tsp
    data/tsp225.tsp
    data/pr226.tsp
    data/gr229.tsp
    data/gil262.tsp
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

seeds=(
    21357
    78986
    86874
    10856
    83214
    13476
)

echo "${#models[@]} ${models[@]}" | tr -s ' ' ',' > $bmfile_times
echo "${#models[@]} ${models[@]}" | tr -s ' ' ',' > $bmfile_nodes

echo "[*] Testing ${#models[@]} models (${models[@]}) on ${#testbed[@]} files (${#seeds} seeds each)"


for file in "${testbed[@]}"; do
    for seed in "${seeds[@]}"; do
        echo -e "\n[*] ================ $file : $seed ================"

        echo -n "$file:$seed" >> $bmfile_times
        echo -n "$file:$seed" >> $bmfile_nodes

        for model in "${models[@]}"; do
            testresult=( $(timeout $timelimit ./bin/tsp $file --model=$model --seed $seed --timelimit $timelimit --nodelimit $nodelimit -j16 --noplot --quiet) )

            if [ $? -ne 0 ]; then
                testresult=( $timelimit $nodelimit )
                echo "[-] $model timed out ($timelimit)"
            else
                echo "[>] $model finished in ${testresult[0]}s (${testresult[1]} nodes)"
            fi

            echo -n ",${testresult[0]}" >> $bmfile_times
            echo -n ",${testresult[1]}" >> $bmfile_nodes

            sleep 1
        done

        echo "" >> $bmfile_times
        echo "" >> $bmfile_nodes

    done
done

sleep 1  # let file streams flush

python2 ./perfprof.py        \
    -D ','                   \
    -T $timelimit            \
    -S 1                     \
    -M 5                     \
    $bmfile_times            \
    $bmfile_times_png        \
    -P "Times, shift 1s"     > /dev/null

python2 ./perfprof.py        \
    -D ','                   \
    -T $nodelimit            \
    -S 100                   \
    -M 5                     \
    $bmfile_nodes            \
    $bmfile_nodes_png        \
    -P "Nodes, shift 100"    > /dev/null

echo -e "\n\n[+] All done"
