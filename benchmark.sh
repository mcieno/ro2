#!/usr/bin/env sh
# TIME CSV ON STDOUT, NODES CSV ON STDERR5

echo "[*] Building with make all"
make all > /dev/null || exit 1

#timelimit=3600      # 60 minutes
timelimit=360       # 6 minutes
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
    #dummy
    #mtz
    #flow1
    #mtzlazy
    #flow1lazy
    #loopBC
    lazyBC
    lazyBCg
    lazyBCc
    lazyBCcg
)

testbed=(
    data/rat195.tsp
    data/d493.tsp
    data/pr152.tsp
    data/u159.tsp
    data/d657.tsp
    data/lin105.tsp
    data/gil262.tsp
    data/gr96.tsp
    data/ch150.tsp
    data/pr124.tsp
    data/lin318.tsp
    data/u724.tsp
    data/a280.tsp
    data/gr202.tsp
    data/pr226.tsp
    data/kroB200.tsp
    data/linhp318.tsp
    data/gr137.tsp
    data/rd100.tsp
    data/pr144.tsp
    data/pcb442.tsp
    data/pr264.tsp
    data/burma14.tsp
    data/kroA100.tsp
    data/att532.tsp
    data/dummy.tsp
    data/ulysses16.tsp
    data/eil101.tsp
    data/p654.tsp
    data/tsp225.tsp
    data/berlin52.tsp
    data/rat575.tsp
    data/ulysses22.tsp
    data/kroB150.tsp
    data/u574.tsp
    data/ts225.tsp
    data/d198.tsp
    data/eil76.tsp
    data/pr439.tsp
    data/kroC100.tsp
    data/att48.tsp
    data/pr76.tsp
    data/rat99.tsp
    data/kroB100.tsp
    data/eil51.tsp
    data/ch130.tsp
    data/pr299.tsp
    data/st70.tsp
    data/pr107.tsp
    data/gr229.tsp
    data/kroD100.tsp
    data/bier127.tsp
    data/rd400.tsp
    data/kroA200.tsp
    data/gr666.tsp
    data/pr136.tsp
    data/rat783.tsp
    data/gr431.tsp
    data/kroA150.tsp
    data/ali535.tsp
    data/kroE100.tsp
    data/fl417.tsp
)

seeds=(
    2222
    3333
    4444
    #5555
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
            testresult=( $(timeout $timelimit ./bin/tsp $file --model=$model --seed $seed --timelimit $timelimit --nodelimit $nodelimit -j4 --noplot --quiet) )

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
    -P "Nodes, shift 1000"   > /dev/null

echo -e "\n\n[+] All done"
