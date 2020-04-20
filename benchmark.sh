#!/usr/bin/env sh
# TIME CSV ON STDOUT, NODES CSV ON STDERR5

echo "[*] Building with make all"
make all > /dev/null || exit 1

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
    #dummy
    #mtz
    #flow1
    #mtzlazy
    #flow1lazy
    loopBC
    lazyBC
    lazyBCg
)

testbed=(
    data/ulysses16.tsp
    data/att48.tsp
    data/berlin52.tsp
    data/st70.tsp
    data/pr76.tsp
    data/rat99.tsp
    data/kroB100.tsp
    data/kroE100.tsp
    data/rd100.tsp
    data/lin105.tsp
    data/bier127.tsp
    data/pr136.tsp
    data/pr144.tsp
    data/kroA150.tsp
    data/pr152.tsp
    data/rat195.tsp
    data/kroA200.tsp
    data/gr202.tsp
    data/tsp225.tsp
    data/gr229.tsp
    #data/pr264.tsp
    #data/pr299.tsp
    #data/fl417.tsp
    #data/pr439.tsp
)

seeds=(
    2222
    3333
    4444
    5555
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
