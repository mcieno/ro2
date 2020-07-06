#!/usr/bin/env sh

echo "[*] Building with make all"
make -f Makefile all > /dev/null || exit 1

timelimit=2400      # 40 minutes
nodelimit=1000000   # 1 million nodes

bmdir="benchmarks"
bmsig="benchmarks/bm_$(date +%F_%T)"

bmfile_nodes="$bmsig.compact.nodes.csv"
bmfile_times="$bmsig.compact.times.csv"

bmfile_nodes_png="$bmsig.compact.nodes.png"
bmfile_times_png="$bmsig.compact.times.png"

mkdir -p $bmdir || exit 1

echo "[*] Saving benchmark to $bmsig.compact.[nodes|times].csv"

models=(
    MTZ
    Flow1
    LazyMTZ
    LazyFlow1
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
)

# Check if seeds are provided as command line args
if [ "$#" -eq 0 ]; then
    seeds=(
        91820
        74125
        80478
    )
else
    seeds=( $@ )
fi


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
    -P "Nodes, shift 100"    > /dev/null

echo -e "\n\n[+] All done"
