#!/usr/bin/env sh

echo "[*] Building with make all"
#make -f Makefile all > /dev/null || exit 1

timelimit=1200      # 20 minutes
nodelimit=1000000   # 1 million nodes

bmdir="benchmarks"
bmsig="benchmarks/bm_$(date +%F_%T)"

bmfile_nodes="$bmsig.th2.nodes.csv"
bmfile_times="$bmsig.th2.times.csv"

bmfile_nodes_png="$bmsig.th2.nodes.png"
bmfile_times_png="$bmsig.th2.times.png"

mkdir -p $bmdir || exit 1

echo "[*] Saving benchmark to $bmsig.th2.[nodes|times].csv"

models=(
    Legacy
    Generic
    LegacyConcorde
    GenericConcorde
    LegacyConcordeShallow
    GenericConcordeShallow
    LegacyConcordeRand
    GenericConcordeRand
    GenericConcordeRandWithPatching
)

testbed=(
    data/att48.tsp
    data/eil51.tsp
    data/berlin52.tsp
    data/st70.tsp
    data/eil76.tsp
    data/pr76.tsp
    data/gr96.tsp
    data/rat99.tsp
    data/kroA100.tsp
    data/lin105.tsp
    data/pr107.tsp
    data/pr124.tsp
    data/bier127.tsp
    data/ch130.tsp
    data/u159.tsp
    data/d198.tsp
    data/kroA200.tsp
    data/kroB200.tsp
    data/gr202.tsp
    data/pr264.tsp
    data/a280.tsp
    data/pr299.tsp
    data/lin318.tsp
    data/linhp318.tsp
    data/rd400.tsp
)

# Check if seeds are provided as command line args
if [ "$#" -eq 0 ]; then
    seeds=(
        91824
        71826
        11924
        27389
        40652
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
            testresult=( $(timeout $timelimit ./bin/tsp $file --model=$model --seed $seed --timelimit $timelimit --nodelimit $nodelimit -j2 --noplot --quiet) )

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

python3 ./perfprof.py        \
    -D ','                   \
    -T $timelimit            \
    -S 1                     \
    -M 5                     \
    $bmfile_times            \
    $bmfile_times_png        \
    -P "Times 2th, shift 1s"     > /dev/null

python3 ./perfprof.py        \
    -D ','                   \
    -T $nodelimit            \
    -S 100                   \
    -M 5                     \
    $bmfile_nodes            \
    $bmfile_nodes_png        \
    -P "Nodes 2th, shift 100"    > /dev/null

echo -e "\n\n[+] All done"
