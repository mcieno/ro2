#!/usr/bin/env sh

echo "[*] Building with make all"
make -f Makefile all > /dev/null || exit 1

timelimit=90       # 1.5 minutes for finding the first solution
heurtime=60        # 1 minutes for optimizing the initial solution

bmdir="benchmarks"
bmsig="$bmdir/bm_$(date +%F_%T)"

bmfile="$bmsig.1min.heur.csv"
bmfile_png="$bmsig.1min.heur.png"

mkdir -p $bmdir || exit 1

echo "[*] Saving benchmark to $bmfile"

models=(
    HeurNearestNeighbor
    HeurGRASP
    HeurInsertion
    HeurConvHullInsertion
    HeurGRASPWith2OPTRefinement
    HeurTabuSearch
    HeurVNS
    HeurGenetic
    HeurSimulatedAnnealing
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
    data/pr264.tsp
    data/a280.tsp
    data/pr299.tsp
    data/lin318.tsp
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
)

testbedopt=(
    33523.708507
    428.871756
    7544.365902
    677.109609
    544.369053
    108159.438274
    510.886315
    1219.243769
    21285.443182
    22139.074615
    20750.762504
    21294.290821
    22068.758669
    7910.396210
    640.211591
    14382.995933
    44301.683677
    59030.735703
    118293.523816
    6110.722200
    96770.924122
    706.289816
    58535.221761
    6530.902722
    26524.863036
    26127.357889
    73683.640628
    42075.670040
    2333.873188
    15808.652051
    29369.407047
    29440.412221
    486.349386
    49135.004963
    2586.769648
    48194.920103
    42042.535089
    15275.984985
    11914.306105
    1924.155132
    107215.301715
    50783.547514
    35019.220247
    86742.422162
    2009.446807
    36934.771414
    6795.967520
    34646.834710
    48917.596353
    3088.465842
    41907.722275
    8842.994960
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
    # Pop true optimal solution
    trueopt="${testbedopt[0]}"
    testbedopt=("${testbedopt[@]:1}")
    echo -e "\n[O] Optimal solution for $file is $trueopt. Starting heuristics..."

    for seed in "${seeds[@]}"; do
        echo -e "\n[*] ================ $file : $seed ================"

        echo -n "$file:$seed,$trueopt" >> $bmfile

        for model in "${models[@]}"; do
            testresult=( $(timeout $timelimit ./bin/tsp $file --model=$model --seed $seed --timelimit $timelimit --heurtime $heurtime -j4 --noplot --quiet) )

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
