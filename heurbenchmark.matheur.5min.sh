#!/usr/bin/env sh

echo "[*] Building with make all"
make -f Makefile all > /dev/null || exit 1

timelimit=60       # 1 minutes for finding the first solution
heurtime=300       # 5 minutes for optimizing the initial solution
totaltimelimit=$(($timelimit+$heurtime+60))

bmdir="benchmarks"
bmsig="$bmdir/bm_$(date +%F_%T)"

bmfile="$bmsig.5min.matheur.csv"
bmfile_png="$bmsig.5min.matheur.png"

mkdir -p $bmdir || exit 1

echo "[*] Saving benchmark to $bmfile"

models=(
    HeurHardfix
    HeurLocalBranching
)

testbed=(
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
            testresult=( $(timeout $totaltimelimit ./bin/tsp $file --model=$model --seed $seed --timelimit $timelimit --heurtime $heurtime -j1 --noplot --quiet) )

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

python3 ./perfprof.py        \
    -D ','                   \
    -S 0                     \
    -M 3                     \
    $bmfile                  \
    $bmfile_png              \
    -X "Cost ratio"          \
    -P "Heuristc solutions ($(($heurtime/60)) min)"  > /dev/null

echo -e "\n\n[+] All done"
