#!/usr/bin/env sh

# Benchmark all TSP solvers

timelimit=1800     # 30 minutes
nodelimit=10000000 # 10 millions

models=(
    #dummy
    #mtz
    #flow1
    #mtzlazy
    #flow1lazy
    dummyBB
	dummyBBf
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
    #data/rat195.tsp
    #data/kroA200.tsp
    #data/gr202.tsp
    #data/tsp225.tsp
    #data/gr229.tsp
    #data/pr264.tsp
    #data/pr299.tsp
    #data/fl417.tsp
    #data/pr439.tsp
)

seeds=(
    1234
    4321
    1111
    9999
)

bmdir="benchmarks"
bmfile="bm_$(date +%F_%T).csv"

mkdir -p $bmdir || exit 1

echo "Saving benchmark to $bmdir/$bmfile"

echo "${#models[@]} ${models[@]}" | tr -s ' ' ',' | tee "$bmdir/$bmfile"

for tspfile in "${testbed[@]}"; do
    for seed in "${seeds[@]}"; do
        echo -n "$tspfile:$seed" | tee -a "$bmdir/$bmfile"
        for model in "${models[@]}"; do
            testresult=( $(./bin/tsp $tspfile --model=$model --seed $seed --timelimit $timelimit --nodelimit $nodelimit --noplot --quiet) )
            if [ $? -ne 0 ]; then
                testresult=( $timelimit $nodelimit )
            fi
            echo -n ",${testresult[0]}" | tee -a "$bmdir/$bmfile"
        done
        echo "" | tee -a "$bmdir/$bmfile"
    done
done

sleep 1  # let file streams flush

python2 perfprof.py                \
    -D ','                         \
    -T $timelimit                  \
    -S 0.5                         \
    -M 3                           \
    $bmdir/$bmfile                 \
    $bmdir/$bmfile.png             \
    -P "all instances, shift 0.5s"
