#!/usr/bin/env sh

# Benchmark all TSP solvers

timelimit=1800  # 30 minutes

models=(
    dummy
    mtz
    flow1
    mtzlazy
    flow1lazy
    dummyBB
)

testbed=(
    data/att48.tsp
    data/eil51.tsp
    data/berlin52.tsp
    data/st70.tsp
    data/eil76.tsp
    data/pr76.tsp
    data/rat99.tsp
    data/rd100.tsp
)

bmdir="benchmarks"
bmfile="bm_$(date +%F_%T).csv"

mkdir -p $bmdir || exit 1

echo "Saving benchmark to $bmdir/$bmfile"

echo "${#models[@]} ${models[@]}" | tr -s ' ' ',' | tee "$bmdir/$bmfile"

for tspfile in "${testbed[@]}"; do
    echo -n $tspfile | tee -a "$bmdir/$bmfile"
    for model in "${models[@]}"; do
        testresult=$(timeout $timelimit ./bin/tsp $tspfile --model=$model  --noplot --quiet)
        if [ $? -ne 0 ]; then
            testresult=$timelimit
        fi
        echo -n ",$testresult" | tee -a "$bmdir/$bmfile"
    done
    echo "" | tee -a "$bmdir/$bmfile"
done

python2 perfprof.py              \
    -D ','                       \
    -T 3600                      \
    -S 2                         \
    -M 20                        \
    $bmdir/$bmfile               \
    $bmdir/$bmfile.pdf           \
    -P "all instances, shift 2s"
