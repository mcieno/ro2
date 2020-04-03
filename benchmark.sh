#!/usr/bin/env sh

# Benchmark all TSP solvers

models=(
    dummy
    flow1
    mtz
)

testset=(
    data/att48.tsp
    data/eil51.tsp
    data/berlin52.tsp
    data/st70.tsp
    data/eil76.tsp
    data/pr76.tsp
    data/rat99.tsp
    data/rd100.tsp
)

for tspfile in "${testset[@]}"; do
    echo "====================== $tspfile ======================"
    for model in "${models[@]}"; do
        echo "[ $model ]"
        ./bin/tsp $tspfile --model=$model  --noplot
    done
done
