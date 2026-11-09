#!/usr/bin/env bash
# Rebuilds the programs and regenerates every output file used in the paper.
# The sampled functions depend on the thread count, so keep THREADS=12 to
# obtain the published files exactly. Takes about 25 minutes on 12 cores.
set -euo pipefail
cd "$(dirname "$0")"

THREADS=12
SEED=2026

cc -O3 -march=native -o tot5 tot5.c -lpthread
cc -O3 -march=native -o totn totn.c -lpthread

: > runtimes.txt
run() { # run OUTPUT_FILE COMMAND...
    local out=$1 start=$SECONDS
    shift
    "$@" >> "$out" 2> /dev/null
    echo "$((SECONDS - start)) s: $*" >> runtimes.txt
}

: > results_n5.txt
run results_n5.txt ./tot5 -t "$THREADS"

: > results_small.txt
for n in 1 2 3 4; do run results_small.txt ./totn -n "$n" -all; done
for n in 1 2 3 4 5; do run results_small.txt ./totn -n "$n" -mon; done

: > results_general.txt
run results_general.txt ./totn -n 6 -mon -t "$THREADS"
while read -r n samples; do
    run results_general.txt ./totn -n "$n" -s "$samples" -seed "$SEED" -t "$THREADS"
done <<SIZES
6 1000000
7 1000000
8 5000000
9 100000
10 30000
11 100000
12 3000
SIZES

: > results_n9_seeds.txt
for seed in 7 31337; do
    run results_n9_seeds.txt ./totn -n 9 -s 100000 -seed "$seed" -t "$THREADS"
done
