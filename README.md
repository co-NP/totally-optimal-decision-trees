# Totally optimal decision trees for Boolean functions

Programs and output files for the computations in *On the prevalence of totally
optimal decision trees for Boolean functions* (S. Hussain, manuscript in
preparation).

A decision tree is totally optimal for a set of measures if it minimizes all of
them at once. These programs decide, for the depth `h`, the average depth `havg`
and the number of nodes `L`, whether a Boolean function has such a tree for each
combination of the measures. They also determine whether the function is evasive.

## Files

| File | Purpose |
| --- | --- |
| `tot5.c` | Exact enumeration of all 2^32 Boolean functions of 5 variables |
| `totn.c` | Dynamic program for n <= 12: random samples, all monotone functions (n <= 6), all functions (n <= 4), or one function |
| `tot.py` | Reference implementation using full Pareto fronts of `(h, tpl, L)` |
| `ev.py` | Non-evasive fractions and failure rates for small n (imports `tot.py`) |
| `reproduce.sh` | Builds the programs and regenerates every output file below |
| `results_n5.txt` | Table 4 and the n = 5 row of Table 3 |
| `results_small.txt` | Tables 3 and 5 for small n |
| `results_general.txt` | Table 5 (n = 6), Table 6, and the sampled values in Table 1 |
| `results_n9_seeds.txt` | Two further n = 9 runs with other seeds |
| `results_n9_mt.txt` | An n = 9 run with `std::mt19937_64`, as in the program [BT](https://github.com/co-NP/BT) |
| `runtimes.txt` | Wall-clock time of each run |

## Reproducing the results

The C programs need a C compiler and POSIX threads. The Python scripts need only
Python 3.

```sh
./reproduce.sh
```

This takes about 25 minutes on a 12-core machine. The published files were produced
on an Apple M3 Pro (12 cores, 18 GB) with Apple clang 21. Release `v1.1` is the
version cited in the paper.

## Sampling

Each random truth table consists of independent uniform bits, taken from 64-bit
words of SplitMix64, and functions are drawn independently (with replacement). All
samples use seed 2026, except for the two extra n = 9 runs. Each thread has its own
generator state, so the sampled functions depend on the thread count as well as the
seed: keep `THREADS=12` to obtain the published files exactly. The exhaustive runs
give the same counts for any thread count, although the example functions they print
may differ.

The option `-mt` of `totn` replaces SplitMix64 by `std::mt19937_64`, taking each
truth-table bit from the low bit of one output and using a separate seed for each
function, truncated to 32 bits. It is used only for the run in `results_n9_mt.txt`,
which compares with the samples of the 2016 paper.

## Other commands

```sh
./tot5 -f 0x000c013b               # classify one 5-variable function
./totn -n 6 -f 0x1ce088fdfee2a0af  # one function of n variables
./totn -n 4 -all                   # all functions of 4 variables
./totn -n 5 -mon                   # all monotone functions of 5 variables
./totn -n 8 -s 100000 -seed 1      # random sample
python3 tot.py 5000                # small exhaustive cases and random samples
python3 ev.py
```

## Truth-table conventions

Bit `x` of a truth table is `f(x)`. In `tot.py` and `tot5.c`, `x1` is the most
significant bit of the input, and in `totn.c`, variable `j` is bit `j`. The
classification is invariant under permuting variables, so the two conventions
give the same answer for the same table.

## License

Released under the MIT License. See `LICENSE`.
