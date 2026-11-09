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
| `results_n5.txt` | Output of `tot5` |
| `results_general.txt` | Output of `totn` for monotone n = 6 and random n = 6 to 12 |

## Building

The C programs need a C compiler and POSIX threads. The Python scripts need only
Python 3.

```sh
cc -O3 -march=native -o tot5 tot5.c -lpthread
cc -O3 -march=native -o totn totn.c -lpthread
```

## Reproducing the results

```sh
./tot5 -t 12 > results_n5.txt

./totn -n 6 -mon -t 12                  >  results_general.txt
./totn -n 6  -s 1000000 -seed 2026 -t 12 >> results_general.txt
./totn -n 7  -s 1000000 -seed 2026 -t 12 >> results_general.txt
./totn -n 8  -s 200000  -seed 2026 -t 12 >> results_general.txt
./totn -n 9  -s 100000  -seed 2026 -t 12 >> results_general.txt
./totn -n 10 -s 30000   -seed 2026 -t 12 >> results_general.txt
./totn -n 11 -s 10000   -seed 2026 -t 12 >> results_general.txt
./totn -n 12 -s 3000    -seed 2026 -t 12 >> results_general.txt
```

Each thread draws from its own random stream, so the sampled functions depend on
the thread count as well as the seed. Use `-t 12` to obtain the published files
exactly. The exhaustive runs give the same counts for any thread count, although
the example functions they print may differ.

The full set takes a few minutes on a 12-core machine, and `tot5` alone takes
about 6 seconds.

Other useful commands:

```sh
./tot5 -f 0x000c013b          # classify one 5-variable function
./totn -n 4 -all              # all functions of 4 variables
./totn -n 5 -mon              # all monotone functions of 5 variables
./totn -n 6 -f 0x44ac6f4fdcf233f2
python3 tot.py 5000           # small exhaustive cases and random samples
python3 ev.py
```

## Truth-table conventions

Bit `x` of a truth table is `f(x)`. In `tot.py` and `tot5.c`, `x1` is the most
significant bit of the input, and in `totn.c`, variable `j` is bit `j`. The
classification is invariant under permuting variables, so the two conventions
give the same answer for the same table.

## License

Released under the MIT License. See `LICENSE`.
