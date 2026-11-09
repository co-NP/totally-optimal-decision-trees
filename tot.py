"""Exact total-optimality check for exact (alpha=0) decision trees of total Boolean functions.

A function on m variables is a truth-table int of 2^m bits (bit index = input as integer,
x1 = most significant). For each subfunction we compute the 3D Pareto front of
(h, tpl, L) over ALL decision trees (not only strictly optimal ones), so the depth
(non-strictly increasing) is handled exactly.
"""
import sys, random
from functools import lru_cache

def restrict(tt, m, i, b):
    # fix variable i (0 = most significant) to b; returns truth table on m-1 vars
    out = 0; k = 0
    shift = m - 1 - i
    for x in range(1 << m):
        if ((x >> shift) & 1) == b:
            if (tt >> x) & 1:
                out |= 1 << k
            k += 1
    return out

def pareto(pts):
    pts = sorted(set(pts))
    res = []
    for p in pts:
        if not any(q[0] <= p[0] and q[1] <= p[1] and q[2] <= p[2] for q in res):
            res.append(p)
    return tuple(res)

@lru_cache(maxsize=None)
def front(tt, m):
    full = (1 << (1 << m)) - 1
    if tt == 0 or tt == full:
        return ((0, 0, 1),)
    rows = 1 << m
    pts = []
    for i in range(m):
        f0 = front(restrict(tt, m, i, 0), m - 1)
        f1 = front(restrict(tt, m, i, 1), m - 1)
        for a in f0:
            for b in f1:
                pts.append((max(a[0], b[0]) + 1, a[1] + b[1] + rows, a[2] + b[2] + 1))
    return pareto(pts)

PAIRS = {"h,havg": (0, 1), "h,L": (0, 2), "havg,L": (1, 2), "h,havg,L": (0, 1, 2)}

def tot_opt(fr):
    mins = [min(p[k] for p in fr) for k in range(3)]
    return {name: any(all(p[k] == mins[k] for k in idx) for p in fr) for name, idx in PAIRS.items()}

def is_monotone(tt, n):
    for x in range(1 << n):
        if (tt >> x) & 1:
            for j in range(n):
                y = x | (1 << j)
                if not (tt >> y) & 1:
                    return False
    return True

def monotone_all(n):
    if n == 0:
        return [0, 1]
    prev = monotone_all(n - 1)
    half = 1 << (n - 1)
    # f = (f0 when x1=0) low half, (f1 when x1=1) high half, need f0 <= f1
    return [f0 | (f1 << half) for f0 in prev for f1 in prev if f0 & ~f1 == 0]

def run(label, funcs, n):
    cnt = {k: 0 for k in PAIRS}
    ex = {}
    tot = 0
    for f in funcs:
        r = tot_opt(front(f, n))
        tot += 1
        for k, v in r.items():
            if v:
                cnt[k] += 1
            elif k not in ex:
                ex[k] = f
    print(f"{label} n={n} total={tot} " + " ".join(f"{k}:{tot-cnt[k]}fail" for k in PAIRS))
    sys.stdout.flush()
    return ex

if __name__ == "__main__":
    for n in range(0, 5):
        run("ALL", range(1 << (1 << n)), n)
    for n in range(0, 6):
        run("MON", monotone_all(n), n)
    random.seed(1)
    for n in (5, 6):
        N = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
        run("RAND", (random.getrandbits(1 << n) for _ in range(N)), n)
