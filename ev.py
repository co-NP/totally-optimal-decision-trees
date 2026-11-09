import random, sys
from tot import front, tot_opt
random.seed(2)
def stats(n, funcs):
    tot=ne=fail_hL=fail_hL_ev=fail_ha=fail_aL=0
    for f in funcs:
        fr=front(f,n); r=tot_opt(fr); tot+=1
        nonev = min(p[0] for p in fr) < n
        ne+=nonev
        fail_hL += not r["h,L"]; fail_ha += not r["h,havg"]; fail_aL += not r["havg,L"]
        if not r["h,L"] and not nonev: fail_hL_ev+=1
    print(f"n={n} N={tot} nonevasive={ne} ({ne/tot:.4%}) fail(h,L)={fail_hL} fail(h,havg)={fail_ha} fail(havg,L)={fail_aL} ({fail_aL/tot:.2%}) sanity(evasive&fail h,L)={fail_hL_ev}  RV-bound={2**(-n/2):.4f}")
stats(4, range(1<<16))
for n,N in ((5,20000),(6,5000),(7,600)):
    stats(n,(random.getrandbits(1<<n) for _ in range(N)))
