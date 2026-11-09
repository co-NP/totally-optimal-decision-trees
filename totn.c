/*
 * Total optimality of exact decision trees for Boolean functions with n <= 12
 * variables: random sampling, exhaustive monotone enumeration (n <= 6), or a
 * single function.
 *
 * Dynamic programming over the 3^n subcubes (ternary index, digit j of the
 * index belongs to variable j: 0/1 = fixed, 2 = free). For every subcube S and
 * every depth bound d = 0..n we compute
 *   L[d]  = min number of nodes        over trees for S of depth <= d
 *   T[d]  = min total path length      over trees for S of depth <= d
 *   TL[d] = lexicographic min (tpl, L) over trees for S of depth <= d
 *   LT[d] = lexicographic min (L, tpl) over trees for S of depth <= d
 * All four decompose over the root split because the children are independent
 * and h, tpl, L are increasing. With hmin = min{d : L[d] finite}:
 *   (h, L)        totally optimal  iff  L[hmin] == L[n]
 *   (h, havg)     totally optimal  iff  T[hmin] == T[n]
 *   (havg, L)     totally optimal  iff  TL[n].L == L[n]
 *   (h, havg, L)  totally optimal  iff  TL[hmin] == (T[n], L[n])
 *
 * Bit x of the truth table is f(x), variable j is bit j of x.
 *
 * Build:  cc -O3 -march=native -o totn totn.c -lpthread
 * Usage:  ./totn -n 8 -s 100000 [-seed 1] [-t threads]   random functions
 *         ./totn -n 6 -mon                               all monotone functions
 *         ./totn -n 4 -all                               all functions (n <= 4)
 *         ./totn -n 6 -f HEX                             one function
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>
#include <time.h>

#define MAXN 12
#define INF 0xFFFF
#define W(n) (((1u << (n)) + 63) / 64)

typedef uint16_t u16;

static int N;
static uint32_t P3[MAXN + 1];

typedef struct {
    uint8_t *cst; /* 0, 1 constant; 2 mixed */
    u16 *L, *T, *TLt, *TLl, *LTl, *LTt; /* [idx * (N + 1) + d] */
} work_t;

typedef struct {
    int nonev, ok_hl, ok_ha, ok_al, ok_all;
    double price_al_L; /* L of a min-tpl tree with fewest nodes / min L */
    double price_al_T; /* tpl of a min-L tree with smallest tpl / min tpl */
    double price_hl_L; /* min L at min depth / min L */
    double price_ha_T; /* min tpl at min depth / min tpl */
} res_t;

static void work_init(work_t *w)
{
    size_t s = (size_t)P3[N] * (N + 1);
    w->cst = malloc(P3[N]);
    w->L = malloc(s * sizeof(u16));
    w->T = malloc(s * sizeof(u16));
    w->TLt = malloc(s * sizeof(u16));
    w->TLl = malloc(s * sizeof(u16));
    w->LTl = malloc(s * sizeof(u16));
    w->LTt = malloc(s * sizeof(u16));
}

static inline int lexlt(u16 a1, u16 a2, u16 b1, u16 b2)
{
    return a1 < b1 || (a1 == b1 && a2 < b2);
}

static void evaluate(work_t *w, const uint64_t *tt, res_t *r)
{
    const int D = N + 1;
    for (uint32_t idx = 0; idx < P3[N]; idx++) {
        /* decode: point address and free variables */
        uint32_t q = idx, x = 0;
        int nfree = 0, fj = -1, fr[MAXN];
        for (int j = 0; j < N; j++, q /= 3) {
            int dgt = q % 3;
            if (dgt == 2) {
                fr[nfree++] = j;
                fj = j;
            } else if (dgt == 1)
                x |= 1u << j;
        }
        u16 *L = w->L + (size_t)idx * D, *T = w->T + (size_t)idx * D;
        u16 *TLt = w->TLt + (size_t)idx * D, *TLl = w->TLl + (size_t)idx * D;
        u16 *LTl = w->LTl + (size_t)idx * D, *LTt = w->LTt + (size_t)idx * D;
        if (nfree == 0)
            w->cst[idx] = (tt[x >> 6] >> (x & 63)) & 1;
        else {
            uint8_t a = w->cst[idx - 2 * P3[fj]], b = w->cst[idx - P3[fj]];
            w->cst[idx] = (a == b && a != 2) ? a : 2;
        }
        if (w->cst[idx] != 2) {
            for (int d = 0; d < D; d++) {
                L[d] = 1; T[d] = 0;
                TLt[d] = 0; TLl[d] = 1;
                LTl[d] = 1; LTt[d] = 0;
            }
            continue;
        }
        u16 rows = (u16)(1u << nfree);
        L[0] = T[0] = TLt[0] = TLl[0] = LTl[0] = LTt[0] = INF;
        for (int d = 1; d < D; d++) {
            u16 bl = INF, bt = INF, b1t = INF, b1l = INF, b2l = INF, b2t = INF;
            for (int k = 0; k < nfree; k++) {
                int j = fr[k];
                size_t c0 = (size_t)(idx - 2 * P3[j]) * D + d - 1;
                size_t c1 = (size_t)(idx - P3[j]) * D + d - 1;
                if (w->L[c0] == INF || w->L[c1] == INF)
                    continue;
                u16 l = w->L[c0] + w->L[c1] + 1;
                u16 t = w->T[c0] + w->T[c1] + rows;
                u16 tlt = w->TLt[c0] + w->TLt[c1] + rows, tll = w->TLl[c0] + w->TLl[c1] + 1;
                u16 ltl = w->LTl[c0] + w->LTl[c1] + 1, ltt = w->LTt[c0] + w->LTt[c1] + rows;
                if (l < bl) bl = l;
                if (t < bt) bt = t;
                if (lexlt(tlt, tll, b1t, b1l)) { b1t = tlt; b1l = tll; }
                if (lexlt(ltl, ltt, b2l, b2t)) { b2l = ltl; b2t = ltt; }
            }
            L[d] = bl; T[d] = bt; TLt[d] = b1t; TLl[d] = b1l; LTl[d] = b2l; LTt[d] = b2t;
        }
    }
    uint32_t root = P3[N] - 1; /* all digits 2 */
    const u16 *L = w->L + (size_t)root * D, *T = w->T + (size_t)root * D;
    const u16 *TLt = w->TLt + (size_t)root * D, *TLl = w->TLl + (size_t)root * D;
    const u16 *LTt = w->LTt + (size_t)root * D;
    int h = 0;
    while (L[h] == INF)
        h++;
    r->nonev = h < N;
    r->ok_hl = L[h] == L[N];
    r->ok_ha = T[h] == T[N];
    r->ok_al = TLl[N] == L[N];
    r->ok_all = TLt[h] == T[N] && TLl[h] == L[N];
    r->price_al_L = (double)TLl[N] / L[N];
    r->price_al_T = T[N] ? (double)LTt[N] / T[N] : 1.0;
    r->price_hl_L = (double)L[h] / L[N];
    r->price_ha_T = T[N] ? (double)T[h] / T[N] : 1.0;
}

/* ---------- statistics ---------- */

typedef struct {
    uint64_t tot, nonev, f_hl, f_ha, f_al, f_all, f_al_evasive;
    double mx_al_L, mx_al_T, mx_hl_L, mx_ha_T, sum_al_L_fail, sum_al_T_fail;
    uint64_t ex_hl[W(MAXN)], ex_ha[W(MAXN)], ex_al[W(MAXN)];
    int has_hl, has_ha, has_al, has_t, has_x;
    uint64_t triple_only, hl_not_ha; /* see report() */
    uint64_t ex_t[W(MAXN)], ex_x[W(MAXN)];
} stat_t;

static void record(stat_t *s, const res_t *r, const uint64_t *tt)
{
    size_t wb = W(N) * sizeof(uint64_t);
    s->tot++;
    s->nonev += r->nonev;
    if (!r->ok_hl) {
        s->f_hl++;
        if (!s->has_hl++) memcpy(s->ex_hl, tt, wb);
    }
    if (!r->ok_ha) {
        s->f_ha++;
        if (!s->has_ha++) memcpy(s->ex_ha, tt, wb);
    }
    if (!r->ok_al) {
        s->f_al++;
        s->f_al_evasive += !r->nonev;
        s->sum_al_L_fail += r->price_al_L;
        s->sum_al_T_fail += r->price_al_T;
        if (!s->has_al++) memcpy(s->ex_al, tt, wb);
    }
    s->f_all += !r->ok_all;
    if (!r->ok_hl && r->ok_ha) {
        s->hl_not_ha++;
        if (!s->has_x++) memcpy(s->ex_x, tt, wb);
    }
    if (r->ok_hl && r->ok_ha && r->ok_al && !r->ok_all) {
        s->triple_only++;
        if (!s->has_t++) memcpy(s->ex_t, tt, wb);
    }
    if (r->price_al_L > s->mx_al_L) s->mx_al_L = r->price_al_L;
    if (r->price_al_T > s->mx_al_T) s->mx_al_T = r->price_al_T;
    if (r->price_hl_L > s->mx_hl_L) s->mx_hl_L = r->price_hl_L;
    if (r->price_ha_T > s->mx_ha_T) s->mx_ha_T = r->price_ha_T;
}

static void merge(stat_t *a, const stat_t *b)
{
    size_t wb = W(N) * sizeof(uint64_t);
    a->tot += b->tot; a->nonev += b->nonev;
    a->f_hl += b->f_hl; a->f_ha += b->f_ha; a->f_al += b->f_al; a->f_all += b->f_all;
    a->f_al_evasive += b->f_al_evasive;
    a->triple_only += b->triple_only; a->hl_not_ha += b->hl_not_ha;
    if (b->has_t && !a->has_t) { a->has_t = 1; memcpy(a->ex_t, b->ex_t, wb); }
    if (b->has_x && !a->has_x) { a->has_x = 1; memcpy(a->ex_x, b->ex_x, wb); }
    a->sum_al_L_fail += b->sum_al_L_fail; a->sum_al_T_fail += b->sum_al_T_fail;
    if (b->mx_al_L > a->mx_al_L) a->mx_al_L = b->mx_al_L;
    if (b->mx_al_T > a->mx_al_T) a->mx_al_T = b->mx_al_T;
    if (b->mx_hl_L > a->mx_hl_L) a->mx_hl_L = b->mx_hl_L;
    if (b->mx_ha_T > a->mx_ha_T) a->mx_ha_T = b->mx_ha_T;
    if (b->has_hl && !a->has_hl) { a->has_hl = 1; memcpy(a->ex_hl, b->ex_hl, wb); }
    if (b->has_ha && !a->has_ha) { a->has_ha = 1; memcpy(a->ex_ha, b->ex_ha, wb); }
    if (b->has_al && !a->has_al) { a->has_al = 1; memcpy(a->ex_al, b->ex_al, wb); }
}

static void print_tt(const uint64_t *tt)
{
    for (int k = (int)W(N) - 1; k >= 0; k--)
        printf(N >= 6 ? "%016llx" : "%llx", (unsigned long long)tt[k]);
}

static void report(const char *label, const stat_t *s)
{
    double t = (double)s->tot;
    printf("%s n=%d N=%llu\n", label, N, (unsigned long long)s->tot);
    printf("  nonevasive      %12llu  %.6f\n", (unsigned long long)s->nonev, s->nonev / t);
    printf("  fail h,havg     %12llu  %.6f   max T(hmin)/Tmin  %.4f\n", (unsigned long long)s->f_ha, s->f_ha / t, s->mx_ha_T);
    printf("  fail h,L        %12llu  %.6f   max L(hmin)/Lmin  %.4f\n", (unsigned long long)s->f_hl, s->f_hl / t, s->mx_hl_L);
    printf("  fail havg,L     %12llu  %.6f   (evasive among them %llu)\n", (unsigned long long)s->f_al, s->f_al / t, (unsigned long long)s->f_al_evasive);
    printf("  fail h,havg,L   %12llu  %.6f\n", (unsigned long long)s->f_all, s->f_all / t);
    if (s->f_al)
        printf("  havg,L failures: mean L(Tmin)/Lmin %.4f (max %.4f), mean T(Lmin)/Tmin %.4f (max %.4f)\n",
               s->sum_al_L_fail / s->f_al, s->mx_al_L, s->sum_al_T_fail / s->f_al, s->mx_al_T);
    printf("  fail h,havg,L but every pair ok %llu\n", (unsigned long long)s->triple_only);
    printf("  fail h,L but h,havg ok          %llu\n", (unsigned long long)s->hl_not_ha);
    if (s->has_t) { printf("  example fail h,havg,L only: 0x"); print_tt(s->ex_t); printf("\n"); }
    if (s->has_x) { printf("  example fail h,L but not h,havg: 0x"); print_tt(s->ex_x); printf("\n"); }
    if (s->has_ha) { printf("  example fail h,havg: 0x"); print_tt(s->ex_ha); printf("\n"); }
    if (s->has_hl) { printf("  example fail h,L:    0x"); print_tt(s->ex_hl); printf("\n"); }
    if (s->has_al) { printf("  example fail havg,L: 0x"); print_tt(s->ex_al); printf("\n"); }
    fflush(stdout);
}

/* ---------- random sampling ---------- */

static uint64_t splitmix(uint64_t *s)
{
    uint64_t z = (*s += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}

typedef struct {
    int id, nth;
    uint64_t count, seed;
    uint64_t *list; /* monotone mode: list of 32-bit tables of n-1 = 5 variables */
    uint32_t nlist;
    stat_t st;
} job_t;

static void *rand_worker(void *arg)
{
    job_t *J = arg;
    work_t w;
    work_init(&w);
    uint64_t rng = J->seed * 1000003ull + (uint64_t)J->id;
    uint64_t tt[W(MAXN)];
    uint64_t mask = N >= 6 ? ~0ull : ((1ull << (1u << N)) - 1);
    for (uint64_t i = (uint64_t)J->id; i < J->count; i += (uint64_t)J->nth) {
        for (uint32_t k = 0; k < W(N); k++)
            tt[k] = splitmix(&rng) & mask;
        res_t r;
        evaluate(&w, tt, &r);
        record(&J->st, &r, tt);
    }
    return NULL;
}

/* monotone n = 6: f = f0 | f1 << 32 with f0 <= f1 monotone on 5 variables */
static void *mon_worker(void *arg)
{
    job_t *J = arg;
    work_t w;
    work_init(&w);
    uint64_t tt[1];
    for (uint32_t a = (uint32_t)J->id; a < J->nlist; a += (uint32_t)J->nth)
        for (uint32_t b = 0; b < J->nlist; b++) {
            uint64_t f0 = J->list[b], f1 = J->list[a];
            if (f0 & ~f1)
                continue;
            tt[0] = f0 | f1 << 32;
            res_t r;
            evaluate(&w, tt, &r);
            record(&J->st, &r, tt);
        }
    return NULL;
}

static uint32_t monotone_list(int n, uint64_t *out)
{
    if (n == 0) { out[0] = 0; out[1] = 1; return 2; }
    uint64_t *prev = malloc(8000 * sizeof(uint64_t));
    uint32_t m = monotone_list(n - 1, prev), k = 0;
    int half = 1 << (n - 1);
    for (uint32_t a = 0; a < m; a++)
        for (uint32_t b = 0; b < m; b++)
            if (!(prev[b] & ~prev[a]))
                out[k++] = prev[b] | prev[a] << half;
    free(prev);
    return k;
}

int main(int argc, char **argv)
{
    int nth = (int)sysconf(_SC_NPROCESSORS_ONLN), mon = 0, all = 0;
    uint64_t samples = 0, seed = 1;
    const char *one = NULL;
    N = 5;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-n")) N = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-s")) samples = strtoull(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "-seed")) seed = strtoull(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "-t")) nth = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-mon")) mon = 1;
        else if (!strcmp(argv[i], "-all")) all = 1;
        else if (!strcmp(argv[i], "-f")) one = argv[++i];
    }
    if (N < 1 || N > MAXN) { fprintf(stderr, "1 <= n <= %d\n", MAXN); return 1; }
    P3[0] = 1;
    for (int j = 1; j <= MAXN; j++) P3[j] = 3 * P3[j - 1];

    if (one) {
        uint64_t tt[W(MAXN)] = {0};
        const char *h = one + (one[0] == '0' && (one[1] == 'x' || one[1] == 'X') ? 2 : 0);
        int len = (int)strlen(h);
        for (int i = 0; i < len; i++) {
            int c = h[len - 1 - i];
            uint64_t v = c <= '9' ? c - '0' : (c | 32) - 'a' + 10;
            tt[i / 16] |= v << (4 * (i % 16));
        }
        work_t w;
        work_init(&w);
        res_t r;
        evaluate(&w, tt, &r);
        printf("nonev=%d ok: h,havg=%d h,L=%d havg,L=%d h,havg,L=%d\n", r.nonev, r.ok_ha, r.ok_hl, r.ok_al, r.ok_all);
        return 0;
    }

    time_t t0 = time(NULL);
    job_t *jobs = calloc(nth, sizeof(job_t));
    pthread_t *th = malloc(nth * sizeof(pthread_t));
    uint64_t *list = NULL;
    uint32_t nlist = 0;
    const char *label = "RAND";
    void *(*fn)(void *) = rand_worker;
    if (mon) {
        if (N > 6) { fprintf(stderr, "-mon supports n <= 6\n"); return 1; }
        if (N <= 5) { /* single-threaded, small */
            list = malloc(8000 * sizeof(uint64_t));
            nlist = monotone_list(N, list);
            work_t w; work_init(&w);
            stat_t s; memset(&s, 0, sizeof s);
            for (uint32_t i = 0; i < nlist; i++) {
                uint64_t tt[1] = {list[i]};
                res_t r; evaluate(&w, tt, &r); record(&s, &r, tt);
            }
            report("MON", &s);
            return 0;
        }
        list = malloc(8000 * sizeof(uint64_t));
        nlist = monotone_list(5, list);
        label = "MON";
        fn = mon_worker;
    } else if (all) {
        if (N > 4) { fprintf(stderr, "-all supports n <= 4 (use tot5 for n = 5)\n"); return 1; }
        work_t w; work_init(&w);
        stat_t s; memset(&s, 0, sizeof s);
        for (uint64_t f = 0; f < (1ull << (1u << N)); f++) {
            uint64_t tt[1] = {f};
            res_t r; evaluate(&w, tt, &r); record(&s, &r, tt);
        }
        report("ALL", &s);
        return 0;
    }
    for (int i = 0; i < nth; i++) {
        jobs[i].id = i; jobs[i].nth = nth; jobs[i].count = samples; jobs[i].seed = seed;
        jobs[i].list = list; jobs[i].nlist = nlist;
        pthread_create(&th[i], NULL, fn, &jobs[i]);
    }
    stat_t s;
    memset(&s, 0, sizeof s);
    for (int i = 0; i < nth; i++) {
        pthread_join(th[i], NULL);
        merge(&s, &jobs[i].st);
    }
    report(label, &s);
    fprintf(stderr, "(%ld s, %d threads)\n", (long)(time(NULL) - t0), nth);
    return 0;
}
