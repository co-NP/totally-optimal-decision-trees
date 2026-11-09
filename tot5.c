/*
 * Exact enumeration of total optimality for exact (alpha = 0) decision trees
 * of ALL Boolean functions with n = 5 variables (2^32 functions).
 *
 * Same conventions as tot.py: a function on m variables is a truth table of
 * 2^m bits, bit index = input read as an integer with x1 the most significant
 * bit. For each function we need the set of achievable (h, tpl, L) triples;
 * the Pareto front of the children suffices because h combines by max+1 and
 * tpl, L by sums (all increasing), so depth is handled exactly.
 *
 * Fronts of all functions with m <= 4 variables are precomputed. A 5-variable
 * function f = f0 | f1 << 16 (f0: x1 = 0, f1: x1 = 1) is then handled by
 * combining children fronts over its 5 splitting variables.
 *
 * Complementing the output maps decision trees bijectively onto decision
 * trees with the same h, tpl, L, so we only enumerate f with bit 31 = 0 and
 * count every function twice.
 *
 * Build:  cc -O3 -march=native -o tot5 tot5.c -lpthread
 * Usage:  ./tot5            full n = 5 enumeration (all cores)
 *         ./tot5 -t 8       use 8 threads
 *         ./tot5 -f HEX     print the front of one 5-variable function
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>
#include <time.h>

#define MAXF 64 /* max Pareto front size for m <= 4 (checked at runtime) */

typedef struct { uint8_t h, t, l; } pt;    /* depth, total path length, nodes */
typedef struct { uint8_t n; pt p[MAXF]; } front_t;

static front_t *F[5];          /* F[m][tt] for m = 0..4 */
static uint8_t R[4][2][65536]; /* R[v][b][g]: 4-var g with x_{v+1} = b -> 3-var */
static uint8_t MONO4[65536];

static uint32_t restrict_tt(uint32_t tt, int m, int i, int b)
{
    uint32_t out = 0;
    int k = 0, shift = m - 1 - i;
    for (uint32_t x = 0; x < (1u << m); x++)
        if (((x >> shift) & 1) == (uint32_t)b) {
            if ((tt >> x) & 1)
                out |= 1u << k;
            k++;
        }
    return out;
}

static int dominated(const pt *a, const pt *b) /* b dominates a */
{
    return b->h <= a->h && b->t <= a->t && b->l <= a->l;
}

static void add_pareto(front_t *fr, pt q)
{
    for (int i = 0; i < fr->n; i++)
        if (dominated(&q, &fr->p[i]))
            return;
    int j = 0;
    for (int i = 0; i < fr->n; i++)
        if (!dominated(&fr->p[i], &q))
            fr->p[j++] = fr->p[i];
    fr->n = j;
    if (fr->n >= MAXF) {
        fprintf(stderr, "MAXF too small\n");
        exit(1);
    }
    fr->p[fr->n++] = q;
}

static void precompute(void)
{
    for (int m = 0; m <= 4; m++) {
        uint32_t cnt = 1u << (1u << m), full = cnt - 1, rows = 1u << m;
        F[m] = calloc(cnt, sizeof(front_t));
        for (uint32_t tt = 0; tt < cnt; tt++) {
            front_t *fr = &F[m][tt];
            if (tt == 0 || tt == full) {
                fr->n = 1;
                fr->p[0] = (pt){0, 0, 1};
                continue;
            }
            for (int i = 0; i < m; i++) {
                const front_t *a = &F[m - 1][restrict_tt(tt, m, i, 0)];
                const front_t *b = &F[m - 1][restrict_tt(tt, m, i, 1)];
                for (int u = 0; u < a->n; u++)
                    for (int v = 0; v < b->n; v++) {
                        pt q = {(uint8_t)((a->p[u].h > b->p[v].h ? a->p[u].h : b->p[v].h) + 1),
                                (uint8_t)(a->p[u].t + b->p[v].t + rows),
                                (uint8_t)(a->p[u].l + b->p[v].l + 1)};
                        add_pareto(fr, q);
                    }
            }
        }
    }
    for (int v = 0; v < 4; v++)
        for (int b = 0; b < 2; b++)
            for (uint32_t g = 0; g < 65536; g++)
                R[v][b][g] = (uint8_t)restrict_tt(g, 4, v, b);
    for (uint32_t g = 0; g < 65536; g++) {
        int ok = 1;
        for (uint32_t x = 0; x < 16 && ok; x++)
            if ((g >> x) & 1)
                for (int j = 0; j < 4; j++)
                    if (!((g >> (x | (1u << j))) & 1))
                        ok = 0;
        MONO4[g] = (uint8_t)ok;
    }
}

/* failure pattern bits */
enum { NONEV = 1, F_HA = 2, F_HL = 4, F_AL = 8, F_HAL = 16, NPAT = 32 };

/* Evaluate an m-variable function given its children for every split variable. */
static int classify(int m, const front_t **c0, const front_t **c1)
{
    const uint8_t rows = (uint8_t)(1u << m);
    int mh = 255, mt = 255, ml = 255;
    for (int i = 0; i < m; i++) {
        const front_t *a = c0[i], *b = c1[i];
        for (int u = 0; u < a->n; u++)
            for (int v = 0; v < b->n; v++) {
                int h = (a->p[u].h > b->p[v].h ? a->p[u].h : b->p[v].h) + 1;
                int t = a->p[u].t + b->p[v].t + rows;
                int l = a->p[u].l + b->p[v].l + 1;
                if (h < mh) mh = h;
                if (t < mt) mt = t;
                if (l < ml) ml = l;
            }
    }
    int ha = 0, hl = 0, al = 0, hal = 0;
    for (int i = 0; i < m && !hal; i++) {
        const front_t *a = c0[i], *b = c1[i];
        for (int u = 0; u < a->n; u++)
            for (int v = 0; v < b->n; v++) {
                int h = (a->p[u].h > b->p[v].h ? a->p[u].h : b->p[v].h) + 1;
                int t = a->p[u].t + b->p[v].t + rows;
                int l = a->p[u].l + b->p[v].l + 1;
                int eh = h == mh, et = t == mt, el = l == ml;
                ha |= eh && et;
                hl |= eh && el;
                al |= et && el;
                hal |= eh && et && el;
            }
    }
    return (mh < m ? NONEV : 0) | (ha ? 0 : F_HA) | (hl ? 0 : F_HL) |
           (al ? 0 : F_AL) | (hal ? 0 : F_HAL);
}

static int classify5(uint32_t f0, uint32_t f1)
{
    if ((f0 == 0 && f1 == 0) || (f0 == 0xFFFF && f1 == 0xFFFF))
        return 0; /* constants: single leaf, totally optimal, h = 0 < 5 */
    const front_t *c0[5], *c1[5];
    c0[0] = &F[4][f0];
    c1[0] = &F[4][f1];
    for (int v = 0; v < 4; v++) {
        c0[v + 1] = &F[4][R[v][0][f0] | (uint32_t)R[v][0][f1] << 8];
        c1[v + 1] = &F[4][R[v][1][f0] | (uint32_t)R[v][1][f1] << 8];
    }
    return classify(5, c0, c1);
}

static int classify4(uint32_t g)
{
    if (g == 0 || g == 0xFFFF)
        return 0;
    const front_t *c0[4], *c1[4];
    for (int i = 0; i < 4; i++) {
        c0[i] = &F[3][restrict_tt(g, 4, i, 0)];
        c1[i] = &F[3][restrict_tt(g, 4, i, 1)];
    }
    return classify(4, c0, c1);
}

/* constants count as non-evasive only for n >= 1; fix up NONEV for them */
static int const_pat(int n) { return n >= 1 ? NONEV : 0; }

typedef struct {
    int id, nthreads;
    uint64_t hist[NPAT];
    uint32_t example[NPAT];
    int has_ex[NPAT];
} job_t;

static volatile uint32_t done_rows;

static void *worker(void *arg)
{
    job_t *J = arg;
    for (uint32_t f1 = (uint32_t)J->id; f1 < 0x8000; f1 += (uint32_t)J->nthreads) {
        for (uint32_t f0 = 0; f0 < 65536; f0++) {
            int p = classify5(f0, f1);
            if (f1 == 0 && f0 == 0)
                p = const_pat(5);
            J->hist[p]++;
            if (!J->has_ex[p]) {
                J->has_ex[p] = 1;
                J->example[p] = f0 | f1 << 16;
            }
        }
        __sync_fetch_and_add(&done_rows, 1);
        if (J->id == 0 && (f1 / J->nthreads) % 256 == 0)
            fprintf(stderr, "\rprogress %5.1f%%", 100.0 * done_rows / 0x8000);
    }
    return NULL;
}

static void print_table(const char *label, int n, const uint64_t *hist, uint64_t total)
{
    uint64_t c[5] = {0};
    for (int p = 0; p < NPAT; p++)
        for (int k = 0; k < 5; k++)
            if (p & (1 << k))
                c[k] += hist[p];
    printf("%s n=%d total=%llu nonevasive=%llu | fail h,havg=%llu  h,L=%llu  havg,L=%llu  h,havg,L=%llu\n",
           label, n, (unsigned long long)total, (unsigned long long)c[0],
           (unsigned long long)c[1], (unsigned long long)c[2], (unsigned long long)c[3],
           (unsigned long long)c[4]);
}

static void print_pattern_hist(const uint64_t *hist, const uint32_t *ex, const int *has_ex)
{
    printf("pattern (nonev h,havg h,L havg,L h,havg,L; 1 = nonevasive / fails) -> count, example\n");
    for (int p = 0; p < NPAT; p++)
        if (hist[p])
            printf("  %d %d %d %d %d  %12llu  0x%08x\n", !!(p & NONEV), !!(p & F_HA),
                   !!(p & F_HL), !!(p & F_AL), !!(p & F_HAL),
                   (unsigned long long)hist[p], has_ex[p] ? ex[p] : 0);
}

static void self_check(void)
{
    uint64_t hist[NPAT] = {0};
    for (uint32_t g = 0; g < 65536; g++)
        hist[(g == 0 || g == 0xFFFF) ? const_pat(4) : classify4(g)]++;
    print_table("ALL", 4, hist, 65536);

    memset(hist, 0, sizeof hist);
    uint64_t tot = 0;
    for (uint32_t f1 = 0; f1 < 65536; f1++)
        if (MONO4[f1])
            for (uint32_t f0 = 0; f0 < 65536; f0++)
                if (MONO4[f0] && (f0 & ~f1) == 0) {
                    int p = classify5(f0, f1);
                    if ((f0 == 0 && f1 == 0) || (f0 == 0xFFFF && f1 == 0xFFFF))
                        p = const_pat(5);
                    hist[p]++;
                    tot++;
                }
    print_table("MON", 5, hist, tot);
}

int main(int argc, char **argv)
{
    int nthreads = (int)sysconf(_SC_NPROCESSORS_ONLN);
    const char *one = NULL;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-t") && i + 1 < argc)
            nthreads = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-f") && i + 1 < argc)
            one = argv[++i];
    }
    precompute();
    int maxn = 0;
    for (uint32_t g = 0; g < 65536; g++)
        if (F[4][g].n > maxn)
            maxn = F[4][g].n;
    printf("precomputed fronts for m <= 4 (max front size %d)\n", maxn);

    if (one) {
        uint32_t f = (uint32_t)strtoul(one, NULL, 16);
        int p = classify5(f & 0xFFFF, f >> 16);
        printf("f=0x%08x nonev=%d fail: h,havg=%d h,L=%d havg,L=%d h,havg,L=%d\n", f,
               !!(p & NONEV), !!(p & F_HA), !!(p & F_HL), !!(p & F_AL), !!(p & F_HAL));
        return 0;
    }

    self_check();

    time_t t0 = time(NULL);
    pthread_t *th = malloc(sizeof(pthread_t) * nthreads);
    job_t *jobs = calloc(nthreads, sizeof(job_t));
    for (int i = 0; i < nthreads; i++) {
        jobs[i].id = i;
        jobs[i].nthreads = nthreads;
        pthread_create(&th[i], NULL, worker, &jobs[i]);
    }
    uint64_t hist[NPAT] = {0};
    uint32_t ex[NPAT] = {0};
    int has_ex[NPAT] = {0};
    for (int i = 0; i < nthreads; i++) {
        pthread_join(th[i], NULL);
        for (int p = 0; p < NPAT; p++) {
            hist[p] += 2 * jobs[i].hist[p]; /* output complement doubles */
            if (jobs[i].has_ex[p] && !has_ex[p]) {
                has_ex[p] = 1;
                ex[p] = jobs[i].example[p];
            }
        }
    }
    fprintf(stderr, "\rprogress 100.0%%  (%ld s, %d threads)\n", (long)(time(NULL) - t0), nthreads);
    print_table("ALL", 5, hist, 1ull << 32);
    print_pattern_hist(hist, ex, has_ex);
    return 0;
}
