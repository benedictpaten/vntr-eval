/* vntreval.c -- fast exact dynamic programs behind tools/evaluate.py (loaded with ctypes).
 *
 * Character model everywhere: two bases match iff they are the same symbol after upper-casing;
 * A, C, G, T are symbols, and every other byte (N, IUPAC codes) is one extra symbol "N".  So N
 * matches N (as two bases stacked on one graph node would) and nothing else.
 *
 *   ed_unit(a, n, b, m)                 exact global unit-cost Levenshtein distance
 *                                       (Myers 1999 bit-vector, 64-bit blocks, O(n*ceil(m/64))).
 *   ed_affine(a, n, b, m, X, GO, GE)    exact global gap-affine cost (Gotoh): mismatch X, a gap of
 *                                       length k costs GO + GE*k.  Banded around the corridor of
 *                                       diagonals between 0 and m-n, doubling the band until the
 *                                       banded optimum is provably global (see the comment there).
 *   dag_align(...)                      min unit-cost edit distance of a query to any source->sink
 *                                       path of an acyclic graph, with the handle path (traceback).
 *
 * Build: see Makefile (cc -O3 -shared -fPIC -o ../bin/libvntreval.so vntreval.c).
 */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static inline int code(unsigned char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;
    }
}

/* ------------------------------------------------------------------ unit cost (Myers) */

long ed_unit(const char *p, long m, const char *t, long n) {
    if (m == 0) return n;
    if (n == 0) return m;
    long W = (m + 63) / 64;
    uint64_t *peq = calloc(5 * W, sizeof(uint64_t));
    uint64_t *Pv = malloc(W * sizeof(uint64_t)), *Mv = calloc(W, sizeof(uint64_t));
    if (!peq || !Pv || !Mv) { free(peq); free(Pv); free(Mv); return -1; }
    for (long i = 0; i < m; i++) peq[code((unsigned char)p[i]) * W + i / 64] |= (uint64_t)1 << (i % 64);
    for (long w = 0; w < W; w++) Pv[w] = ~(uint64_t)0;
    int lastbit = (int)((m - 1) % 64);
    long score = m;
    for (long j = 0; j < n; j++) {
        const uint64_t *eqr = peq + code((unsigned char)t[j]) * W;
        int hin = 1;                               /* global: D[0][j] = j, so the top delta is +1 */
        for (long w = 0; w < W; w++) {
            uint64_t Eq = eqr[w], pv = Pv[w], mv = Mv[w];
            uint64_t Xv = Eq | mv;
            if (hin < 0) Eq |= 1;
            uint64_t Xh = (((Eq & pv) + pv) ^ pv) | Eq;
            uint64_t Ph = mv | ~(Xh | pv);
            uint64_t Mh = pv & Xh;
            int hb = (w == W - 1) ? lastbit : 63;
            int hout = 0;
            if ((Ph >> hb) & 1) hout = 1;
            else if ((Mh >> hb) & 1) hout = -1;
            Ph <<= 1; Mh <<= 1;
            if (hin < 0) Mh |= 1; else if (hin > 0) Ph |= 1;
            Pv[w] = Mh | ~(Xv | Ph);
            Mv[w] = Ph & Xv;
            hin = hout;
        }
        score += hin;
    }
    free(peq); free(Pv); free(Mv);
    return score;
}

/* ------------------------------------------------------------------ gap-affine (Gotoh) */

#define AINF (INT32_MAX / 4)

/* Banded Gotoh over diagonals d = j - i in [dlo, dhi]; returns the banded optimum. */
static long gotoh_band(const char *a, long n, const char *b, long m, int X, int GO, int GE, long dlo, long dhi) {
    int32_t *H = malloc((m + 1) * sizeof(int32_t)), *F = malloc((m + 1) * sizeof(int32_t));
    if (!H || !F) { free(H); free(F); return -1; }
    for (long j = 0; j <= m; j++) { H[j] = AINF; F[j] = AINF; }
    /* row 0 */
    long jhi0 = dhi < m ? dhi : m;
    H[0] = 0;
    for (long j = 1; j <= jhi0; j++) H[j] = GO + GE * (int32_t)j;
    long pjlo = 0, pjhi = jhi0;                   /* previous row's valid range */
    for (long i = 1; i <= n; i++) {
        long jlo = i + dlo; if (jlo < 0) jlo = 0;
        long jhi = i + dhi; if (jhi > m) jhi = m;
        if (jlo > jhi) { free(H); free(F); return AINF; }
        int ca = code((unsigned char)a[i - 1]);
        /* H/F currently hold row i-1 over [pjlo, pjhi]; cells outside are AINF */
        int32_t diag = (jlo >= 1 && jlo - 1 >= pjlo && jlo - 1 <= pjhi) ? H[jlo - 1] : AINF;
        int32_t E = AINF, hleft = AINF;
        /* invalidate row i-1 cells left of this row's range (they cannot be read again) */
        for (long j = pjlo; j < jlo && j <= pjhi; j++) { H[j] = AINF; F[j] = AINF; }
        for (long j = jlo; j <= jhi; j++) {
            int32_t hup = (j >= pjlo && j <= pjhi) ? H[j] : AINF;
            int32_t fup = (j >= pjlo && j <= pjhi) ? F[j] : AINF;
            int32_t f = fup + GE, t = hup + GO + GE;
            if (t < f) f = t;
            int32_t h;
            if (j == 0) {
                h = GO + GE * (int32_t)i;          /* first column: a[0..i) deleted */
                f = h;                            /* the vertical gap state equals it */
                E = AINF;
            } else {
                int32_t e = E + GE; t = hleft + GO + GE;
                if (t < e) e = t;
                E = e;
                int32_t s = diag + ((ca == code((unsigned char)b[j - 1])) ? 0 : X);
                h = s; if (e < h) h = e; if (f < h) h = f;
            }
            if (h > AINF) h = AINF;
            if (f > AINF) f = AINF;
            if (E > AINF) E = AINF;
            diag = hup;                           /* H[i-1][j] becomes the diagonal of j+1 */
            H[j] = h; F[j] = f; hleft = h;
        }
        for (long j = jhi + 1; j <= pjhi; j++) { H[j] = AINF; F[j] = AINF; }
        pjlo = jlo; pjhi = jhi;
    }
    long r = (m >= pjlo && m <= pjhi) ? H[m] : AINF;
    free(H); free(F);
    return r;
}

/* Exactness of the band: an alignment that leaves the corridor [min(0,D)-w, max(0,D)+w], D = m-n,
 * reaches a diagonal at distance w+1 beyond it; the diagonal starts at 0 and ends at D and moves by
 * one per gap base, so such an alignment has >= 2(w+1) + |D| gap bases and costs
 * >= GO + GE*(2(w+1) + |D|).  If the banded optimum is <= that bound, it is the global optimum. */
long ed_affine(const char *a, long n, const char *b, long m, int X, int GO, int GE) {
    if (n == 0 && m == 0) return 0;
    if (n == 0) return GO + (long)GE * m;
    if (m == 0) return GO + (long)GE * n;
    long D = m - n, aD = D < 0 ? -D : D;
    long w = 64;
    long maxw = (n > m ? n : m) + 1;
    for (;;) {
        long dlo = (D < 0 ? D : 0) - w, dhi = (D > 0 ? D : 0) + w;
        long c = gotoh_band(a, n, b, m, X, GO, GE, dlo, dhi);
        if (c < 0) return -1;
        long lb = GO + (long)GE * (2 * (w + 1) + aD);
        if (c <= lb || w >= maxw) return c;
        w *= 2;
        if (w > maxw) w = maxw;
    }
}

/* ------------------------------------------------------------------ query vs DAG */

/* Handles 0..H-1 are in topological order.  Handle h spells seq[off[h] .. off[h+1]).
 * Its predecessors are padj[pstart[h] .. pstart[h+1]), every one < h.
 * src[h]: a path may start at h (entry column j -> j: query prefix inserted before the path).
 * snk[h]: a path may end at h.
 * prio[h]: tie-break of the traceback. Among equally good sinks and predecessors the one with the
 * highest priority is taken (evaluate.py: CHM13 handles first, then handles on more panel paths), so
 * the path returned does not depend on the order of the GFA's lines.
 * ub: an upper bound on the answer (evaluate.py passes the distance to the closest panel path, which
 * is itself a source->sink path), or -1 for none.  Cells whose value exceeds ub cannot lie on an
 * optimal alignment (values never decrease along an alignment), so they are pruned: each handle keeps
 * only the range of query columns [lo, hi] where its exit column is <= ub.  This is exact, and it
 * shrinks memory and time from H x m to H x (the band the near-optimal paths use).
 * Returns the distance, and the handle path in path_out[0..*path_len) (source first), or
 *   -2 when the stored columns would exceed max_bytes, -3 when the optimum is >= 65535 (uint16
 *   columns; only without a smaller ub), -4 on a traceback inconsistency (a bug), -5 when no source
 *   reaches a sink, -6 when ub is not an upper bound (no sink within it).
 * Traceback re-runs the DP inside each node on the path over columns [j - L - best - 1, j] carrying,
 * per cell, the entry column its optimal path came from, so it needs O(L + best) memory per node. */
#define SAT 65535
#define DINF (1 << 30)

typedef struct { long lo, hi; uint16_t *v; } Col;      /* v[j - lo] for j in [lo, hi]; v == NULL: empty */

static inline int32_t colv(const Col *c, long j) {
    return (c->v && j >= c->lo && j <= c->hi && c->v[j - c->lo] < SAT) ? c->v[j - c->lo] : DINF;
}

long dag_align(int H, const char *seq, const long *off, const int *pstart, const int *padj,
               const unsigned char *src, const unsigned char *snk, const int *prio, const char *q, long m,
               long ub, double max_bytes, int *path_out, int *path_len) {
    *path_len = 0;
    long cap = (ub < 0 || ub >= SAT - 1) ? SAT - 1 : ub;    /* prune cells > cap */
    Col *out = calloc(H, sizeof(Col));
    int32_t *c = malloc((m + 2) * sizeof(int32_t)), *tmp = malloc((m + 2) * sizeof(int32_t));
    unsigned char *qc = malloc(m + 1);
    if (!out || !c || !tmp || !qc) { free(out); free(c); free(tmp); free(qc); return -2; }
    for (long j = 0; j < m; j++) qc[j] = (unsigned char)code((unsigned char)q[j]);
    double used = 0;
    long best = -5;
    int cur = -1;
    int err = 0;
    for (int h = 0; h < H && !err; h++) {
        /* entry column: union of the predecessors' ranges (and [0, cap] for a source) */
        long clo = m + 1, chi = -1;
        if (src[h]) { clo = 0; chi = cap < m ? cap : m; }
        for (int k = pstart[h]; k < pstart[h + 1]; k++) {
            const Col *g = &out[padj[k]];
            if (!g->v) continue;
            if (g->lo < clo) clo = g->lo;
            if (g->hi > chi) chi = g->hi;
        }
        if (chi < clo) continue;                        /* unreachable within the bound */
        for (long j = clo; j <= chi; j++) c[j] = DINF;
        if (src[h]) for (long j = 0; j <= (cap < m ? cap : m); j++) c[j] = (int32_t)j;
        for (int k = pstart[h]; k < pstart[h + 1]; k++) {
            const Col *g = &out[padj[k]];
            if (!g->v) continue;
            for (long j = g->lo; j <= g->hi; j++) { int32_t v = g->v[j - g->lo]; if (v < c[j]) c[j] = v; }
        }
        long L = off[h + 1] - off[h];
        const char *s = seq + off[h];
        int empty = 0;
        for (long k = 0; k < L; k++) {
            int ch = code((unsigned char)s[k]);
            long first = -1, last = -1;
            int32_t left = DINF;
            for (long j = clo; j <= m; j++) {
                int32_t v = DINF, d;
                if (j <= chi && c[j] < DINF) v = c[j] + 1;                       /* node base deleted */
                if (j - 1 >= clo && j - 1 <= chi && c[j - 1] < DINF) {             /* (mis)match */
                    d = c[j - 1] + (qc[j - 1] != ch); if (d < v) v = d;
                }
                if (left < DINF && left + 1 < v) v = left + 1;                     /* query base inserted */
                if (v > cap) v = DINF;
                tmp[j] = v; left = v;
                if (v < DINF) { if (first < 0) first = j; last = j; }
                if (j > chi && v == DINF) break;          /* nothing further right can be <= cap */
            }
            if (first < 0) { empty = 1; break; }
            int32_t *sw = c; c = tmp; tmp = sw;
            clo = first; chi = last;
        }
        if (empty) continue;
        long w = chi - clo + 1;
        used += 2.0 * w;
        if (used > max_bytes) { err = -2; break; }
        out[h].v = malloc(w * sizeof(uint16_t));
        if (!out[h].v) { err = -2; break; }
        out[h].lo = clo; out[h].hi = chi;
        for (long j = clo; j <= chi; j++) out[h].v[j - clo] = c[j] >= SAT ? SAT : (uint16_t)c[j];
        if (snk[h]) {
            int32_t v = colv(&out[h], m);
            if (v < DINF && (best < 0 || v < best || (v == best && prio[h] > prio[cur]))) { best = v; cur = h; }
        }
    }
    if (err) best = err;
    else if (best < 0) {
        int any_sink = 0;
        for (int h = 0; h < H; h++) if (snk[h]) any_sink = 1;
        best = !any_sink ? -5 : (cap < SAT - 1 ? -6 : -3);
    }
    if (best >= 0) {
        /* traceback */
        long j = m;
        int np = 0;
        long W0 = 0;
        int32_t *in = NULL, *val = NULL, *org = NULL, *nval = NULL, *norg = NULL;
        for (;;) {
            path_out[np++] = cur;
            long L = off[cur + 1] - off[cur];
            const char *s = seq + off[cur];
            long lo = j - L - best - 1; if (lo < 0) lo = 0;
            long W = j - lo + 1;
            if (W > W0) {
                free(in); free(val); free(org); free(nval); free(norg);
                in = malloc(W * sizeof(int32_t)); val = malloc(W * sizeof(int32_t)); org = malloc(W * sizeof(int32_t));
                nval = malloc(W * sizeof(int32_t)); norg = malloc(W * sizeof(int32_t));
                W0 = W;
            }
            for (long x = 0; x < W; x++) in[x] = DINF;
            if (src[cur]) for (long x = 0; x < W; x++) if (lo + x <= cap) in[x] = (int32_t)(lo + x);
            for (int k = pstart[cur]; k < pstart[cur + 1]; k++) {
                const Col *g = &out[padj[k]];
                for (long x = 0; x < W; x++) { int32_t v = colv(g, lo + x); if (v < in[x]) in[x] = v; }
            }
            for (long x = 0; x < W; x++) { val[x] = in[x]; org[x] = (int32_t)(lo + x); }
            for (long k = 0; k < L; k++) {
                int ch = code((unsigned char)s[k]);
                nval[0] = val[0] < DINF ? val[0] + 1 : DINF; norg[0] = org[0];
                for (long x = 1; x < W; x++) {
                    long jj = lo + x;
                    int32_t v = val[x - 1] < DINF ? val[x - 1] + (qc[jj - 1] != ch) : DINF; int32_t o = org[x - 1];
                    int32_t d = val[x] < DINF ? val[x] + 1 : DINF; if (d < v) { v = d; o = org[x]; }
                    d = nval[x - 1] < DINF ? nval[x - 1] + 1 : DINF; if (d < v) { v = d; o = norg[x - 1]; }
                    nval[x] = v; norg[x] = o;
                }
                int32_t *sw = val; val = nval; nval = sw;
                sw = org; org = norg; norg = sw;
            }
            int32_t exitv = val[W - 1];
            long j0 = org[W - 1];
            if (exitv != colv(&out[cur], j)) { best = -4; break; }
            int32_t target = in[j0 - lo];
            if (src[cur] && target == (int32_t)j0) break;   /* the path starts here */
            int nxt = -1;                                   /* among optimal predecessors, the highest priority */
            for (int k = pstart[cur]; k < pstart[cur + 1]; k++) {
                int g = padj[k];
                if (colv(&out[g], j0) == target && (nxt < 0 || prio[g] > prio[nxt])) nxt = g;
            }
            if (nxt < 0) { best = -4; break; }
            cur = nxt; j = j0;
        }
        free(in); free(val); free(org); free(nval); free(norg);
        if (best >= 0) {
            for (int a = 0, b = np - 1; a < b; a++, b--) { int t = path_out[a]; path_out[a] = path_out[b]; path_out[b] = t; }
            *path_len = np;
        }
    }
    for (int h = 0; h < H; h++) free(out[h].v);
    free(out); free(c); free(tmp); free(qc);
    return best;
}
