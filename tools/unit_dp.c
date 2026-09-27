/* unit_dp.c -- dynamic programs for tools/realign_units.py (repeat-unit-aware realignment).
 *
 * Built on demand into tools/bin/libunitdp.so (cc -O3 -shared -fPIC) and called through ctypes.
 *
 * unit_align(): align one sequence to a model made of
 *     [free prefix] -> FL (linear, k1 bases) -> cycle of the motif (p bases, traversed 0..many times)
 *                   -> FR (linear, k2 bases) -> [free suffix]
 * with gap-affine scores (maximised). FL and FR are aligned globally (every position is matched
 * or deleted); the prefix / suffix of the sequence outside them is skipped at no cost when the
 * corresponding flag is set, otherwise the alignment is global in the sequence too. With k1 = k2 = 0
 * and both flags set this is a local "wrap-around" alignment of the sequence to the motif (Fischetti
 * et al. 1993), which is how the array is found in the reference.
 *
 * The cycle is entered at motif position 0 only (a unit starting later begins with deletions) unless
 * flag 4 (free entry: any position, used for the local search in the reference) is set; it may
 * be left after any position (the last unit may be partial), and may be bypassed (no array). A
 * chain of deletions around the cycle is resolved with the usual two passes per row (a full turn of
 * deletions never scores better than none, because every deletion costs).
 *
 * Path output (forward order), one entry per event: op[] and node[]
 *   op 0 prefix base skipped, 1 M (base aligned to node), 2 I (base inserted after node),
 *   3 D (node deleted), 4 suffix base skipped, 5 wrap (cycle p-1 -> 0: a new unit starts),
 *   6 entry (FL -> cycle: the first unit starts), 7 exit (cycle -> FR), 8 bypass (no array).
 *   node: L = 0..k1-1, cycle = k1..k1+p-1, R = k1+p..k1+p+k2-1 (-1 for markers / skips).
 * Returns the path length, -(needed length) when cap is too small, or INT32_MIN on allocation failure.
 *
 * edit_distance(): unit-cost Levenshtein distance (two rows), used for the symbol score matrix.
 */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define NEG (-1000000000)

static inline int subst(char a, char b, int ma, int mi, int nsc)
{
    if (a == 'N' || b == 'N')
        return nsc;
    return a == b ? ma : -mi;
}

static inline int max3c(int a, int b, int c, int *code)
{
    /* ties prefer M (0), then I (1), then D (2) */
    int best = a, k = 0;
    if (b > best) { best = b; k = 1; }
    if (c > best) { best = c; k = 2; }
    *code = k;
    return best;
}

int32_t unit_align(const char *seq, int32_t n,
                   const char *fl, int32_t k1,
                   const char *motif, int32_t p,
                   const char *fr, int32_t k2,
                   int32_t ma, int32_t mi, int32_t go, int32_t ge, int32_t nsc,
                   int32_t flags,
                   uint8_t *op_out, int32_t *node_out, int32_t cap, int32_t *score_out)
{
    const int free_prefix = flags & 1, free_suffix = (flags >> 1) & 1, free_entry = (flags >> 2) & 1;
    const int m = k1 + p + k2;
    const int C0 = k1, R0 = k1 + p;
    int i, j, v;
    if (m <= 0 || p <= 0)
        return INT32_MIN;
    size_t W = (size_t)m;
    int *Mp = malloc(sizeof(int) * W), *Ip = malloc(sizeof(int) * W), *Dp = malloc(sizeof(int) * W);
    int *Mc = malloc(sizeof(int) * W), *Ic = malloc(sizeof(int) * W), *Dc = malloc(sizeof(int) * W);
    uint8_t *tb = malloc((size_t)(n + 1) * W);
    uint8_t *j1arg = malloc((size_t)(n + 1));
    int32_t *j2arg = malloc(sizeof(int32_t) * (size_t)(n + 1));
    uint8_t *endarg = malloc((size_t)(n + 1));
    int *J1 = malloc(sizeof(int) * (size_t)(n + 1)), *J2 = malloc(sizeof(int) * (size_t)(n + 1));
    int *END = malloc(sizeof(int) * (size_t)(n + 1));
    if (!Mp || !Ip || !Dp || !Mc || !Ic || !Dc || !tb || !j1arg || !j2arg || !endarg || !J1 || !J2 || !END) {
        free(Mp); free(Ip); free(Dp); free(Mc); free(Ic); free(Dc); free(tb); free(j1arg); free(j2arg);
        free(endarg); free(J1); free(J2); free(END);
        return INT32_MIN;
    }
    for (i = 0; i <= n; i++) {
        const int START = (i == 0 || free_prefix) ? 0 : NEG;
        const char b = i > 0 ? seq[i - 1] : 0;
        uint8_t *t = tb + (size_t)i * W;
        int c1, c2, c3, best;
        /* ---- FL */
        for (v = 0; v < k1; v++) {
            int mc = 3, ic = 0, dc = 3;
            if (i > 0) {
                if (v == 0) {
                    best = (i - 1 == 0 || free_prefix) ? 0 : NEG;
                    mc = 3;
                } else {
                    best = max3c(Mp[v - 1], Ip[v - 1], Dp[v - 1], &mc);
                }
                Mc[v] = best > NEG / 2 ? best + subst(b, fl[v], ma, mi, nsc) : NEG;
                Ic[v] = max3c(Mp[v] - go, Ip[v] - ge, Dp[v] - go, &ic);
                if (Ic[v] < NEG / 2) Ic[v] = NEG;
            } else {
                Mc[v] = NEG; Ic[v] = NEG;
            }
            if (v == 0) {
                Dc[v] = START > NEG / 2 ? START - go : NEG;
                dc = 3;
            } else {
                Dc[v] = max3c(Mc[v - 1] - go, Ic[v - 1] - go, Dc[v - 1] - ge, &c1);
                dc = c1;
                if (Dc[v] < NEG / 2) Dc[v] = NEG;
            }
            t[v] = (uint8_t)(mc | (ic << 2) | (dc << 4));
        }
        if (k1 > 0) {
            J1[i] = max3c(Mc[k1 - 1], Ic[k1 - 1], Dc[k1 - 1], &c1);
            j1arg[i] = (uint8_t)c1;
        } else {
            J1[i] = START;
            j1arg[i] = 3;
        }
        /* ---- cycle: M and I from the previous row */
        for (j = 0; j < p; j++) {
            v = C0 + j;
            int mc = 0, ic = 0;
            if (i > 0) {
                if (j == 0) {
                    best = max3c(Mp[C0 + p - 1], Ip[C0 + p - 1], Dp[C0 + p - 1], &mc);
                    if (J1[i - 1] >= best) { best = J1[i - 1]; mc = 3; }   /* entry preferred on ties */
                } else {
                    best = max3c(Mp[v - 1], Ip[v - 1], Dp[v - 1], &mc);
                    if (free_entry && J1[i - 1] > best) { best = J1[i - 1]; mc = 3; }
                }
                Mc[v] = best > NEG / 2 ? best + subst(b, motif[j], ma, mi, nsc) : NEG;
                Ic[v] = max3c(Mp[v] - go, Ip[v] - ge, Dp[v] - go, &ic);
                if (Ic[v] < NEG / 2) Ic[v] = NEG;
            } else {
                Mc[v] = NEG; Ic[v] = NEG;
            }
            t[v] = (uint8_t)(mc | (ic << 2));
        }
        /* ---- cycle: D, two passes around the ring */
        {
            int dc;
            /* pass 1: D[0] from entry and from M/I at p-1 (D at p-1 not known yet) */
            v = C0;
            best = J1[i] > NEG / 2 ? J1[i] - go : NEG;
            dc = 3;
            c2 = Mc[C0 + p - 1] - go;
            if (c2 > best) { best = c2; dc = 0; }
            c3 = Ic[C0 + p - 1] - go;
            if (c3 > best) { best = c3; dc = 1; }
            Dc[v] = best < NEG / 2 ? NEG : best;
            t[v] = (uint8_t)((t[v] & 0x0f) | (dc << 4));
            for (j = 1; j < p; j++) {
                v = C0 + j;
                Dc[v] = max3c(Mc[v - 1] - go, Ic[v - 1] - go, Dc[v - 1] - ge, &c1);
                if (Dc[v] < NEG / 2) Dc[v] = NEG;
                t[v] = (uint8_t)((t[v] & 0x0f) | (c1 << 4));
            }
            /* pass 2: D at p-1 feeds D at 0; propagate while it improves */
            if (p > 1 && Dc[C0 + p - 1] > NEG / 2) {
                int cand = Dc[C0 + p - 1] - ge;
                if (cand > Dc[C0]) {
                    Dc[C0] = cand;
                    t[C0] = (uint8_t)((t[C0] & 0x0f) | (2 << 4));
                    for (j = 1; j < p; j++) {
                        v = C0 + j;
                        cand = Dc[v - 1] - ge;
                        if (cand > Dc[v]) {
                            Dc[v] = cand;
                            t[v] = (uint8_t)((t[v] & 0x0f) | (2 << 4));
                        } else {
                            break;
                        }
                    }
                }
            }
        }
        /* ---- J2: bypass or exit from M / I of any cycle position */
        best = J1[i];
        {
            int arg = -1;
            for (j = 0; j < p; j++) {
                v = C0 + j;
                if (Mc[v] > best) { best = Mc[v]; arg = j * 4 + 0; }
                if (Ic[v] > best) { best = Ic[v]; arg = j * 4 + 1; }
            }
            J2[i] = best;
            j2arg[i] = arg;
        }
        /* ---- FR */
        for (int k = 0; k < k2; k++) {
            v = R0 + k;
            int mc = 3, ic = 0, dc = 3;
            if (i > 0) {
                if (k == 0) {
                    best = J2[i - 1];
                    mc = 3;
                } else {
                    best = max3c(Mp[v - 1], Ip[v - 1], Dp[v - 1], &mc);
                }
                Mc[v] = best > NEG / 2 ? best + subst(b, fr[k], ma, mi, nsc) : NEG;
                Ic[v] = max3c(Mp[v] - go, Ip[v] - ge, Dp[v] - go, &ic);
                if (Ic[v] < NEG / 2) Ic[v] = NEG;
            } else {
                Mc[v] = NEG; Ic[v] = NEG;
            }
            if (k == 0) {
                Dc[v] = J2[i] > NEG / 2 ? J2[i] - go : NEG;
                dc = 3;
            } else {
                Dc[v] = max3c(Mc[v - 1] - go, Ic[v - 1] - go, Dc[v - 1] - ge, &c1);
                dc = c1;
                if (Dc[v] < NEG / 2) Dc[v] = NEG;
            }
            t[v] = (uint8_t)(mc | (ic << 2) | (dc << 4));
        }
        if (k2 > 0) {
            END[i] = max3c(Mc[R0 + k2 - 1], Ic[R0 + k2 - 1], Dc[R0 + k2 - 1], &c1);
            endarg[i] = (uint8_t)c1;
        } else {
            END[i] = J2[i];
            endarg[i] = 3;
        }
        /* roll rows */
        int *x;
        x = Mp; Mp = Mc; Mc = x;
        x = Ip; Ip = Ic; Ic = x;
        x = Dp; Dp = Dc; Dc = x;
    }
    /* ---- best end row */
    int iend = n, bestend = END[n];
    if (free_suffix) {
        for (i = n; i >= 0; i--)
            if (END[i] > bestend) { bestend = END[i]; iend = i; }
    }
    *score_out = bestend;
    /* ---- traceback (reverse), then flip */
    int32_t len = 0, need = 0;
    /* first pass counts, second pass writes: simpler to write into a growing local buffer */
    size_t bcap = (size_t)n + (size_t)m + 64, blen = 0;
    uint8_t *bop = malloc(bcap);
    int32_t *bnode = malloc(sizeof(int32_t) * bcap);
    if (!bop || !bnode) {
        free(bop); free(bnode);
        len = INT32_MIN;
        goto done;
    }
#define EMIT(o, nd) do { if (blen == bcap) { bcap *= 2; bop = realloc(bop, bcap); bnode = realloc(bnode, sizeof(int32_t) * bcap); } \
                         bop[blen] = (uint8_t)(o); bnode[blen] = (nd); blen++; } while (0)
    for (i = n; i > iend; i--)
        EMIT(4, -1);
    {
        /* where: 0 = node state, 1 = START, 2 = J1, 3 = J2, 4 = END */
        int where = 4, st = 0;
        i = iend;
        v = 0;
        while (1) {
            if (where == 4) {
                if (k2 > 0) { where = 0; v = R0 + k2 - 1; st = endarg[i]; }
                else where = 3;
                continue;
            }
            if (where == 3) {
                int arg = j2arg[i];
                if (arg < 0) { EMIT(8, -1); where = 2; }
                else { EMIT(7, -1); where = 0; v = C0 + arg / 4; st = arg % 4; }
                continue;
            }
            if (where == 2) {
                if (k1 > 0) { where = 0; v = k1 - 1; st = j1arg[i]; }
                else where = 1;
                continue;
            }
            if (where == 1) {
                for (; i > 0; i--)
                    EMIT(0, -1);
                break;
            }
            /* node state */
            uint8_t code = tb[(size_t)i * W + v];
            int c;
            if (st == 0) {            /* M: consumes base i-1 */
                EMIT(1, v);
                c = code & 3;
                i -= 1;
                if (c == 3) {
                    if (v == 0 && k1 > 0) where = 1;
                    else if (v >= C0 && v < C0 + p) { EMIT(6, -1); where = 2; }   /* entry (any position if free) */
                    else if (v == R0) where = 3;
                    continue;
                }
                if (v == C0) { EMIT(5, -1); v = C0 + p - 1; }
                else v = v - 1;
                st = c;
            } else if (st == 1) {     /* I: consumes base i-1, stays on v */
                EMIT(2, v);
                c = (code >> 2) & 3;
                i -= 1;
                st = c;
            } else {                  /* D: no base */
                EMIT(3, v);
                c = (code >> 4) & 3;
                if (c == 3) {
                    if (v == 0 && k1 > 0) where = 1;
                    else if (v == C0) { EMIT(6, -1); where = 2; }
                    else if (v == R0) where = 3;
                    continue;
                }
                if (v == C0) { EMIT(5, -1); v = C0 + p - 1; }
                else v = v - 1;
                st = c;
            }
            if (i < 0) { len = INT32_MIN; goto done; }
        }
    }
    need = (int32_t)blen;
    if (need > cap) {
        len = -need;
    } else {
        for (size_t q = 0; q < blen; q++) {
            op_out[q] = bop[blen - 1 - q];
            node_out[q] = bnode[blen - 1 - q];
        }
        len = need;
    }
done:
    free(bop); free(bnode);
    free(Mp); free(Ip); free(Dp); free(Mc); free(Ic); free(Dc); free(tb); free(j1arg); free(j2arg);
    free(endarg); free(J1); free(J2); free(END);
    return len;
}

int32_t edit_distance(const char *a, int32_t la, const char *b, int32_t lb)
{
    if (la == 0) return lb;
    if (lb == 0) return la;
    int32_t *prev = malloc(sizeof(int32_t) * (size_t)(lb + 1)), *cur = malloc(sizeof(int32_t) * (size_t)(lb + 1));
    if (!prev || !cur) { free(prev); free(cur); return -1; }
    for (int j = 0; j <= lb; j++) prev[j] = j;
    for (int i = 1; i <= la; i++) {
        cur[0] = i;
        for (int j = 1; j <= lb; j++) {
            int32_t x = prev[j - 1] + (a[i - 1] != b[j - 1]);
            int32_t y = prev[j] + 1, z = cur[j - 1] + 1;
            if (y < x) x = y;
            if (z < x) x = z;
            cur[j] = x;
        }
        int32_t *t = prev; prev = cur; cur = t;
    }
    int32_t r = prev[lb];
    free(prev); free(cur);
    return r;
}

/* All-pairs distances among ns strings (concatenated in buf, offsets off[ns+1]) into out[ns*ns]. */
void edit_matrix(const char *buf, const int32_t *off, int32_t ns, int32_t *out)
{
    for (int a = 0; a < ns; a++) {
        out[(size_t)a * ns + a] = 0;
        for (int b = a + 1; b < ns; b++) {
            int32_t d = edit_distance(buf + off[a], off[a + 1] - off[a], buf + off[b], off[b + 1] - off[b]);
            out[(size_t)a * ns + b] = d;
            out[(size_t)b * ns + a] = d;
        }
    }
}

/* Distances from each of nq query strings to each of nr reference strings: out[nq*nr]. */
void edit_cross(const char *qbuf, const int32_t *qoff, int32_t nq,
                const char *rbuf, const int32_t *roff, int32_t nr, int32_t *out)
{
    for (int a = 0; a < nq; a++)
        for (int b = 0; b < nr; b++)
            out[(size_t)a * nr + b] = edit_distance(qbuf + qoff[a], qoff[a + 1] - qoff[a],
                                                    rbuf + roff[b], roff[b + 1] - roff[b]);
}

/* profile_align(): optimal alignment of two unit-level profiles under the sum-of-pairs cost with
 * a symbol cost matrix C ((S+1) x (S+1), index S = gap; C[S][S] = 0) -- linear gap costs, so the
 * sum of pairs of the merged alignment is exact. A profile of L columns is given sparsely: column
 * c holds entries off[c]..off[c+1]-1 of (sym, weight). w is the profile's total row weight (a gap
 * column inserted into a profile gaps every row). ops_out (forward): 0 = column of A with column
 * of B, 1 = column of A against a new gap column in B, 2 = gap column in A against column of B.
 * Returns the cost of the merge (pairs across A and B only); *nops gets the path length (<= la+lb).
 * Ties prefer the diagonal, then a gap column in B. */
double profile_align(int32_t la, const int32_t *offA, const int32_t *symA, const double *wtA, double wA,
                     int32_t lb, const int32_t *offB, const int32_t *symB, const double *wtB, double wB,
                     const double *C, int32_t S1, uint8_t *ops_out, int32_t *nops)
{
    const int32_t G = S1 - 1;
    size_t W = (size_t)lb + 1;
    double *prev = malloc(sizeof(double) * W), *cur = malloc(sizeof(double) * W);
    double *gapA = malloc(sizeof(double) * (size_t)(lb + 1)), *gapB = malloc(sizeof(double) * (size_t)(la + 1));
    uint8_t *tb = malloc((size_t)(la + 1) * W);
    if (!prev || !cur || !gapA || !gapB || !tb) {
        free(prev); free(cur); free(gapA); free(gapB); free(tb);
        *nops = -1;
        return -1.0;
    }
    for (int i = 0; i < la; i++) {           /* column i of A against an all-gap column of B */
        double s = 0;
        for (int e = offA[i]; e < offA[i + 1]; e++)
            s += wtA[e] * C[(size_t)symA[e] * S1 + G];
        gapB[i] = s * wB;
    }
    for (int j = 0; j < lb; j++) {
        double s = 0;
        for (int e = offB[j]; e < offB[j + 1]; e++)
            s += wtB[e] * C[(size_t)G * S1 + symB[e]];
        gapA[j] = s * wA;
    }
    prev[0] = 0;
    tb[0] = 0;
    for (int j = 1; j <= lb; j++) {
        prev[j] = prev[j - 1] + gapA[j - 1];
        tb[j] = 2;
    }
    for (int i = 1; i <= la; i++) {
        uint8_t *t = tb + (size_t)i * W;
        cur[0] = prev[0] + gapB[i - 1];
        t[0] = 1;
        const int ea0 = offA[i - 1], ea1 = offA[i];
        for (int j = 1; j <= lb; j++) {
            double s = 0;
            for (int e = ea0; e < ea1; e++) {
                const double *Crow = C + (size_t)symA[e] * S1;
                double x = 0;
                for (int f = offB[j - 1]; f < offB[j]; f++)
                    x += wtB[f] * Crow[symB[f]];
                s += wtA[e] * x;
            }
            double d = prev[j - 1] + s, u = prev[j] + gapB[i - 1], l = cur[j - 1] + gapA[j - 1];
            double best = d;
            uint8_t k = 0;
            const double eps = 1e-9;
            if (u < best - eps) { best = u; k = 1; }
            if (l < best - eps) { best = l; k = 2; }
            cur[j] = best;
            t[j] = k;
        }
        double *x = prev; prev = cur; cur = x;
    }
    double total = prev[lb];
    /* traceback */
    int i = la, j = lb, n = 0;
    while (i > 0 || j > 0) {
        uint8_t k = tb[(size_t)i * W + j];
        ops_out[n++] = k;
        if (k == 0) { i--; j--; }
        else if (k == 1) i--;
        else j--;
    }
    for (int a = 0, b = n - 1; a < b; a++, b--) {
        uint8_t x = ops_out[a]; ops_out[a] = ops_out[b]; ops_out[b] = x;
    }
    *nops = n;
    free(prev); free(cur); free(gapA); free(gapB); free(tb);
    return total;
}

/* Pairwise unit-level distances among n symbol strings (seqs concatenated, offsets off[n+1]),
 * linear-gap DP with cost matrix C (as profile_align with one row each); out[n*n]. */
void unit_pair_costs(const int32_t *seqs, const int32_t *off, int32_t n, const double *C, int32_t S1, double *out)
{
    const int32_t G = S1 - 1;
    int32_t maxl = 0;
    for (int a = 0; a < n; a++)
        if (off[a + 1] - off[a] > maxl) maxl = off[a + 1] - off[a];
    double *prev = malloc(sizeof(double) * (size_t)(maxl + 1)), *cur = malloc(sizeof(double) * (size_t)(maxl + 1));
    for (int a = 0; a < n; a++) {
        out[(size_t)a * n + a] = 0;
        const int32_t *A = seqs + off[a];
        int la = off[a + 1] - off[a];
        for (int b = a + 1; b < n; b++) {
            const int32_t *B = seqs + off[b];
            int lb = off[b + 1] - off[b];
            prev[0] = 0;
            for (int j = 1; j <= lb; j++) prev[j] = prev[j - 1] + C[(size_t)G * S1 + B[j - 1]];
            for (int i = 1; i <= la; i++) {
                cur[0] = prev[0] + C[(size_t)A[i - 1] * S1 + G];
                const double *Crow = C + (size_t)A[i - 1] * S1;
                const double gi = Crow[G];
                for (int j = 1; j <= lb; j++) {
                    double d = prev[j - 1] + Crow[B[j - 1]];
                    double u = prev[j] + gi;
                    double l = cur[j - 1] + C[(size_t)G * S1 + B[j - 1]];
                    cur[j] = d < u ? (d < l ? d : l) : (u < l ? u : l);
                }
                double *x = prev; prev = cur; cur = x;
            }
            out[(size_t)a * n + b] = out[(size_t)b * n + a] = prev[lb];
        }
    }
    free(prev); free(cur);
}

/* Guide-tree distances for many strings: L1 distance between hashed symbol k-mer count vectors
 * (4096 bins), normalised by the summed k-mer counts; out[n*n]. */
void kmer_dists(const int32_t *seqs, const int32_t *off, int32_t n, int32_t k, double *out)
{
    const int B = 4096;
    float *v = calloc((size_t)n * B, sizeof(float));
    double *tot = calloc((size_t)n, sizeof(double));
    if (!v || !tot) { free(v); free(tot); return; }
    for (int a = 0; a < n; a++) {
        const int32_t *s = seqs + off[a];
        int L = off[a + 1] - off[a];
        for (int i = 0; i + k <= L; i++) {
            uint32_t h = 2166136261u;
            for (int q = 0; q < k; q++) { h ^= (uint32_t)s[i + q] + 1; h *= 16777619u; }
            v[(size_t)a * B + (h % B)] += 1;
            tot[a] += 1;
        }
    }
    for (int a = 0; a < n; a++) {
        out[(size_t)a * n + a] = 0;
        for (int b = a + 1; b < n; b++) {
            double d = 0;
            const float *x = v + (size_t)a * B, *y = v + (size_t)b * B;
            for (int q = 0; q < B; q++) { double t = x[q] - y[q]; d += t < 0 ? -t : t; }
            double den = tot[a] + tot[b];
            out[(size_t)a * n + b] = out[(size_t)b * n + a] = den > 0 ? d / den : 0;
        }
    }
    free(v); free(tot);
}

/* profile_align_affine(): as profile_align, plus a gap-open cost go per pair of rows whenever a
 * run of new gap columns opens in one profile against residues of the other (the usual
 * approximation: every row of the gapped profile is taken to open a gap, weighted by the
 * residue weight of the other profile's column). go = 0 gives profile_align's result. */
double profile_align_affine(int32_t la, const int32_t *offA, const int32_t *symA, const double *wtA, double wA,
                            int32_t lb, const int32_t *offB, const int32_t *symB, const double *wtB, double wB,
                            const double *C, int32_t S1, double go, uint8_t *ops_out, int32_t *nops)
{
    const int32_t G = S1 - 1;
    const double INF = 1e300;
    size_t W = (size_t)lb + 1;
    double *M0 = malloc(sizeof(double) * W), *X0 = malloc(sizeof(double) * W), *Y0 = malloc(sizeof(double) * W);
    double *M1 = malloc(sizeof(double) * W), *X1 = malloc(sizeof(double) * W), *Y1 = malloc(sizeof(double) * W);
    double *gapA = malloc(sizeof(double) * (size_t)(lb + 1)), *gapB = malloc(sizeof(double) * (size_t)(la + 1));
    double *openA = malloc(sizeof(double) * (size_t)(lb + 1)), *openB = malloc(sizeof(double) * (size_t)(la + 1));
    uint8_t *tb = malloc((size_t)(la + 1) * W);
    if (!M0 || !X0 || !Y0 || !M1 || !X1 || !Y1 || !gapA || !gapB || !openA || !openB || !tb) {
        free(M0); free(X0); free(Y0); free(M1); free(X1); free(Y1); free(gapA); free(gapB); free(openA);
        free(openB); free(tb);
        *nops = -1;
        return -1.0;
    }
    for (int i = 0; i < la; i++) {
        double s = 0, res = 0;
        for (int e = offA[i]; e < offA[i + 1]; e++) {
            s += wtA[e] * C[(size_t)symA[e] * S1 + G];
            if (symA[e] != G) res += wtA[e];
        }
        gapB[i] = s * wB;
        openB[i] = go * res * wB;
    }
    for (int j = 0; j < lb; j++) {
        double s = 0, res = 0;
        for (int e = offB[j]; e < offB[j + 1]; e++) {
            s += wtB[e] * C[(size_t)G * S1 + symB[e]];
            if (symB[e] != G) res += wtB[e];
        }
        gapA[j] = s * wA;
        openA[j] = go * res * wA;
    }
    /* tb bits: M pred (0-1), X pred (2-3), Y pred (4-5); pred codes 0 = M, 1 = X, 2 = Y */
    M0[0] = 0; X0[0] = INF; Y0[0] = INF;
    tb[0] = 0;
    for (int j = 1; j <= lb; j++) {
        M0[j] = INF; X0[j] = INF;
        Y0[j] = (j == 1 ? M0[0] + openA[0] : Y0[j - 1]) + gapA[j - 1];
        tb[j] = (uint8_t)((j == 1 ? 0 : 2) << 4);
    }
    for (int i = 1; i <= la; i++) {
        uint8_t *t = tb + (size_t)i * W;
        M1[0] = INF; Y1[0] = INF;
        X1[0] = (i == 1 ? M0[0] + openB[0] : X0[0]) + gapB[i - 1];
        t[0] = (uint8_t)((i == 1 ? 0 : 1) << 2);
        const int ea0 = offA[i - 1], ea1 = offA[i];
        const double gb = gapB[i - 1], ob = openB[i - 1];
        for (int j = 1; j <= lb; j++) {
            double s = 0;
            for (int e = ea0; e < ea1; e++) {
                const double *Crow = C + (size_t)symA[e] * S1;
                double x = 0;
                for (int f = offB[j - 1]; f < offB[j]; f++)
                    x += wtB[f] * Crow[symB[f]];
                s += wtA[e] * x;
            }
            int cm = 0, cx = 0, cy = 0;
            double v = M0[j - 1];
            if (X0[j - 1] < v - 1e-9) { v = X0[j - 1]; cm = 1; }
            if (Y0[j - 1] < v - 1e-9) { v = Y0[j - 1]; cm = 2; }
            M1[j] = v + s;
            /* X: A column i against a gap column in B (comes from row i-1, same j) */
            double vx = M0[j] + ob;
            if (X0[j] < vx - 1e-9) { vx = X0[j]; cx = 1; }
            if (Y0[j] + ob < vx - 1e-9) { vx = Y0[j] + ob; cx = 2; }
            X1[j] = vx + gb;
            /* Y: gap column in A against B column j (same row i, j-1) */
            const double oa = openA[j - 1];
            double vy = M1[j - 1] + oa;
            if (Y1[j - 1] < vy - 1e-9) { vy = Y1[j - 1]; cy = 2; }
            if (X1[j - 1] + oa < vy - 1e-9) { vy = X1[j - 1] + oa; cy = 1; }
            Y1[j] = vy + gapA[j - 1];
            t[j] = (uint8_t)(cm | (cx << 2) | (cy << 4));
        }
        double *x;
        x = M0; M0 = M1; M1 = x;
        x = X0; X0 = X1; X1 = x;
        x = Y0; Y0 = Y1; Y1 = x;
    }
    int st = 0;
    double total = M0[lb];
    if (X0[lb] < total - 1e-9) { total = X0[lb]; st = 1; }
    if (Y0[lb] < total - 1e-9) { total = Y0[lb]; st = 2; }
    int i = la, j = lb, n = 0;
    while (i > 0 || j > 0) {
        uint8_t code = tb[(size_t)i * W + j];
        if (st == 0) {
            ops_out[n++] = 0;
            st = code & 3;
            i--; j--;
        } else if (st == 1) {
            ops_out[n++] = 1;
            st = (code >> 2) & 3;
            i--;
        } else {
            ops_out[n++] = 2;
            st = (code >> 4) & 3;
            j--;
        }
        if (i < 0 || j < 0) { n = -1; break; }
    }
    if (n >= 0)
        for (int a = 0, b = n - 1; a < b; a++, b--) {
            uint8_t x = ops_out[a]; ops_out[a] = ops_out[b]; ops_out[b] = x;
        }
    *nops = n;
    free(M0); free(X0); free(Y0); free(M1); free(X1); free(Y1); free(gapA); free(gapB); free(openA); free(openB);
    free(tb);
    return total;
}

/* Exact weighted sum-of-pairs cost of an MSA of symbol rows (rows[n*L], -1 = gap): per pair,
 * the columns where both are gaps are skipped, every other column costs C[a][b] (gap = G), and
 * each run of gaps in one row against residues of the other costs go more. */
double sp_affine(const int32_t *rows, int32_t n, int32_t L, const double *w, const double *C, int32_t S1, double go)
{
    const int32_t G = S1 - 1;
    double tot = 0;
    for (int a = 0; a < n; a++) {
        const int32_t *ra = rows + (size_t)a * L;
        for (int b = a + 1; b < n; b++) {
            const int32_t *rb = rows + (size_t)b * L;
            double c = 0;
            int opens = 0, st = 0;
            for (int k = 0; k < L; k++) {
                int x = ra[k], y = rb[k];
                if (x < 0 && y < 0) continue;
                if (x < 0) {
                    c += C[(size_t)G * S1 + y];
                    if (st != 1) { opens++; st = 1; }
                } else if (y < 0) {
                    c += C[(size_t)x * S1 + G];
                    if (st != 2) { opens++; st = 2; }
                } else {
                    c += C[(size_t)x * S1 + y];
                    st = 0;
                }
            }
            tot += w[a] * w[b] * (c + go * opens);
        }
    }
    return tot;
}

/* ---- flat-row helpers for the progressive unit aligner: an MSA is n rows x L columns of int32
 * symbol ids (-1 = gap), row-major. */

/* sparse column profile of a flat MSA: off[L+1], sym[], wt[] (gap = G); returns entries used */
static int32_t build_profile(const int32_t *R, int32_t n, int32_t L, const double *w, int32_t G,
                             int32_t *off, int32_t *sym, double *wt, double *scratch, int32_t *touched)
{
    int32_t e = 0;
    for (int c = 0; c < L; c++) {
        off[c] = e;
        int nt = 0;
        for (int r = 0; r < n; r++) {
            int x = R[(size_t)r * L + c];
            if (x < 0) x = G;
            if (scratch[x] == 0) touched[nt++] = x;
            scratch[x] += w[r];
        }
        for (int t = 0; t < nt; t++) {
            sym[e] = touched[t];
            wt[e] = scratch[touched[t]];
            scratch[touched[t]] = 0;
            e++;
        }
    }
    off[L] = e;
    return e;
}

/* Align flat MSAs A (nA x La) and B (nB x Lb); ops_out needs La + Lb entries. Returns the path
 * length (or -1); *cost_out = merge cost from the DP (approximate when go > 0). */
int32_t profile_align_flat(const int32_t *A, int32_t nA, int32_t La, const double *wA,
                           const int32_t *B, int32_t nB, int32_t Lb, const double *wB,
                           const double *C, int32_t S1, double go, uint8_t *ops_out, double *cost_out)
{
    const int32_t G = S1 - 1;
    int32_t *offA = malloc(sizeof(int32_t) * (size_t)(La + 1)), *offB = malloc(sizeof(int32_t) * (size_t)(Lb + 1));
    size_t eA = (size_t)La * (nA < S1 ? nA : S1) + 1, eB = (size_t)Lb * (nB < S1 ? nB : S1) + 1;
    int32_t *symA = malloc(sizeof(int32_t) * eA), *symB = malloc(sizeof(int32_t) * eB);
    double *wtA = malloc(sizeof(double) * eA), *wtB = malloc(sizeof(double) * eB);
    double *scratch = calloc((size_t)S1, sizeof(double));
    int32_t *touched = malloc(sizeof(int32_t) * (size_t)S1);
    int32_t n = -1;
    if (!offA || !offB || !symA || !symB || !wtA || !wtB || !scratch || !touched)
        goto out;
    build_profile(A, nA, La, wA, G, offA, symA, wtA, scratch, touched);
    build_profile(B, nB, Lb, wB, G, offB, symB, wtB, scratch, touched);
    double sa = 0, sb = 0;
    for (int r = 0; r < nA; r++) sa += wA[r];
    for (int r = 0; r < nB; r++) sb += wB[r];
    if (go > 0)
        *cost_out = profile_align_affine(La, offA, symA, wtA, sa, Lb, offB, symB, wtB, sb, C, S1, go, ops_out, &n);
    else
        *cost_out = profile_align(La, offA, symA, wtA, sa, Lb, offB, symB, wtB, sb, C, S1, ops_out, &n);
out:
    free(offA); free(offB); free(symA); free(symB); free(wtA); free(wtB); free(scratch); free(touched);
    return n;
}

/* Merged flat MSA (A rows then B rows, nops columns) from an ops path. */
void merge_flat(const int32_t *A, int32_t nA, int32_t La, const int32_t *B, int32_t nB, int32_t Lb,
                const uint8_t *ops, int32_t nops, int32_t *out)
{
    int i = 0, j = 0;
    for (int k = 0; k < nops; k++) {
        const uint8_t o = ops[k];
        for (int r = 0; r < nA; r++)
            out[(size_t)r * nops + k] = (o != 2) ? A[(size_t)r * La + i] : -1;
        for (int r = 0; r < nB; r++)
            out[(size_t)(nA + r) * nops + k] = (o != 1) ? B[(size_t)r * Lb + j] : -1;
        if (o != 2) i++;
        if (o != 1) j++;
    }
}

/* Drop columns that are gaps in every row; returns the new column count (out: n x newL). */
int32_t drop_gap_cols(const int32_t *R, int32_t n, int32_t L, int32_t *out)
{
    int32_t k = 0;
    for (int c = 0; c < L; c++) {
        int any = 0;
        for (int r = 0; r < n && !any; r++)
            if (R[(size_t)r * L + c] >= 0) any = 1;
        if (any) k++;
    }
    int32_t q = 0;
    for (int c = 0; c < L; c++) {
        int any = 0;
        for (int r = 0; r < n && !any; r++)
            if (R[(size_t)r * L + c] >= 0) any = 1;
        if (!any) continue;
        for (int r = 0; r < n; r++)
            out[(size_t)r * k + q] = R[(size_t)r * L + c];
        q++;
    }
    return k;
}

/* Sum over pairs (a in rows ia, b in rows ib) of an MSA (n x L) of the exact pair cost
 * (as sp_affine). */
double cross_cost(const int32_t *R, int32_t L, const int32_t *ia, int32_t na, const int32_t *ib, int32_t nb,
                  const double *w, const double *C, int32_t S1, double go)
{
    const int32_t G = S1 - 1;
    double tot = 0;
    for (int p = 0; p < na; p++) {
        const int32_t *ra = R + (size_t)ia[p] * L;
        for (int q = 0; q < nb; q++) {
            const int32_t *rb = R + (size_t)ib[q] * L;
            double c = 0;
            int opens = 0, st = 0;
            for (int k = 0; k < L; k++) {
                int x = ra[k], y = rb[k];
                if (x < 0 && y < 0) continue;
                if (x < 0) {
                    c += C[(size_t)G * S1 + y];
                    if (st != 1) { opens++; st = 1; }
                } else if (y < 0) {
                    c += C[(size_t)x * S1 + G];
                    if (st != 2) { opens++; st = 2; }
                } else {
                    c += C[(size_t)x * S1 + y];
                    st = 0;
                }
            }
            tot += w[ia[p]] * w[ib[q]] * (c + go * opens);
        }
    }
    return tot;
}

/* row_to_profile_banded(): align one ungapped row r[0..m-1] (symbols) to a profile (flat MSA
 * P, n x L, weights w) within a band: after profile column i (i = 0..L) the row prefix length j
 * must lie in [lo[i], hi[i]] (both non-decreasing, lo[0] = 0 allowed, hi[L] = m required).
 * Gotoh on the profile with the gap-open approximation of profile_align_affine (go per pair).
 * ops_out (L + m entries): 0 = profile column with row base, 1 = profile column with a gap in the
 * row, 2 = row base in a new column (gaps in the profile). Returns the path length or -1. */
int32_t row_to_profile_banded(const int32_t *P, int32_t n, int32_t L, const double *w,
                              const int32_t *r, int32_t m, const int32_t *lo, const int32_t *hi,
                              const double *C, int32_t S1, double go, uint8_t *ops_out, double *cost_out)
{
    const int32_t G = S1 - 1;
    const double INF = 1e300;
    double W = 0;
    for (int k = 0; k < n; k++) W += w[k];
    /* per column: sparse profile entries and residue weight */
    int32_t *off = malloc(sizeof(int32_t) * (size_t)(L + 1));
    size_t cap = (size_t)L * (n < S1 ? n : S1) + 1;
    int32_t *sym = malloc(sizeof(int32_t) * cap);
    double *wt = malloc(sizeof(double) * cap);
    double *scratch = calloc((size_t)S1, sizeof(double));
    int32_t *touched = malloc(sizeof(int32_t) * (size_t)S1);
    size_t *boff = malloc(sizeof(size_t) * (size_t)(L + 2));
    double *colgap = malloc(sizeof(double) * (size_t)(L + 1)), *colres = malloc(sizeof(double) * (size_t)(L + 1));
    if (!off || !sym || !wt || !scratch || !touched || !boff || !colgap || !colres) {
        free(off); free(sym); free(wt); free(scratch); free(touched); free(boff); free(colgap); free(colres);
        return -1;
    }
    build_profile(P, n, L, w, G, off, sym, wt, scratch, touched);
    for (int i = 0; i < L; i++) {
        double g = 0, res = 0;
        for (int e = off[i]; e < off[i + 1]; e++) {
            g += wt[e] * C[(size_t)sym[e] * S1 + G];
            if (sym[e] != G) res += wt[e];
        }
        colgap[i] = g;       /* row gapped against column i: sum over profile rows */
        colres[i] = res;
    }
    size_t tot = 0;
    for (int i = 0; i <= L; i++) {
        boff[i] = tot;
        if (hi[i] >= lo[i]) tot += (size_t)(hi[i] - lo[i] + 1);
    }
    boff[L + 1] = tot;
    double *M = malloc(sizeof(double) * tot), *X = malloc(sizeof(double) * tot), *Y = malloc(sizeof(double) * tot);
    uint8_t *tb = malloc(tot);
    if (!M || !X || !Y || !tb) {
        free(off); free(sym); free(wt); free(scratch); free(touched); free(boff); free(colgap); free(colres);
        free(M); free(X); free(Y); free(tb);
        return -1;
    }
#define CELL(i, j) (boff[i] + (size_t)((j) - lo[i]))
#define INB(i, j) ((j) >= lo[i] && (j) <= hi[i])
    /* M: column i with base j-1; X: column i with a gap in the row (from (i-1, j));
       Y: base j-1 in a new column (from (i, j-1)). Gap-open costs: X opens go * colres[i-1]
       (pairs row x profile rows with residues), Y opens go * W (all profile rows gapped). */
    for (int i = 0; i <= L; i++) {
        for (int j = lo[i]; j <= hi[i]; j++) {
            size_t c = CELL(i, j);
            double vm = INF, vx = INF, vy = INF;
            uint8_t cm = 0, cx = 0, cy = 0;
            if (i == 0 && j == 0) {
                M[c] = 0; X[c] = INF; Y[c] = INF; tb[c] = 0;
                continue;
            }
            if (i > 0 && j > 0 && INB(i - 1, j - 1)) {
                size_t p = CELL(i - 1, j - 1);
                double s = 0;
                const int32_t b = r[j - 1];
                for (int e = off[i - 1]; e < off[i]; e++)
                    s += wt[e] * C[(size_t)sym[e] * S1 + b];
                vm = M[p]; cm = 0;
                if (X[p] < vm - 1e-9) { vm = X[p]; cm = 1; }
                if (Y[p] < vm - 1e-9) { vm = Y[p]; cm = 2; }
                vm += s;
            }
            if (i > 0 && INB(i - 1, j)) {
                size_t p = CELL(i - 1, j);
                const double o = go * colres[i - 1];
                vx = M[p] + o; cx = 0;
                if (X[p] < vx - 1e-9) { vx = X[p]; cx = 1; }
                if (Y[p] + o < vx - 1e-9) { vx = Y[p] + o; cx = 2; }
                vx += colgap[i - 1];
            }
            if (j > 0 && INB(i, j - 1)) {
                size_t p = CELL(i, j - 1);
                const double o = go * W;
                vy = M[p] + o; cy = 0;
                if (Y[p] < vy - 1e-9) { vy = Y[p]; cy = 2; }
                if (X[p] + o < vy - 1e-9) { vy = X[p] + o; cy = 1; }
                vy += W * C[(size_t)G * S1 + r[j - 1]];
            }
            M[c] = vm; X[c] = vx; Y[c] = vy;
            tb[c] = (uint8_t)(cm | (cx << 2) | (cy << 4));
        }
    }
    int32_t nout = -1;
    if (INB(L, m)) {
        size_t c = CELL(L, m);
        int st = 0;
        double best = M[c];
        if (X[c] < best - 1e-9) { best = X[c]; st = 1; }
        if (Y[c] < best - 1e-9) { best = Y[c]; st = 2; }
        *cost_out = best;
        if (best < INF / 2) {
            int i = L, j = m, k = 0;
            while (i > 0 || j > 0) {
                uint8_t code = tb[CELL(i, j)];
                if (st == 0) { ops_out[k++] = 0; st = code & 3; i--; j--; }
                else if (st == 1) { ops_out[k++] = 1; st = (code >> 2) & 3; i--; }
                else { ops_out[k++] = 2; st = (code >> 4) & 3; j--; }
                if (i < 0 || j < 0 || !INB(i, j)) { k = -1; break; }
            }
            if (k >= 0) {
                for (int a = 0, b = k - 1; a < b; a++, b--) {
                    uint8_t x = ops_out[a]; ops_out[a] = ops_out[b]; ops_out[b] = x;
                }
            }
            nout = k;
        }
    }
#undef CELL
#undef INB
    free(off); free(sym); free(wt); free(scratch); free(touched); free(boff); free(colgap); free(colres);
    free(M); free(X); free(Y); free(tb);
    return nout;
}

/* For row x of a flat MSA (n x L): rest_out gets the other n-1 rows without the columns that are
 * gaps in all of them (returns that column count Lr); row_out the ungapped row (*m); lo/hi
 * (Lr + 1 entries) the band of half-width B around the row's current path (row prefix length
 * after each kept column). */
int32_t band_extract(const int32_t *R, int32_t n, int32_t L, int32_t x, int32_t B,
                     int32_t *rest_out, int32_t *row_out, int32_t *m_out, int32_t *lo, int32_t *hi)
{
    int32_t Lr = 0, m = 0;
    for (int c = 0; c < L; c++) {
        int any = 0;
        for (int r = 0; r < n && !any; r++)
            if (r != x && R[(size_t)r * L + c] >= 0) any = 1;
        if (any) Lr++;
    }
    for (int c = 0; c < L; c++)
        if (R[(size_t)x * L + c] >= 0) m++;
    int32_t q = 0, pre = 0;
    lo[0] = 0;
    hi[0] = 0;
    int32_t *pos = malloc(sizeof(int32_t) * (size_t)(Lr + 1));
    pos[0] = 0;
    for (int c = 0; c < L; c++) {
        int any = 0;
        for (int r = 0; r < n && !any; r++)
            if (r != x && R[(size_t)r * L + c] >= 0) any = 1;
        if (R[(size_t)x * L + c] >= 0)
            row_out[pre++] = R[(size_t)x * L + c];
        if (!any) continue;
        int k = 0;
        for (int r = 0; r < n; r++) {
            if (r == x) continue;
            rest_out[(size_t)k * Lr + q] = R[(size_t)r * L + c];
            k++;
        }
        q++;
        pos[q] = pre;
    }
    /* band: [pos - B, pos + B] clipped, and the start of column 0 / end of column Lr */
    for (int i = 0; i <= Lr; i++) {
        int a = pos[i] - B, b = pos[i] + B;
        lo[i] = a < 0 ? 0 : a;
        hi[i] = b > m ? m : b;
    }
    lo[0] = 0;
    hi[Lr] = m;
    if (lo[Lr] > m) lo[Lr] = m;
    for (int i = 1; i <= Lr; i++)       /* keep both bounds non-decreasing */
        if (lo[i] < lo[i - 1]) lo[i] = lo[i - 1];
    for (int i = Lr - 1; i >= 0; i--)
        if (hi[i] > hi[i + 1]) hi[i] = hi[i + 1];
    free(pos);
    *m_out = m;
    return Lr;
}

/* Unit-cost global edit distance, Myers/Hyyro bit-vector (the same algorithm as ed_unit in
 * tools/c/vntreval.c, copied so this library stands alone). */
static inline int bcode(unsigned char c)
{
    switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return 4;
    }
}

int64_t ed_bitpar(const char *p, int64_t m, const char *t, int64_t n)
{
    if (m == 0) return n;
    if (n == 0) return m;
    int64_t W = (m + 63) / 64;
    uint64_t *peq = calloc((size_t)(5 * W), sizeof(uint64_t));
    uint64_t *Pv = malloc((size_t)W * sizeof(uint64_t)), *Mv = calloc((size_t)W, sizeof(uint64_t));
    if (!peq || !Pv || !Mv) { free(peq); free(Pv); free(Mv); return -1; }
    for (int64_t i = 0; i < m; i++) peq[bcode((unsigned char)p[i]) * W + i / 64] |= (uint64_t)1 << (i % 64);
    for (int64_t w = 0; w < W; w++) Pv[w] = ~(uint64_t)0;
    int lastbit = (int)((m - 1) % 64);
    int64_t score = m;
    for (int64_t j = 0; j < n; j++) {
        const uint64_t *eqr = peq + bcode((unsigned char)t[j]) * W;
        int hin = 1;
        for (int64_t w = 0; w < W; w++) {
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

/* Induced unit cost of an aligned pair of rows (chars, '-' gaps): columns with a difference,
 * both-gap columns skipped. */
int64_t induced_cost(const char *a, const char *b, int64_t L)
{
    int64_t c = 0;
    for (int64_t k = 0; k < L; k++)
        if (a[k] != b[k]) c++;          /* a both-gap column has a[k] == b[k] */
    return c;
}
