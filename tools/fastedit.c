/* fastedit.c -- exact global unit-cost Levenshtein distance, Myers (1999) / Hyyro bit-vector
 * algorithm with 64-bit blocks. Characters are compared as raw bytes (so 'N' matches only 'N'),
 * exactly as region.edit_distance does in Python.
 *
 * region.py compiles this on first use into config.CACHE_DIR (cc -O3 -shared) and loads it with
 * ctypes; without a C compiler the Python implementation is used. The results are identical.
 */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

long edit_distance(const unsigned char *p, long m, const unsigned char *t, long n) {
    if (m == 0) return n;
    if (n == 0) return m;
    long W = (m + 63) / 64;
    int slot[256];
    int nsym = 0;
    for (int c = 0; c < 256; c++) slot[c] = -1;
    for (long i = 0; i < m; i++) if (slot[p[i]] < 0) slot[p[i]] = nsym++;
    uint64_t *peq = calloc((size_t)(nsym + 1) * W, sizeof(uint64_t)); /* last row: all zero */
    uint64_t *Pv = malloc(W * sizeof(uint64_t)), *Mv = calloc(W, sizeof(uint64_t));
    if (!peq || !Pv || !Mv) { free(peq); free(Pv); free(Mv); return -1; }
    for (long i = 0; i < m; i++) peq[(size_t)slot[p[i]] * W + i / 64] |= (uint64_t)1 << (i % 64);
    for (long w = 0; w < W; w++) Pv[w] = ~(uint64_t)0;
    int lastbit = (int)((m - 1) % 64);
    long score = m;
    for (long j = 0; j < n; j++) {
        int s = slot[t[j]];
        const uint64_t *eqr = peq + (size_t)(s < 0 ? nsym : s) * W;
        int hin = 1;                    /* global alignment: D[0][j] = j */
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
