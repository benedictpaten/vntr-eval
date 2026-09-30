#include <stdint.h>
#include <stdlib.h>
#include <string.h>
static long semiglobal(const unsigned char *p, long m, const unsigned char *t, long n) {
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
    long score = m, best = m;
    for (long j = 0; j < n; j++) {
        int s = slot[t[j]];
        const uint64_t *eqr = peq + (size_t)(s < 0 ? nsym : s) * W;
        int hin = 0;                    /* semi-global: D[0][j] = 0 (free start in text) */
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
        if (score < best) best = score;
    }
    free(peq); free(Pv); free(Mv);
    return best;
}

/* semig haps.fa reads.fa -> one line per read: name, then the semi-global edit distance of the read
 * (and of its reverse complement; the minimum is reported) to each haplotype, tab-separated. */
#include <stdio.h>
static char *slurp_fa(const char *path, char ***names, char ***seqs, long **lens, long *n) {
    FILE *f = fopen(path, "r"); if (!f) { perror(path); exit(1); }
    long cap = 1024; *names = malloc(cap * sizeof(char*)); *seqs = malloc(cap * sizeof(char*)); *lens = malloc(cap * sizeof(long)); *n = 0;
    char *line = NULL; size_t lc = 0; ssize_t r; char *cur = NULL; long cl = 0, cc = 0;
    while ((r = getline(&line, &lc, f)) > 0) {
        while (r > 0 && (line[r-1] == '\n' || line[r-1] == '\r')) line[--r] = 0;
        if (line[0] == '>') {
            if (*n == cap) { cap *= 2; *names = realloc(*names, cap*sizeof(char*)); *seqs = realloc(*seqs, cap*sizeof(char*)); *lens = realloc(*lens, cap*sizeof(long)); }
            if (*n > 0) { (*seqs)[*n-1] = cur; (*lens)[*n-1] = cl; }
            char *sp = strchr(line + 1, ' '); if (sp) *sp = 0;
            (*names)[(*n)++] = strdup(line + 1); cur = malloc(1024); cl = 0; cc = 1024;
        } else {
            if (cl + r + 1 > cc) { cc = (cl + r + 1) * 2; cur = realloc(cur, cc); }
            for (ssize_t i = 0; i < r; i++) { char c = line[i]; if (c >= 'a' && c <= 'z') c -= 32; cur[cl++] = c; }
        }
    }
    if (*n > 0) { (*seqs)[*n-1] = cur; (*lens)[*n-1] = cl; }
    fclose(f); free(line); return NULL;
}
int main(int argc, char **argv) {
    if (argc != 3) { fprintf(stderr, "usage: semig haps.fa reads.fa\n"); return 2; }
    char **hn, **hs, **rn, **rs; long *hl, *rl, nh, nr;
    slurp_fa(argv[1], &hn, &hs, &hl, &nh); slurp_fa(argv[2], &rn, &rs, &rl, &nr);
    char *rc = NULL; long rcc = 0;
    for (long i = 0; i < nr; i++) {
        if (rl[i] + 1 > rcc) { rcc = rl[i] + 1; rc = realloc(rc, rcc); }
        for (long k = 0; k < rl[i]; k++) { char c = rs[i][rl[i]-1-k];
            rc[k] = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : c == 'T' ? 'A' : 'N'; }
        printf("%s", rn[i]);
        for (long h = 0; h < nh; h++) {
            long a = semiglobal((unsigned char*)rs[i], rl[i], (unsigned char*)hs[h], hl[h]);
            long b = semiglobal((unsigned char*)rc, rl[i], (unsigned char*)hs[h], hl[h]);
            printf("\t%ld", a < b ? a : b);
        }
        printf("\n");
    }
    return 0;
}
