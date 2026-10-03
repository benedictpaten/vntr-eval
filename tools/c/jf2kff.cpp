// jf2kff: write k-mer counts from `jellyfish dump -c` (lines "KMER COUNT") as a KFF file that
// `vg haplotypes -k` reads, without KMC.
//
//   jellyfish count -m 29 -C -s 500M -t 8 -o reads.jf <(gzip -dc reads.fq.gz)
//   jellyfish dump -c -L 2 reads.jf | jf2kff out.kff 29
//   jf2kff --check out.kff        # read it back: k-mers, total count, max count
//
// Encoding A=0 C=1 G=3 T=2 (KMC's), one k-mer per block (max = 1), counts as 2-byte big-endian
// integers capped at MAX_COUNT (65535; KMC's default cap is 255). The k-mers are whatever jellyfish dumped; with -C they are canonical,
// and `vg haplotypes` looks every k-mer up in both orientations, so that is enough.
//
// Build against vg's bundled kff-cpp-api:
//   c++ -O2 -std=c++17 -I $VG/deps/kff-cpp-api/build tools/c/jf2kff.cpp \
//       $VG/lib/libkff.a -o work/bin/jf2kff
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
#include "kff_io.hpp"

static const uint8_t ENC[4] = {0, 1, 3, 2};  // A C G T

static int code(char c) {
    switch (c) {
        case 'A': case 'a': return ENC[0];
        case 'C': case 'c': return ENC[1];
        case 'G': case 'g': return ENC[2];
        case 'T': case 't': return ENC[3];
        default: return -1;
    }
}

// KFF packs 4 bases per byte, most significant first, padding the first byte when k % 4 != 0.
static bool encode(const char* s, size_t k, std::vector<uint8_t>& out) {
    out.assign((k + 3) / 4, 0);
    size_t rem = k & 3, byte = 0, i = 0;
    if (rem) {
        uint8_t v = 0;
        for (; i < rem; i++) {
            int c = code(s[i]);
            if (c < 0) return false;
            v = (v << 2) | c;
        }
        out[byte++] = v;
    }
    for (; i < k; i += 4) {
        uint8_t v = 0;
        for (size_t j = i; j < i + 4; j++) {
            int c = code(s[j]);
            if (c < 0) return false;
            v = (v << 2) | c;
        }
        out[byte++] = v;
    }
    return true;
}

static int check(const char* path) {
    Kff_reader reader(path);
    uint64_t n = 0, total = 0, mx = 0, k = 0, dsz = 0;
    uint8_t* kmer;
    uint8_t* data;
    while (reader.has_next()) {
        if (n == 0) {
            k = reader.get_var("k");
            dsz = reader.get_var("data_size");
        }
        reader.next_kmer(kmer, data);
        uint64_t c = (uint64_t(data[0]) << 8) | data[1];
        n++;
        total += c;
        if (c > mx) mx = c;
    }
    std::cout << "k " << k << "\tdata_size " << dsz << "\tkmers " << n << "\ttotal_count " << total << "\tmax_count " << mx << std::endl;
    return 0;
}

int main(int argc, char** argv) {
    if (argc == 3 && std::strcmp(argv[1], "--check") == 0) {
        return check(argv[2]);
    }
    if (argc != 3 && argc != 4) {
        std::cerr << "usage: jellyfish dump -c ... | jf2kff out.kff K [MAX_COUNT]\n       jf2kff --check out.kff" << std::endl;
        return 1;
    }
    const size_t k = std::stoul(argv[2]);
    // KMC caps counts at 255 unless told otherwise (-cs); pass 255 to match a KMC-made file
    const unsigned long cap = argc == 4 ? std::stoul(argv[3]) : 65535;
    Kff_file file(argv[1], "w");
    file.write_encoding(const_cast<uint8_t*>(ENC));
    Section_GV gv(&file);
    gv.write_var("k", k);
    gv.write_var("max", 1);
    gv.write_var("data_size", 2);
    gv.close();
    Section_Raw raw(&file);
    std::string line;
    std::vector<uint8_t> seq;
    uint8_t data[2];
    uint64_t n = 0, skipped = 0;
    while (std::getline(std::cin, line)) {
        size_t sp = line.find_first_of(" \t");
        if (sp != k) {
            skipped++;
            continue;
        }
        if (!encode(line.c_str(), k, seq)) {
            skipped++;
            continue;
        }
        unsigned long c = std::stoul(line.substr(sp + 1));
        if (c > cap) c = cap;
        data[0] = (c >> 8) & 0xff;
        data[1] = c & 0xff;
        raw.write_compacted_sequence(seq.data(), k, data);
        n++;
    }
    raw.close();
    file.close();
    std::cerr << "jf2kff: wrote " << n << " " << k << "-mers to " << argv[1] << " (" << skipped << " lines skipped)" << std::endl;
    return 0;
}
