#!/bin/bash
# Run after the st_amed alignment: call and score it on chr20, then compare it with the k-mer medoid star.
cd /Users/benedictpaten/PycharmProjects/vntr-eval
until grep -q ALIGN_DONE work/iterate/logs/st_amed.chr20.out; do sleep 30; done
work/iterate/mst_chr20.sh st_amed
until grep -q PAIRS_DONE work/iterate/logs/amed_pairs.out 2>/dev/null; do sleep 30; done
source work/stage4/chr20/env.sh
cd work/iterate/perregion
grep -q "'st_amed': 'st_amed__all'" region_metrics.py || sed -i '' "s/'st_maj_med': 'st_maj_med__all'}/'st_maj_med': 'st_maj_med__all', 'st_amed': 'st_amed__all'}/" region_metrics.py
python3 region_metrics.py > /dev/null && python3 boot.py st_amed:st_medoid,st_amed:all,st_amed:st_chm13
echo AMED_ALL_DONE
