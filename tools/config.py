"""config.py -- where vntr-eval finds its data and tools.

Every location has a default and an environment variable that overrides it. Nothing here
points at a scratch directory: the large inputs live in the vg-call-eval repository (read only),
and everything vntr-eval writes goes under this repository.

    VNTR_EVAL_DIR    the vg-call-eval repository holding graphs, reads, truth and calls
                     (default ~/PycharmProjects/vg-call-eval)
    VNTR_REGIONS     region packages (default <repo>/regions)
    VNTR_CANDIDATES  candidate graphs (default <repo>/candidates)
    VNTR_RESULTS     results tables (default <repo>/results)
    VNTR_CENSUS      census tables: loci.tsv, strata.tsv, vntr_regions.tsv, catalogue.tsv
                     (default <repo>/census; optional, used only for annotation)
    VNTR_WORK        caches, toolkit output and viewer intermediates, not tracked
                     (default <repo>/work)
    VNTR_GBZ_BASE, VNTR_VG, VNTR_MAFFT, VNTR_BANDAGE, VNTR_BCFTOOLS, VNTR_SAMTOOLS,
    VNTR_BGZIP, VNTR_CC  individual tool binaries

Use from a tool:

    import config
    p = config.data_paths('chr4')      # dict of per-contig data files
    config.have_graph_data()           # False on a machine with only the region packages
"""
import os
import shutil

HOME = os.path.expanduser('~')
TOOLS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(TOOLS_DIR)


def _env(name, default):
    v = os.environ.get(name)
    return os.path.expanduser(v) if v else default


EVAL_DIR = _env('VNTR_EVAL_DIR', os.path.join(HOME, 'PycharmProjects/vg-call-eval'))
REGIONS_DIR = _env('VNTR_REGIONS', os.path.join(REPO, 'regions'))
CANDIDATES_DIR = _env('VNTR_CANDIDATES', os.path.join(REPO, 'candidates'))
RESULTS_DIR = _env('VNTR_RESULTS', os.path.join(REPO, 'results'))
CENSUS_DIR = _env('VNTR_CENSUS', os.path.join(REPO, 'census'))
WORK_DIR = _env('VNTR_WORK', os.path.join(REPO, 'work'))
CACHE_DIR = os.path.join(WORK_DIR, 'cache')

# Search path for tools: the user's local bin and Homebrew first, then $PATH.
TOOL_PATH = os.pathsep.join([os.path.join(HOME, '.local/bin'), '/opt/homebrew/bin',
                             '/usr/local/bin', os.environ.get('PATH', '')])


def _tool(env, name, fallback=None):
    v = os.environ.get(env)
    if v:
        return os.path.expanduser(v)
    return shutil.which(name, path=TOOL_PATH) or fallback or name


GBZ_BASE = _tool('VNTR_GBZ_BASE', 'gbz-base', os.path.join(HOME, '.local/bin/gbz-base'))
_VG_BUILD = os.path.join(HOME, 'CLionProjects/vg/bin/vg')   # the vg build the evaluation used
VG = (os.path.expanduser(os.environ['VNTR_VG']) if os.environ.get('VNTR_VG')
      else _VG_BUILD if os.path.exists(_VG_BUILD) else shutil.which('vg', path=TOOL_PATH) or 'vg')
MAFFT = _tool('VNTR_MAFFT', 'mafft')
BANDAGE = _tool('VNTR_BANDAGE', 'Bandage', os.path.join(HOME, '.local/bin/Bandage'))
BCFTOOLS = _tool('VNTR_BCFTOOLS', 'bcftools')
SAMTOOLS = _tool('VNTR_SAMTOOLS', 'samtools')
BGZIP = _tool('VNTR_BGZIP', 'bgzip')
CC = _tool('VNTR_CC', 'cc')

TRUTH_PREFIX = 'CHM13v2.0_HG2-T2TQ100-V1.1'


def data_paths(contig=None):
    """Per-contig data files in the vg-call-eval repository (they may not exist)."""
    e = EVAL_DIR
    c = contig or 'CONTIG'
    t = os.path.join(e, 'data/truth', TRUTH_PREFIX)
    return {
        'gbz_db': os.path.join(e, 'work/graph.hap32.gbz.db'),
        'gaf_db': os.path.join(e, 'work/reads.hap32.gaf.db'),
        'hprc_gbz': os.path.join(e, 'data/hprc-v2.1-mc-chm13-eval.gref.gbz'),
        'contig_gbz': os.path.join(e, 'work/wgs', c, c + '.gbz'),
        'ref_fa': os.path.join(e, 'work/wgs', c, c + '.fa'),
        'stvar': t + '_stvar.vcf.gz',
        'smvar': t + '_smvar.vcf.gz',
        'stvar_bed': t + '_stvar.benchmark.bed',
        'smvar_bed': t + '_smvar.benchmark.bed',
        'vg_vcf': os.path.join(e, 'work/wgs-mm095', c, c + '.vcf.gz'),
        'vg_truvari': os.path.join(e, 'work/wgs-mm095/score', c + '.truvari'),
        'vg_aardvark': os.path.join(e, 'work/wgs-mm095/score', c + '.aardvark'),
        'pg_norm': os.path.join(e, 'work/pangenie/score', c + '.norm.vcf.gz'),
        'pg_truvari': os.path.join(e, 'work/pangenie/score', c + '.truvari'),
        'pg_aardvark': os.path.join(e, 'work/pangenie/score', c + '.aardvark'),
        'rmsk': os.path.join(e, 'data/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed'),
    }


def have_graph_data():
    """True when the hap32 graph database is readable (the big-data machine)."""
    return os.path.exists(data_paths()['gbz_db']) and os.path.exists(GBZ_BASE)


def census_file(name):
    """Path of a census table (loci.tsv, strata.tsv, vntr_regions.tsv, catalogue.tsv) or None."""
    p = os.path.join(CENSUS_DIR, name)
    return p if os.path.exists(p) else None


def tool_env():
    """Environment for subprocesses, with the tool search path."""
    return dict(os.environ, PATH=TOOL_PATH)


def describe():
    return {k: v for k, v in globals().items()
            if k.isupper() and isinstance(v, str)}


if __name__ == '__main__':
    import json
    d = describe()
    d['have_graph_data'] = have_graph_data()
    print(json.dumps(d, indent=1))
