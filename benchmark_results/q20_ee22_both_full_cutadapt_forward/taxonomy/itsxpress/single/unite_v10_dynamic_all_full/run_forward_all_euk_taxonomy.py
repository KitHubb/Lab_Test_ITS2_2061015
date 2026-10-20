"""Classify existing forward-only ASVs; preserve fungi-only results."""
import datetime
import json
import os
from pathlib import Path
import shutil
import subprocess

ROOT = Path('/data/home2/ksy/261016_Test_ITS2')
OUT = ROOT / 'benchmark_results/q20_ee22_both_full_cutadapt_forward'
CLASSIFIER = Path('/data/Reference/QIIME2 2027.5/Fungi/unite_ver10_dynamic_all_19.02.2025-Q2-2024.10.qza')
IMAGE = '/data/software/singularity/qiime2_amplicon_2025.7.sif'

def main():
    for method in ('cutadapt', 'itsxpress'):
        target = OUT / 'taxonomy' / method / 'single/unite_v10_dynamic_all_full'
        if target.exists():
            raise FileExistsError(f'Refusing to overwrite {target}')
        assert (OUT / 'dada2' / method / 'single/repseq.qza').is_file()
    assert CLASSIFIER.is_file()
    for method in ('cutadapt', 'itsxpress'):
        target = OUT / 'taxonomy' / method / 'single/unite_v10_dynamic_all_full'
        target.mkdir(parents=True)
        logs = target / 'logs'
        logs.mkdir()
        runtime = OUT / 'runtime/all_euk_taxonomy' / method
        runtime.mkdir(parents=True, exist_ok=True)
        env = os.environ.copy()
        for key in ('TMPDIR', 'TMP', 'TEMP', 'NUMBA_CACHE_DIR', 'MPLCONFIGDIR', 'XDG_CACHE_HOME', 'XDG_CONFIG_HOME'):
            path = runtime / key
            path.mkdir(exist_ok=True)
            env['SINGULARITYENV_' + key] = str(path)
        for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
            env['SINGULARITYENV_' + key] = '4'
        reads = OUT / 'dada2' / method / 'single/repseq.qza'
        (target / 'settings.json').write_text(json.dumps({'reads': str(reads), 'classifier': str(CLASSIFIER), 'confidence': 0.7, 'jobs': 4, 'reads_per_batch': 'auto', 'read_orientation': 'auto'}, indent=2))
        shutil.copy2(__file__, target / Path(__file__).name)
        def run(name, args):
            cmd = ['singularity', 'exec', '--bind', '/data:/data', IMAGE, 'qiime'] + list(map(str, args))
            print(method, name, flush=True)
            with (logs / 'commands.jsonl').open('a') as handle:
                handle.write(json.dumps(cmd) + '\n')
            with (logs / (name + '.log')).open('w') as handle:
                subprocess.run(cmd, cwd=runtime, env=env, stdout=handle, stderr=subprocess.STDOUT, check=True)
        try:
            run('classify', ['feature-classifier', 'classify-sklearn', '--i-reads', reads, '--i-classifier', CLASSIFIER, '--p-confidence', '0.7', '--p-n-jobs', '4', '--p-reads-per-batch', 'auto', '--p-read-orientation', 'auto', '--o-classification', target / 'taxonomy.qza'])
            run('tabulate', ['metadata', 'tabulate', '--m-input-file', target / 'taxonomy.qza', '--o-visualization', target / 'taxonomy.qzv'])
            run('export', ['tools', 'export', '--input-path', target / 'taxonomy.qza', '--output-path', target / 'taxonomy_export'])
            shutil.copy2(target / 'taxonomy_export/taxonomy.tsv', target / 'taxonomy.tsv')
            run('validate', ['tools', 'validate', target / 'taxonomy.qza'])
            (target / 'COMPLETE.json').write_text(json.dumps({'completed': str(datetime.datetime.now())}))
        except Exception as exc:
            (target / 'FAILED.txt').write_text(repr(exc))
            raise
    print('BOTH COMPLETE', flush=True)

if __name__ == '__main__':
    main()
