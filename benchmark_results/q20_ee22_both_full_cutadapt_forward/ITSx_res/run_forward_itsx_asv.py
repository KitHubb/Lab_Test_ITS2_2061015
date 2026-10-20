"""ITSx on existing R1-only ASVs, followed by two UNITE classifications."""
import datetime
import json
import os
from pathlib import Path
import shutil
import subprocess

ROOT = Path('/data/home2/ksy/261016_Test_ITS2')
OUT = ROOT / 'benchmark_results/q20_ee22_both_full/ITSx_res'
SOURCE = ROOT / 'benchmark_results/q20_ee22_both_full_cutadapt_forward/dada2/cutadapt/single/repseq.qza'
QIIME = '/data/software/singularity/qiime2_amplicon_2025.7.sif'
ITSX = '/data/software/singularity/itsx_1.1.3.sif'
CLASSIFIERS = {
    'unite_v10_99_full': '/data/Reference/QIIME2-2024.10/Fungi/UNITE/10.0/unite_ver10_99_19.02.2025-Q2-2024.10.qza',
    'unite_v10_dynamic_all_full': '/data/Reference/QIIME2 2027.5/Fungi/unite_ver10_dynamic_all_19.02.2025-Q2-2024.10.qza',
}

def main():
    OUT.mkdir(parents=True, exist_ok=False)
    logs = OUT / 'logs'
    logs.mkdir()
    runtime = OUT / 'runtime'
    runtime.mkdir()
    env = os.environ.copy()
    for key in ('TMPDIR', 'TMP', 'TEMP', 'NUMBA_CACHE_DIR', 'MPLCONFIGDIR', 'XDG_CACHE_HOME', 'XDG_CONFIG_HOME'):
        path = runtime / key
        path.mkdir()
        env['SINGULARITYENV_' + key] = str(path)
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        env['SINGULARITYENV_' + key] = '4'
    shutil.copy2(__file__, OUT / Path(__file__).name)
    (OUT / 'settings.json').write_text(json.dumps({'input': str(SOURCE), 'itsx': {'version': '1.1.3', 'taxa': 'F', 'save_regions': 'ITS2', 'partial_cutoff': 50, 'preserve': True, 'cpu': 4, 'other_options': 'ITSx defaults'}, 'classifiers': CLASSIFIERS, 'confidence': 0.7, 'note': 'ASV IDs preserved; trimmed sequences are not re-denoised or collapsed'}, indent=2))
    def run(name, image, args):
        cmd = ['singularity', 'exec', '--bind', '/data:/data', image] + list(map(str, args))
        print(name, flush=True)
        with (logs / 'commands.jsonl').open('a') as handle:
            handle.write(json.dumps(cmd) + '\n')
        with (logs / (name + '.log')).open('w') as handle:
            subprocess.run(cmd, env=env, cwd=runtime, stdout=handle, stderr=subprocess.STDOUT, check=True)
    def qi(name, args):
        run(name, QIIME, ['qiime'] + args)
    try:
        qi('export_input', ['tools', 'export', '--input-path', SOURCE, '--output-path', OUT / 'input'])
        summary = {}
        for branch, extra, suffix in [('default', [], '.ITS2.fasta'), ('partial50', ['--partial', '50'], '.ITS2.full_and_partial.fasta')]:
            dest = OUT / branch
            dest.mkdir()
            prefix = dest / 'itsx'
            run(branch + '_itsx', ITSX, ['ITSx', '-i', OUT / 'input/dna-sequences.fasta', '-o', prefix, '-t', 'F', '--save_regions', 'ITS2', '--preserve', 'T', '--cpu', '4'] + extra)
            fasta = Path(str(prefix) + suffix)
            # ITSx may place outputs in a directory named after its prefix.
            if not fasta.exists():
                candidates = list(dest.rglob('itsx' + suffix))
                if len(candidates) != 1:
                    raise RuntimeError(f'Cannot resolve output FASTA: {fasta}, candidates={candidates}')
                fasta = candidates[0]
            headers = [line[1:].split()[0] for line in fasta.read_text().splitlines() if line.startswith('>')]
            if not headers or len(headers) != len(set(headers)):
                raise RuntimeError(f'Empty or duplicate-ID FASTA: {fasta}')
            source_ids = {line[1:].split()[0] for line in (OUT / 'input/dna-sequences.fasta').read_text().splitlines() if line.startswith('>')}
            assert set(headers) <= source_ids, 'ITSx changed ASV IDs'
            summary[branch] = {'asvs': len(headers), 'fasta': str(fasta)}
            (OUT / 'summary.json').write_text(json.dumps(summary, indent=2))
            qi(branch + '_import', ['tools', 'import', '--type', 'FeatureData[Sequence]', '--input-path', fasta, '--output-path', dest / 'repseq.qza'])
            qi(branch + '_validate_repseq', ['tools', 'validate', dest / 'repseq.qza'])
            for label, classifier in CLASSIFIERS.items():
                tax = dest / 'taxonomy' / label
                tax.mkdir(parents=True)
                tag = branch + '_' + label
                qi(tag + '_classify', ['feature-classifier', 'classify-sklearn', '--i-reads', dest / 'repseq.qza', '--i-classifier', classifier, '--p-confidence', '0.7', '--p-n-jobs', '4', '--p-reads-per-batch', 'auto', '--o-classification', tax / 'taxonomy.qza'])
                qi(tag + '_tabulate', ['metadata', 'tabulate', '--m-input-file', tax / 'taxonomy.qza', '--o-visualization', tax / 'taxonomy.qzv'])
                qi(tag + '_export', ['tools', 'export', '--input-path', tax / 'taxonomy.qza', '--output-path', tax / 'taxonomy_export'])
                shutil.copy2(tax / 'taxonomy_export/taxonomy.tsv', tax / 'taxonomy.tsv')
                qi(tag + '_validate', ['tools', 'validate', tax / 'taxonomy.qza'])
        (OUT / 'COMPLETE.json').write_text(json.dumps({'completed': str(datetime.datetime.now()), 'summary': summary}, indent=2))
        print('COMPLETE', flush=True)
    except Exception as exc:
        (OUT / 'FAILED.txt').write_text(repr(exc))
        raise

if __name__ == '__main__':
    main()
