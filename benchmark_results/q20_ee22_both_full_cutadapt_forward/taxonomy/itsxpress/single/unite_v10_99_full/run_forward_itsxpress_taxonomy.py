"""Classify true R1-only ITSxpress ASVs using the existing UNITE settings."""
import datetime
import json
import shutil
from run_forward_itsxpress_single import O, D, L, Q, run

TAX = O / 'taxonomy/itsxpress/single/unite_v10_99_full'
CLASS = '/data/Reference/QIIME2-2024.10/Fungi/UNITE/10.0/unite_ver10_99_19.02.2025-Q2-2024.10.qza'

def main():
    assert (L / 'COMPLETE.json').exists(), 'DADA2 is not complete'
    assert not (TAX / 'taxonomy.qza').exists(), 'Refusing to overwrite existing taxonomy'
    TAX.mkdir(parents=True, exist_ok=True)
    shutil.copy2(__file__, TAX / 'run_forward_itsxpress_taxonomy.py')
    run(Q, ['feature-classifier', 'classify-sklearn', '--i-reads', D / 'repseq.qza', '--i-classifier', CLASS, '--p-n-jobs', '4', '--p-reads-per-batch', 'auto', '--p-confidence', '0.7', '--o-classification', TAX / 'taxonomy.qza'], 'taxonomy')
    run(Q, ['metadata', 'tabulate', '--m-input-file', TAX / 'taxonomy.qza', '--o-visualization', TAX / 'taxonomy.qzv'], 'taxonomy_summary')
    run(Q, ['tools', 'export', '--input-path', TAX / 'taxonomy.qza', '--output-path', TAX / 'taxonomy_export'], 'taxonomy_export')
    shutil.copy2(TAX / 'taxonomy_export/taxonomy.tsv', TAX / 'taxonomy.tsv')
    run(Q, ['tools', 'validate', TAX / 'taxonomy.qza'], 'validate_taxonomy')
    (TAX / 'COMPLETE.json').write_text(json.dumps({'completed': str(datetime.datetime.now()), 'classifier': CLASS, 'confidence': 0.7}))
    print('TAXONOMY COMPLETE', flush=True)

if __name__ == '__main__':
    main()
