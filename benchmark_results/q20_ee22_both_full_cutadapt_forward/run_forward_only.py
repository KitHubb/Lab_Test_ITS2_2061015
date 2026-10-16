"""Reproducible R1-only ITS3 required pipeline; no R2 filtering."""
import csv, json, os, pathlib, subprocess, datetime, shutil
R=pathlib.Path('/data/home2/ksy/261016_Test_ITS2')
O=R/'benchmark_results/q20_ee22_both_full_cutadapt_forward'
QI='/data/software/singularity/qiime2_amplicon_2025.7.sif'
CU='/data/software/singularity/read_cleanup_cutadapt-5.2.sif'
FQ='/data/software/singularity/quay.io-biocontainers-fastqc-0.12.1--hdfd78af_0.img'
MQ='/data/software/singularity/quay.io-biocontainers-multiqc-1.27--pyhdfd78af_0.img'
CLASS='/data/Reference/QIIME2-2024.10/Fungi/UNITE/10.0/unite_ver10_99_19.02.2025-Q2-2024.10.qza'
for d in ['logs','input_validation','read_cleanup','cleanup_qc','post_cleanup_qc/fastqc','post_cleanup_qc/multiqc','qiime_manifest','qiime_import','qiime_preprocess/primer/single','post_primer_qc/fastqc','post_primer_qc/multiqc','dada2/cutadapt/single','taxonomy/cutadapt/single/unite_v10_99_full','runtime']:
    (O/d).mkdir(parents=True,exist_ok=True)
env=os.environ.copy()
for key in ['TMPDIR','TMP','TEMP','NUMBA_CACHE_DIR','MPLCONFIGDIR','XDG_CACHE_HOME']:
    p=O/'runtime'/key; p.mkdir(exist_ok=True); env['SINGULARITYENV_'+key]=str(p)
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
    env['SINGULARITYENV_'+key]='4'
def run(img,args,log):
    cmd=['singularity','exec','--bind','/data:/data',img]+list(map(str,args))
    with open(O/'logs/commands.jsonl','a') as f: f.write(json.dumps(cmd)+'\n')
    with open(O/'logs'/log,'w') as f: subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env,check=True)
def qi(args,name): run(QI,['qiime']+args,name+'.log')
def count(p):
    with open(p) as f: return json.load(f)['read_counts']['output']
def main():
    rows=list(csv.DictReader(open(R/'samplesheet.validated.csv')))
    assert len({r['sample_id'] for r in rows})==len(rows)
    for r in rows: assert pathlib.Path(r['fastq_1']).is_file()
    for p in [QI,CU,FQ,MQ,CLASS]: assert pathlib.Path(p).is_file(),p
    shutil.copy2(__file__,O/'run_forward_only.py')
    shutil.copy2(R/'samplesheet.validated.csv',O/'input_validation/source_samplesheet.csv')
    settings={'cleanup':'R1 only; Nextera core + G{10} + A{10}; -q 20,20; minimum length 50', 'primer':'^GCATCGATGAAGAACGCAGC; error-rate 0.1; overlap 19 (full 19nt ITS3); discard-untrimmed', 'readthrough':'GCATATCAATAAGCGGAGGA; error-rate 0.1; overlap 3; no discard; final minimum length 50','dada2':{'max_ee':2,'trunc_q':2,'trim_left':0,'trunc_len':0,'pooling':'independent','chimera':'consensus','n_reads_learn':1000000},'classifier':CLASS,'confidence':0.7}
    (O/'settings.json').write_text(json.dumps(settings,indent=2))
    summaries=[]; retained=[]; cleaned=[]; final=[]
    pre=O/'qiime_preprocess/primer/single'
    for i,r in enumerate(rows):
        sid=r['sample_id']; clean=O/'read_cleanup'/(sid+'.R1.clean.fastq.gz'); p5=pre/(sid+'.ITS3.fastq.gz'); end=pre/(sid+'.R1.primer.fastq.gz')
        j0=O/'read_cleanup'/(sid+'.cutadapt.json'); j1=pre/(sid+'.ITS3.json'); j2=pre/(sid+'.readthrough.json')
        run(CU,['cutadapt','-j','4','-a','CTGTCTCTTATACACATCT','-a','G{10}','-a','A{10}','-q','20,20','-m','50','--json',j0,'-o',clean,r['fastq_1']],sid+'.cleanup.log')
        run(CU,['cutadapt','-j','4','-g','^GCATCGATGAAGAACGCAGC','-e','0.1','-O','19','--discard-untrimmed','--json',j1,'-o',p5,clean],sid+'.ITS3.log')
        run(CU,['cutadapt','-j','4','-a','GCATATCAATAAGCGGAGGA','-e','0.1','-O','3','-m','50','--json',j2,'-o',end,p5],sid+'.readthrough.log')
        with open(j0) as f: raw=json.load(f)['read_counts']['input']
        n=count(j2); summaries.append([sid,raw,count(j0),count(j1),n])
        if count(j0): cleaned.append(clean)
        if n: retained.append([sid,str(end)]); final.append(end)
        print('TRIM {}/{} {} final={}'.format(i+1,len(rows),sid,n),flush=True)
    with open(O/'cleanup_qc/read_retention.tsv','w') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(['sample-id','raw_R1','cleanup','ITS3_detected','final_min50']); w.writerows(summaries)
    with open(O/'input_validation/excluded_zero_reads.tsv','w') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(['sample-id']); w.writerows([[r[0]] for r in summaries if r[-1]==0])
    assert retained,'No samples survived'
    manifest=O/'qiime_manifest/qiime_manifest_r1.tsv'
    with open(manifest,'w') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(['sample-id','absolute-filepath']); w.writerows(retained)
    for stage,files in [('post_cleanup_qc',cleaned),('post_primer_qc',final)]:
        run(FQ,['fastqc','-t','4','-o',O/stage/'fastqc']+files,stage+'.fastqc.log')
        # Compatibility shim matches ITSdetector's MultiQC module.
        run(MQ,['bash','-c','mkdir -p compat_bin; printf "#!/bin/sh\\nexit 0\\n" > compat_bin/ps; chmod +x compat_bin/ps; export PATH="$PWD/compat_bin:$PATH"; multiqc --force --outdir "$1" "$2"','bash',O/stage/'multiqc',O/stage/'fastqc'],stage+'.multiqc.log')
    demux=O/'qiime_import/demux_r1.qza'
    qi(['tools','import','--type','SampleData[SequencesWithQuality]','--input-format','SingleEndFastqManifestPhred33V2','--input-path',manifest,'--output-path',demux],'import')
    qi(['demux','summarize','--i-data',demux,'--o-visualization',O/'qiime_import/demux_r1.qzv'],'demux_summary')
    d=O/'dada2/cutadapt/single'
    qi(['dada2','denoise-single','--i-demultiplexed-seqs',demux,'--p-trim-left','0','--p-trunc-len','0','--p-trunc-q','2','--p-max-ee','2','--p-chimera-method','consensus','--p-pooling-method','independent','--p-n-reads-learn','1000000','--p-n-threads','4','--o-table',d/'table.qza','--o-representative-sequences',d/'repseq.qza','--o-denoising-stats',d/'denoising_stats.qza'],'dada2')
    qi(['metadata','tabulate','--m-input-file',d/'denoising_stats.qza','--o-visualization',d/'denoising_stats.qzv'],'dada2_stats')
    t=O/'taxonomy/cutadapt/single/unite_v10_99_full'
    qi(['feature-classifier','classify-sklearn','--i-reads',d/'repseq.qza','--i-classifier',CLASS,'--p-n-jobs','4','--p-reads-per-batch','auto','--p-confidence','0.7','--o-classification',t/'taxonomy.qza'],'taxonomy')
    qi(['metadata','tabulate','--m-input-file',t/'taxonomy.qza','--o-visualization',t/'taxonomy.qzv'],'taxonomy_summary')
    qi(['tools','export','--input-path',t/'taxonomy.qza','--output-path',t/'taxonomy_export'],'taxonomy_export')
    shutil.copy2(t/'taxonomy_export/taxonomy.tsv',t/'taxonomy.tsv')
    for p in [d/'table.qza',d/'repseq.qza',t/'taxonomy.qza']:
        qi(['tools','validate',p],p.parent.name+'_'+p.stem+'_validate')
    (O/'COMPLETE.json').write_text(json.dumps({'completed':str(datetime.datetime.now()),'samples_retained':len(retained),'samples_input':len(rows)}))
    print('COMPLETE',flush=True)
if __name__=='__main__':
    os.chdir(O/'runtime')
    try: main()
    except Exception as e:
        (O/'FAILED.txt').write_text(repr(e)); raise
