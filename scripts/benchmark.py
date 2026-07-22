#!/usr/bin/env python3
"""Exact-match benchmarking for pre-normalized, split VCFs within a BED.

This intentionally has no third-party Python dependency. It is not a replacement
for hap.py/vcfeval in complex representations; bcftools norm is mandatory first.
"""
import gzip, sys

calls_path, truth_path, bed_path, overall_out, strata_out = sys.argv[1:]

def op(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)

regions = {}
with op(bed_path) as fh:
    for line in fh:
        if line.startswith(('#','track','browser')) or not line.strip(): continue
        c,s,e,*_ = line.rstrip().split('\t'); regions.setdefault(c, []).append((int(s),int(e)))
for c in regions: regions[c].sort()

def inside(c, pos1):
    p = pos1 - 1
    return any(s <= p < e for s,e in regions.get(c,()))

def read_vcf(path):
    variants = {}
    sample_i = None
    with op(path) as fh:
        for line in fh:
            if line.startswith('#CHROM'):
                cols=line.rstrip().split('\t'); sample_i=9 if len(cols)>9 else None
            if line.startswith('#'): continue
            f=line.rstrip().split('\t')
            if len(f)<8 or not inside(f[0],int(f[1])): continue
            filt=f[6]
            if filt not in ('PASS','.') and 'truth' not in path.lower(): continue
            info=dict(x.split('=',1) if '=' in x else (x,True) for x in f[7].split(';'))
            depth=vaf=None
            if sample_i is not None and len(f)>sample_i:
                fmt=f[8].split(':'); val=f[sample_i].split(':'); d=dict(zip(fmt,val))
                try: depth=int(d.get('DP','.') or '.')
                except ValueError: pass
                try:
                    ad=[int(x) for x in d.get('AD','').split(',')]
                    if len(ad)>1 and sum(ad)>0: vaf=ad[1]/sum(ad)
                except ValueError: pass
                try:
                    if vaf is None: vaf=float(d.get('AF','').split(',')[0])
                except ValueError: pass
            key=(f[0],int(f[1]),f[3],f[4])
            variants[key]=(depth,vaf)
    return variants

calls, truth = read_vcf(calls_path), read_vcf(truth_path)

def vartype(k): return 'SNV' if len(k[2]) == len(k[3]) == 1 else 'INDEL'
def metrics(cset,tset):
    tp=len(cset&tset); fp=len(cset-tset); fn=len(tset-cset)
    precision=tp/(tp+fp) if tp+fp else 0.0
    recall=tp/(tp+fn) if tp+fn else 0.0
    f1=2*precision*recall/(precision+recall) if precision+recall else 0.0
    return tp,fp,fn,precision,recall,f1

with open(overall_out,'w') as out:
    out.write('variant_type\ttp\tfp\tfn\tprecision\trecall\tf1\n')
    for typ in ('SNV','INDEL'):
        c={k for k in calls if vartype(k)==typ}; t={k for k in truth if vartype(k)==typ}
        out.write(typ+'\t'+'\t'.join(map(lambda x:f'{x:.6f}' if isinstance(x,float) else str(x),metrics(c,t)))+'\n')

vaf_bands=[('0-0.05',0,.05),('0.05-0.10',.05,.10),('0.10-0.20',.10,.20),('0.20-1.00',.20,1.000001)]
depth_bands=[('0-19',0,20),('20-49',20,50),('50-99',50,100),('100+',100,10**12)]
with open(strata_out,'w') as out:
    out.write('stratifier\tband\tvariant_type\ttp\tfp\tfn\tprecision\trecall\tf1\n')
    for label,bands,idx in [('VAF',vaf_bands,1),('DEPTH',depth_bands,0)]:
        for band,lo,hi in bands:
            for typ in ('SNV','INDEL'):
                c={k for k,v in calls.items() if vartype(k)==typ and v[idx] is not None and lo<=v[idx]<hi}
                # Query-based stratification: FN truth sites lack query DP/AF and are
                # counted only if caller emitted the locus. Missing values remain explicit.
                t={k for k,v in truth.items() if vartype(k)==typ and v[idx] is not None and lo<=v[idx]<hi}
                vals=metrics(c,t)
                out.write(f'{label}\t{band}\t{typ}\t'+'\t'.join(f'{x:.6f}' if isinstance(x,float) else str(x) for x in vals)+'\n')

