#!/usr/bin/env python3
"""Create a research-only shortlist from a VEP-annotated Mutect2 VCF."""
import gzip, re, sys

vcf, genes_path, out_path = sys.argv[1:]
op = lambda p: gzip.open(p,'rt') if p.endswith('.gz') else open(p)
genes=set()
with open(genes_path) as fh:
    for line in fh:
        if line.startswith(('#','symbol')) or not line.strip(): continue
        genes.add(line.split('\t',1)[0].strip())

protein_terms={'missense_variant','stop_gained','stop_lost','start_lost','frameshift_variant',
               'inframe_insertion','inframe_deletion','splice_acceptor_variant','splice_donor_variant'}
csq_fields=[]; records=[]
with op(vcf) as fh:
    for line in fh:
        if line.startswith('##INFO=<ID=CSQ'):
            m=re.search(r'Format: ([^">]+)',line)
            if m: csq_fields=m.group(1).strip().split('|')
        if line.startswith('#'): continue
        f=line.rstrip().split('\t'); info={}
        for x in f[7].split(';'):
            k,_,v=x.partition('='); info[k]=v
        fmt=f[8].split(':') if len(f)>9 else []
        sample=dict(zip(fmt,f[9].split(':'))) if fmt else {}
        dp=sample.get('DP',''); af=sample.get('AF','').split(',')[0]
        for raw in info.get('CSQ','').split(','):
            a=dict(zip(csq_fields,raw.split('|')))
            consequence=a.get('Consequence',''); symbol=a.get('SYMBOL','')
            try: rare=(not a.get('gnomADe_AF')) or float(a['gnomADe_AF']) < .01
            except ValueError: rare=True
            clin=a.get('CLIN_SIG','').lower()
            damaging=('deleterious' in a.get('SIFT','').lower() or 'damaging' in a.get('PolyPhen','').lower())
            flags=[]
            if any(t in protein_terms for t in consequence.split('&')): flags.append('protein_altering')
            if rare: flags.append('rare')
            if symbol in genes: flags.append('cancer_gene')
            if 'pathogenic' in clin: flags.append('clinvar_pathogenic')
            if damaging: flags.append('predicted_damaging')
            relevant={'cancer_gene','clinvar_pathogenic','predicted_damaging'}
            if {'protein_altering','rare'} <= set(flags) and set(flags) & relevant:
                records.append([f[0],f[1],f[3],f[4],symbol,a.get('Feature',''),consequence,
                    a.get('HGVSc',''),a.get('HGVSp',''),a.get('CANONICAL',''),a.get('Existing_variation',''),
                    a.get('gnomADe_AF',''),a.get('CLIN_SIG',''),a.get('SIFT',''),a.get('PolyPhen',''),dp,af,','.join(flags)])
cols=['chrom','pos','ref','alt','gene','transcript','consequence','hgvsc','hgvsp','canonical','existing_id',
      'gnomade_af','clinvar','sift','polyphen','tumor_dp','tumor_af','prioritization_flags']
with open(out_path,'w') as out:
    out.write('\t'.join(cols)+'\n')
    for row in records: out.write('\t'.join(row)+'\n')
