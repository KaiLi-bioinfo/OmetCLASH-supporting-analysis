#!/usr/bin/env python3
"""Minimal per-library sequence-weighted QC summary, excluding annotated non-piRNA small RNAs."""
import argparse
import pandas as pd
p = argparse.ArgumentParser()
p.add_argument('table'); p.add_argument('output'); p.add_argument('--min-length',type=int,default=26); p.add_argument('--max-length',type=int,default=33)
a=p.parse_args()
d=pd.read_csv(a.table,sep='\t',usecols=['sequence','count','class'])
d=d[d['class'].isin(['genome_0mm','genome_1mm'])].copy()
d['length']=d.sequence.str.len()
d=d[d.length.between(a.min_length,a.max_length)]
total=d['count'].sum()
metrics={'reads':int(total),'unique_sequences':len(d),'fraction_1U':float(d.loc[d.sequence.str.startswith('T'),'count'].sum()/total) if total else 0,'fraction_10A':float(d.loc[d.sequence.str[9]=='A','count'].sum()/total) if total else 0}
for length, group in d.groupby('length'):
    metrics[f'length_{length}_fraction']=float(group['count'].sum()/total) if total else 0
pd.Series(metrics,name='value').rename_axis('metric').to_csv(a.output,sep='\t')
