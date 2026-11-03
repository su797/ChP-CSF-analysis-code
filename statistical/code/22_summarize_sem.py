"""Summarize valid bootstrap fits and paired serial/single-path ratios."""
import sys, csv
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from paths import arguments
args=arguments()
folder=args.results/'sem'
def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f))
def write(name,rows):
    with (folder/name).open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
def number(row,key):
    try: return float(row[key])
    except (ValueError,KeyError): return np.nan
point=read(folder/'point_estimates.csv')
boot={};summary=[]
for file in sorted((folder/'bootstrap').glob('*.csv')):
    rows=read(file)
    if len(rows)!=1000: raise ValueError(f'Incomplete bootstrap: {file.name}')
    valid=np.array([r['valid']=='TRUE' for r in rows])
    ident=file.stem;definition,kind,model,group=ident.split('__')
    pe=[r for r in point if r['id']==ident]
    if not pe: raise ValueError(f'Missing point estimates: {ident}')
    boot[ident]=(rows,valid)
    for column in rows[0]:
        if not column.endswith('__std'): continue
        effect=column[:-5]
        matches=[r for r in pe if (r['op']==':=' and r['lhs']==effect) or (r['op']=='~' and r['label']==effect)]
        if len(matches)!=1: raise ValueError((ident,effect,len(matches)))
        p=matches[0]
        raw=np.array([number(r,effect+'__est') for r in rows])
        std=np.array([number(r,column) for r in rows])
        low,high=np.quantile(raw[valid],[.025,.975]);sl,sh=np.quantile(std[valid],[.025,.975])
        summary.append(dict(id=ident,definition=definition,kind=kind,model=model,group=group,effect=effect,
            point_std=number(p,'std.all'),ci_std_low=sl,ci_std_high=sh,point_unstd=number(p,'est'),
            ci_unstd_low=low,ci_unstd_high=high,MLR_p=number(p,'pvalue'),B=1000,valid=int(valid.sum()),invalid_rate=1-valid.mean()))
write('bootstrap_summary.csv',summary)
ratios=[]
for marker in ['PTAU','TAU']:
    si=f'raw__serial__{marker}_to_sTREM2_to_ChP_to_Cognition__Overall'
    ti=f'raw__single__{marker}_to_ChP_to_Cognition__Overall'
    s,sv=boot[si];t,tv=boot[ti];valid=sv&tv
    ps=next(r for r in point if r['id']==si and r['lhs']=='serial_indirect' and r['op']==':=')
    pt=next(r for r in point if r['id']==ti and r['lhs']=='indirect' and r['op']==':=')
    # Matching iteration numbers refer to the same resampled participants.
    for scale,col in [('est','est'),('std','std.all')]:
        numerator=np.array([number(r,'serial_indirect__'+scale) for r in s])
        denominator=np.array([number(r,'indirect__'+scale) for r in t])
        values=100*numerator[valid]/denominator[valid]
        low,high=np.quantile(values,[.025,.975])
        ratios.append(dict(definition='raw',marker=marker,scale=scale,point_pct=100*number(ps,col)/number(pt,col),
                           ci_low=low,ci_high=high,valid=int(valid.sum())))
write('serial_crossmodel_ratio_diagnostic.csv',ratios)
print(f'{len(boot)} bootstrap models; {len(summary)} parameter summaries')
