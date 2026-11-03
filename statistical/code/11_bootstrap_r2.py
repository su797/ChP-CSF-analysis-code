import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from paths import arguments, font_path
from pathlib import Path
import csv,json,hashlib
import numpy as np
args=arguments(needs_data=True)
W=args.results/'incremental'
W.mkdir(parents=True,exist_ok=True)
src=args.data
rows=list(csv.DictReader(src.open(encoding='utf-8-sig')))
assert rows and len({r['PTID'] for r in rows})==len(rows)
required=['S_AGE','S_PTGENDER','PTEDUCAT','APOE401','S_DX_label','ChPICV','MSD_STREM2CORRECTED','TAU','PTAU','LV_SUM','Hippocampus_SUM','WholeBrain','EstimatedTotalIntraCranialVol']
missing={'','NA','NaN','nan','N/A'}
a=[r for r in rows if all(r[c] not in missing for c in required)]
assert a
v=lambda c:np.array([float(r[c]) for r in a])
cat=lambda c:np.column_stack([np.array([r[c]==level for r in a],dtype=float) for level in sorted({r[c] for r in a})[1:]])
icv=v('EstimatedTotalIntraCranialVol')
x0=np.column_stack([np.ones(len(a)),v('S_AGE'),cat('S_PTGENDER'),v('PTEDUCAT'),cat('APOE401'),np.log(v('LV_SUM')/icv),np.log(v('Hippocampus_SUM')/icv),np.log(v('WholeBrain')/icv)])
x1=np.column_stack([x0,np.log(v('ChPICV'))])
y=np.column_stack([np.log(v(c)) for c in ['MSD_STREM2CORRECTED','TAU','PTAU']])
# Scale non-intercept columns for numerical conditioning; fitted values and R2 are invariant.
for x in [x0,x1]:x[:,1:]=(x[:,1:]-x[:,1:].mean(axis=0))/x[:,1:].std(axis=0)
def rsq(x,y):
 beta,_,rank,_=np.linalg.lstsq(x,y,rcond=None)
 assert rank==x.shape[1]
 return 1-((y-x@beta)**2).sum(axis=0)/((y-y.mean(axis=0))**2).sum(axis=0)
point=np.stack([rsq(x0,y),rsq(x1,y)],axis=1)
saved=[r for r in csv.DictReader((W/'incremental_results.csv').open(encoding='utf-8-sig')) if r['specification']=='primary']
expected=np.array([[float(r['base_r2']),float(r['extended_r2'])] for r in saved])
assert np.allclose(point,expected,atol=1e-12,rtol=0),(point,expected)
seed=20261102;B=10000;rng=np.random.default_rng(seed)
samples=np.empty((B,3,2))
for k in range(B):
 idx=rng.integers(0,len(a),len(a));yb=y[idx]
 samples[k,:,0]=rsq(x0[idx],yb);samples[k,:,1]=rsq(x1[idx],yb)
assert np.all(samples[:,:,1]>=samples[:,:,0]-1e-12)
lo,hi=np.quantile(samples,[.025,.975],axis=0)
out=[]
for j,r in enumerate(saved):
 for m,model in enumerate(['base','extended']):out.append({'biomarker':r['biomarker'],'model':model,'r2':float(point[j,m]),'lower':float(lo[j,m]),'upper':float(hi[j,m])})
with (W/'figure4_r2_bootstrap_ci.csv').open('w',newline='',encoding='utf8') as f:
 wr=csv.DictWriter(f,fieldnames=out[0].keys());wr.writeheader();wr.writerows(out)
np.save(W/'bootstrap_r2_replicates.npy',samples)
meta={'n':len(a),'replicates':B,'seed':seed,'method':'Participant-level nonparametric pairs bootstrap; same resampled participant indices for both models and all outcomes; 95% percentile intervals (2.5th and 97.5th percentiles).','interpretation':'Intervals describe apparent within-sample R2, not validated prediction performance.','maximum_point_estimate_error':float(abs(point-expected).max()),'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'intervals':out}
(W/'bootstrap_validation.json').write_text(json.dumps(meta,indent=2),encoding='utf8')
print(f'R2 bootstrap: {B} resamples, n = {len(a)}')
