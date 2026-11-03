"""Supplementary plots from model estimates and the analysis dataset."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import FancyBboxPatch
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from paths import arguments, ROOT
args=arguments(needs_data=True)
A=args.results/'association'; S=args.results/'sem'
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
    'font.size':10,'axes.labelsize':10,'axes.titlesize':11,'axes.spines.top':False,'axes.spines.right':False,
    'svg.fonttype':'none','pdf.fonttype':42,'savefig.facecolor':'white'})
COLORS=['#1F5A8A','#E8751A','#B91D3A'];GROUPS=['CN','MCI','AD']
NAMES={'MSD_STREM2CORRECTED':'sTREM2','TAU':'Tau','PTAU':'P-Tau','S_ABETA':'CSF Aβ42',
       'ChPICV':'ChP/ICV','ChP_SUM':'Choroid plexus volume','Hippocampus_SUM':'Hippocampal volume',
       'Amygdala_SUM':'Amygdalar volume','LV_SUM':'Lateral ventricular volume','MMSE':'MMSE','MOCA':'MoCA','S_AGE':'Age'}
cmap=LinearSegmentedColormap.from_list('effects',['#4BA3D5','white','#F2A44A'])
def read(name,sem=False): return pd.read_csv((S if sem else A)/name)
def label(s):
    for a,b in [('MSD_STREM2CORRECTED','sTREM2'),('ChPICV','ChP/ICV'),('PTAU','P-Tau'),('TAU','Tau'),('ABETA','CSF Aβ42'),('_to_',' → '),('_on_',' ~ '),('_vs_',' vs '),('AD','Dementia')]: s=s.replace(a,b)
    return s
def stars(p): return '***' if p<.001 else '**' if p<.01 else '*' if p<.05 else ''
def save(fig,n):
    stem='Supplementary_Figure_S'+str(n)
    fig.savefig(args.output/(stem+'.png'),dpi=300,bbox_inches='tight')
    fig.savefig(args.output/(stem+'.svg'),bbox_inches='tight')
    plt.close(fig)
def heat(values,rows,cols,n,limit=.3,symbols=None,q=None,bar='Unstandardized indirect effect (b)'):
    wide=len(cols)>4
    fig,ax=plt.subplots(figsize=(10 if wide else 8.3,max(3.5,len(rows)*.43+1.4)))
    im=ax.imshow(values,cmap=cmap,norm=TwoSlopeNorm(vcenter=0,vmin=-limit,vmax=limit),aspect='auto')
    ax.set_xticks(range(len(cols)),cols,rotation=25 if wide else 0);ax.xaxis.tick_top()
    ax.set_yticks(range(len(rows)),rows);ax.tick_params(length=0,pad=8)
    for i in range(len(rows)):
        for j in range(len(cols)):
            v=values[i,j]
            ax.text(j,i,'' if not np.isfinite(v) else f'{v:.3f}'+(symbols[i][j] if symbols else ''),ha='center',va='center',fontsize=9)
    for x in np.arange(-.5,len(cols),1):ax.axvline(x,color='white',lw=.6)
    for y in np.arange(-.5,len(rows),1):ax.axhline(y,color='white',lw=.6)
    if q is not None:
        ax.set_xlim(-.5,len(cols)+.8)
        ax.text(len(cols)+.12,-.9,'Omnibus $q$',ha='center',fontweight='bold')
        for i,v in enumerate(q):ax.text(len(cols)+.12,i,f'{v:.3f}'.lstrip('0'),ha='center',va='center')
    fig.colorbar(im,ax=ax,fraction=.025,pad=.04,label=bar)
    if symbols:fig.text(.5,.015,'* $P$ < .05, ** $P$ < .01, *** $P$ < .001',ha='center',fontsize=9)
    fig.subplots_adjust(left=.18 if wide else .42,right=.88,bottom=.13,top=.8 if wide else .86)
    save(fig,n)

# Partial correlations and comparator-structure coefficients have different scales.
t=read('association_heatmap_partial_correlations.csv')
vs=['ChPICV','MSD_STREM2CORRECTED','TAU','PTAU','S_ABETA','MMSE','MOCA']
matrix=np.eye(len(vs));sym=[['' for _ in vs] for _ in vs]
for i,x in enumerate(vs):
    for j,y in enumerate(vs):
        if i!=j:
            r=t[(t.row_var==x)&(t.col_var==y)].iloc[0];matrix[i,j]=r.estimate;sym[i][j]=stars(r.p_value)
heat(matrix,[NAMES[v] for v in vs],[NAMES[v] for v in vs],1,1,sym,bar='Partial correlation (r)')
for n,file,samples in [(2,'matching_smd_summary.csv',['Unmatched','Matched']),
                       (3,'overlap_weighting_balance.csv',['Unweighted','Overlap weighted'])]:
    t=read(file);t['pair']=t['pair'].str.replace('AD','Dementia',regex=False)
    keys=t[['pair','variable']].drop_duplicates().apply(tuple,axis=1).tolist()
    fig,ax=plt.subplots(figsize=(8,4.4))
    for sample,marker,color in zip(samples,['o','^'],COLORS):
        tt=t[t['sample']==sample]
        vals=[abs(tt[(tt.pair==p)&(tt.variable==v)].iloc[0].smd) for p,v in keys]
        ax.scatter(vals,range(len(keys)),s=70,marker=marker,color=color,label=sample,zorder=3)
    ax.set_yticks(range(len(keys)),[p+' | '+v for p,v in keys]);ax.invert_yaxis()
    ax.axvline(.1,color='#cccccc',ls='--',lw=.8);ax.set_xlim(-.02,.44)
    ax.set_xlabel('Absolute standardized mean difference')
    ax.legend(loc='upper right',frameon=True,facecolor='white',edgecolor='#888888',framealpha=1)
    fig.tight_layout();save(fig,n)

data=pd.read_csv(args.data)
def distributions(variables,n):
    fig,axes=plt.subplots(2 if len(variables)>2 else 1,3 if len(variables)>2 else 1,figsize=(11,7) if len(variables)>2 else (5.5,4),squeeze=False)
    rng=np.random.default_rng(20260414)
    for ax,var in zip(axes.flat,variables):
        vals=[data.loc[data.S_DX_label==g,var].dropna().to_numpy() for g in GROUPS]
        # Use the same analysis-scale transformation as the regression tables.
        trans=pd.read_csv(ROOT/'config/transformation_plan.csv')
        tr=trans[trans.variable==var].iloc[0].transform_method
        if tr=='log':vals=[np.log(v) for v in vals]
        elif tr=='log1p':vals=[np.log1p(v) for v in vals]
        violins=ax.violinplot(vals,positions=[0,1,2],widths=.7,showextrema=False)
        for i,(body,v,col) in enumerate(zip(violins['bodies'],vals,COLORS)):
            vertices=body.get_paths()[0].vertices;vertices[:,0]=np.minimum(vertices[:,0],i)
            body.set_facecolor(col);body.set_alpha(.35)
            ax.scatter(i+rng.uniform(.08,.32,len(v)),v,s=5,alpha=.3,c=col,edgecolors='none')
        ax.boxplot(vals,positions=[0,1,2],widths=.12,showfliers=False,medianprops={'color':'black'})
        ax.set_xticks([0,1,2],['CN','MCI','Dementia'])
        ax.set_ylabel(('ln(' if tr=='log' else '')+NAMES[var]+(')' if tr=='log' else ''))
        comparisons=read('group_comparisons_pairwise.csv')
        # Pairwise tests are saved separately; the plot displays those adjusted values.
        selected=comparisons[comparisons.variable==var]
        if len(selected):
            pcol=next((c for c in ['p_adjusted','p_adj','p.adj','p_value'] if c in selected),None)
            if pcol:
                caption='\n'.join(str(r['comparison']).replace('AD','Dementia')+f': {r[pcol]:.3g}' for _,r in selected.iterrows())
                ax.set_title('Pairwise $P$\n'+caption,fontsize=8,pad=12)
    for ax in list(axes.flat)[len(variables):]:ax.set_visible(False)
    fig.tight_layout();save(fig,n)
distributions(['ChPICV','MSD_STREM2CORRECTED','TAU','PTAU','S_ABETA'],4)
distributions(['S_AGE'],5)

def forest(ax,t,estimate='std_beta',lo='std_conf.low',hi='std_conf.high',labels=None,colors=None,xlabel='Standardized β'):
    if lo not in t:
        ratio=t[estimate]/t['estimate'];t=t.copy()
        t[lo]=np.minimum(t['conf.low']*ratio,t['conf.high']*ratio);t[hi]=np.maximum(t['conf.low']*ratio,t['conf.high']*ratio)
    for i,(_,r) in enumerate(t.iterrows()):
        color=colors[i%len(colors)] if colors else COLORS[0]
        ax.plot([r[lo],r[hi]],[i,i],color=color,lw=1.4)
        ax.scatter(r[estimate],i,color=color,s=35,zorder=3)
    ax.axvline(0,color='#aaaaaa',lw=.8,ls='--');ax.set_yticks(range(len(t)),labels if labels is not None else range(len(t)))
    ax.invert_yaxis();ax.set_xlabel(xlabel)

t=read('overlap_weighting_compare_unweighted.csv')
fig,axes=plt.subplots(1,2,figsize=(11,5))
for letter,ax,vars in zip(['A','B'],axes,[['ChPICV_on_ABETA','ChPICV_on_PTAU'],['ChPICV_on_sTREM2','ChPICV_on_TAU']]):
    tt=t[t.model_name.isin(vars)].copy();labs=[label(r.model_name)+' | '+r.analysis_type for r in tt.itertuples()]
    types=list(t.analysis_type.unique()); palette=dict(zip(types,COLORS))
    forest(ax,tt,labels=labs,colors=[palette[v] for v in tt.analysis_type])
    ax.set_title(letter,loc='left',fontweight='bold')
fig.tight_layout();save(fig,8)
t=read('matched_key_model_sensitivity.csv');fig,axes=plt.subplots(1,3,figsize=(12,5))
for ax,(pair,tt) in zip(axes,t.groupby('pair_name',sort=False)):
    forest(ax,tt,lo='std_low',hi='std_high',labels=[label(r.model_name)+' | '+r.sample_set.replace('_pair','') for r in tt.itertuples()],colors=[COLORS[0] if v.startswith('unmatched') else COLORS[1] for v in tt.sample_set])
    ax.set_title(label(pair))
fig.tight_layout();save(fig,9)
t=read('abeta_truncation_regression_sensitivity.csv');t=t[(t.term=='exposure')&(t.group=='Overall')]
fig,axes=plt.subplots(1,2,figsize=(10,4))
forest(axes[0],t,labels=[r.model_name+' | '+r.sample_set for r in t.itertuples()])
tt=read('abeta_censored_sensitivity.csv');tt=tt[(tt.term=='chp')&(tt.group=='Overall')].copy()
tt['low']=(tt.estimate-1.96*tt['std.error'])/1e5;tt['high']=(tt.estimate+1.96*tt['std.error'])/1e5;tt['scaled']=tt.estimate/1e5
forest(axes[1],tt,'scaled','low','high',['Right-censored model'],xlabel='Slope / 10⁵ (pg/mL per unit ChP/ICV)')
for letter,ax in zip(['A','B'],axes):ax.set_title(letter,loc='left',fontweight='bold')
fig.tight_layout();save(fig,10)
t=read('alternative_chp_definition_models.csv');fig,ax=plt.subplots(figsize=(8,5.5))
forest(ax,t,labels=[r.model_name+' | '+label(r.group) for r in t.itertuples()],colors=COLORS)
fig.tight_layout();save(fig,11)
t=read('structure_specificity_models.csv');structures=['ChP_SUM','Hippocampus_SUM','Amygdala_SUM','LV_SUM'];markers=['MSD_STREM2CORRECTED','TAU','PTAU','S_ABETA']
matrix=[];sym=[]
for st in structures:
    rs=[t[(t.structure==st)&(t.biomarker==m)&(t.group=='Overall')].iloc[0] for m in markers]
    matrix.append([r.std_beta for r in rs]);sym.append([stars(r['p.value']) for r in rs])
heat(np.array(matrix),[NAMES[v] for v in structures],[NAMES[v] for v in markers],12,.3,sym,bar='Standardized β')

points=read('point_estimates.csv',True);boot=read('bootstrap_summary.csv',True);omnibus=read('group_omnibus_tests.csv',True)
single=['ChP_to_sTREM2_to_Cognition','sTREM2_to_ChP_to_Cognition','ChP_to_TAU_to_Cognition','TAU_to_ChP_to_Cognition',
        'ChP_to_PTAU_to_Cognition','PTAU_to_ChP_to_Cognition','ChP_to_ABETA_to_Cognition','ABETA_to_ChP_to_Cognition']
def effect_row(model,effect,group):
    return boot[(boot.model==model)&(boot.definition=='raw_fullmetric')&(boot.effect==effect+'_'+group)].iloc[0]
for n,models in [(13,single[:4]),(14,single[4:])]:
    fig,axes=plt.subplots(2,2,figsize=(10,6))
    for ax,model in zip(axes.flat,models):
        for i,g in enumerate(GROUPS):
            r=effect_row(model,'indirect',g)
            ax.plot([r.ci_unstd_low,r.ci_unstd_high],[i,i],color=COLORS[i])
            ax.scatter(r.point_unstd,i,s=55,edgecolors=COLORS[i],facecolors=COLORS[i] if r.MLR_p<.05 else 'white',zorder=3)
        ax.axvline(0,color='#cccccc',ls='--');ax.set_yticks(range(3),['CN','MCI','Dementia']);ax.invert_yaxis()
        ax.set_title(label(model),fontsize=9);ax.set_xlabel('Unstandardized indirect effect (b)')
    fig.tight_layout();save(fig,n)

parallel=['ChP_to_sTREM2_PTAU_ABETA_parallel','ChP_to_sTREM2_TAU_ABETA_parallel']
def diagram(ax,model,group):
    overall=group=='Overall';ident=('raw' if overall else 'raw_fullmetric')+'__parallel__'+model+'__'+('Overall' if overall else 'Multigroup')
    p=points[points.id==ident];suffix='' if overall else '_'+group
    def value(name):
        rr=p[(p.label==name+suffix)&(p.op=='~')]
        if len(rr)!=1: raise ValueError((ident,name,len(rr)))
        return rr.iloc[0]
    positions={'x':(.1,.5),'mediator_1':(.5,.85),'mediator_2':(.5,.5),'mediator_3':(.5,.15),'Cog':(.9,.5)}
    names={'x':'ChP/ICV','mediator_1':'sTREM2','mediator_2':'P-Tau' if 'PTAU' in model else 'Tau','mediator_3':'CSF Aβ42','Cog':'Cognition'}
    for j in range(1,4):
        for source,target,param,offset in [('x','mediator_'+str(j),'a'+str(j),.025),('mediator_'+str(j),'Cog','b'+str(j),-.025)]:
            r=value(param);a=np.array(positions[source]);b=np.array(positions[target]);v=b-a
            ax.annotate('',xy=b-v*.19,xytext=a+v*.19,arrowprops={'arrowstyle':'->','lw':1.1,'color':'#555555'})
            mid=(a+b)/2;ax.text(mid[0],mid[1]+offset,f'{r["std.all"]:.3f}'+stars(r.pvalue),ha='center',fontsize=8,bbox={'facecolor':'white','edgecolor':'none','pad':1})
    r=value('c_prime')
    ax.annotate('',xy=(.85,.44),xytext=(.15,.44),arrowprops={'arrowstyle':'->','connectionstyle':'arc3,rad=.65','color':'#777777'})
    ax.text(.5,-.04,'Direct: '+f'{r["std.all"]:.3f}'+stars(r.pvalue),ha='center',fontsize=8)
    for key,pos in positions.items():ax.text(*pos,names[key],ha='center',va='center',fontsize=9,bbox={'boxstyle':'round,pad=.35','fc':'#E3EFF9','ec':'#6E98B6'})
    effect='total_indirect'+suffix;r=p[(p.lhs==effect)&(p.op==':=')].iloc[0]
    ax.text(.5,-.14,f'Total indirect: {r["std.all"]:.3f}'+stars(r.pvalue),ha='center',fontsize=9)
    ax.set_xlim(-.03,1.03);ax.set_ylim(-.2,1.03);ax.axis('off');ax.set_title(group.replace('AD','Dementia')+' | '+('P-Tau' if 'PTAU' in model else 'Tau'),fontsize=11)
for n,groups in [(15,['Overall','CN']),(16,['MCI','AD'])]:
    fig,axes=plt.subplots(2,2,figsize=(11,7))
    for i,g in enumerate(groups):
        for j,model in enumerate(parallel):diagram(axes[i,j],model,g)
    fig.text(.5,.015,'* $P$ < .05, ** $P$ < .01, *** $P$ < .001',ha='center',fontsize=9);fig.tight_layout(rect=(0,.04,1,1));save(fig,n)
for n,model in [(17,parallel[0]),(18,parallel[1])]:
    fig,ax=plt.subplots(figsize=(8,5));diagram(ax,model,'Overall')
    fig.text(.5,.015,'* $P$ < .05, ** $P$ < .01, *** $P$ < .001',ha='center',fontsize=9);save(fig,n)
matrix=[];qs=[]
for model in single:
    matrix.append([effect_row(model,'indirect',g).point_unstd for g in GROUPS])
    qs.append(omnibus[(omnibus.model==model)&(omnibus.effect=='indirect')].iloc[0].q_BH16)
heat(np.array(matrix),[label(m) for m in single],['CN','MCI','Dementia'],19,q=qs)
matrix=[];sym=[];qs=[];labels=[]
for model in parallel:
    for effect,title in [('indirect_1','sTREM2'),('indirect_2','P-Tau' if 'PTAU' in model else 'Tau'),('indirect_3','CSF Aβ42'),('total_indirect','Total indirect')]:
        rs=[effect_row(model,effect,g) for g in GROUPS]
        matrix.append([r.point_unstd for r in rs]);sym.append([stars(r.MLR_p) for r in rs])
        qs.append(omnibus[(omnibus.model==model)&(omnibus.effect==effect)].iloc[0].q_BH16)
        labels.append(('P-Tau' if 'PTAU' in model else 'Tau')+' model | '+title)
heat(np.array(matrix),labels,['CN','MCI','Dementia'],20,symbols=sym)
heat(np.array(matrix),labels,['CN','MCI','Dementia'],21,q=qs)
t=read('serial_crossmodel_ratio_diagnostic.csv',True);t=t[(t.definition=='raw')&(t.scale=='est')]
fig,ax=plt.subplots(figsize=(6,4));x=np.arange(len(t));v=t.point_pct.to_numpy()
ax.bar(x,v,color=COLORS[:2],width=.5);ax.errorbar(x,v,yerr=np.array([v-t.ci_low,t.ci_high-v]),fmt='none',color='black',capsize=4)
ax.set_xticks(x,[label(v) for v in t.marker]);ax.set_ylabel('Serial / single-path indirect effect (%)')
for i,value in enumerate(v):ax.text(i,value/2,f'{value:.1f}%',ha='center')
fig.tight_layout();save(fig,22)
