import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from paths import arguments, font_path
from pathlib import Path
import csv,json,math,shutil,hashlib
import numpy as np
from PIL import Image,ImageDraw,ImageFont
args=arguments(needs_data=True)
FIG=args.output
read=lambda p:list(csv.DictReader(p.open(encoding='utf-8-sig')))
rows=read(args.data)
counts={g:sum(r['S_DX_label']==g for r in rows) for g in ['CN','MCI','AD']}
counts['Overall']=len(rows)
stats=read(args.results/'association/combined_regression_stage_panel_stats.csv')
COLORS={'CN':(31,90,138),'MCI':(232,117,26),'AD':(185,29,58)}
CONFIG=[('A','sTREM2','MSD_STREM2CORRECTED'),('B','Total Tau','TAU'),('C','P-Tau','PTAU'),('D','Aβ42','S_ABETA')]
font=lambda s,b=False:ImageFont.truetype(font_path(bold=b),int(s))
# Continued-fraction evaluation of the regularized incomplete beta for exact t intervals.
def betacf(a,b,x):
 qab=a+b;qap=a+1;qam=a-1;c=1.;d=1-qab*x/qap;d=1/max(d,1e-30);h=d
 for m in range(1,501):
  m2=2*m;aa=m*(b-m)*x/((qam+m2)*(a+m2));d=1+aa*d;c=1+aa/c
  if abs(d)<1e-30:d=1e-30
  if abs(c)<1e-30:c=1e-30
  d=1/d;h*=d*c;aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2));d=1+aa*d;c=1+aa/c
  if abs(d)<1e-30:d=1e-30
  if abs(c)<1e-30:c=1e-30
  d=1/d;v=d*c;h*=v
  if abs(v-1)<3e-14:break
 return h
def ibeta(a,b,x):
 bt=math.exp(math.lgamma(a+b)-math.lgamma(a)-math.lgamma(b)+a*math.log(x)+b*math.log1p(-x))
 return bt*betacf(a,b,x)/a if x<(a+1)/(a+b+2) else 1-bt*betacf(b,a,1-x)/b
def tcrit(df):
 lo=0.;hi=10.
 for _ in range(70):
  t=(lo+hi)/2;cdf=1-.5*ibeta(df/2,.5,df/(df+t*t))
  if cdf<.975:lo=t
  else:hi=t
 return (lo+hi)/2
assert abs(tcrit(128)-1.978670)<1e-6
records=[]
def centered(d,xy,text,size,b=False,fill='black'):
 f=font(size,b);bb=d.textbbox((0,0),text,font=f);d.text((xy[0]-(bb[2]-bb[0])/2,xy[1]),text,font=f,fill=fill)
def right(d,xy,text,size,fill='black'):
 f=font(size)
 if text.startswith('P '):
  fi=ImageFont.truetype(font_path(italic=True),int(size));x=xy[0]-fi.getlength('P')-f.getlength(text[1:])
  d.text((x,xy[1]),'P',font=fi,fill=fill);d.text((x+fi.getlength('P'),xy[1]),text[1:],font=f,fill=fill)
 else:
  bb=d.textbbox((0,0),text,font=f);d.text((xy[0]-(bb[2]-bb[0]),xy[1]),text,font=f,fill=fill)
def fit(bio,group):
 rr=[r for r in rows if group=='Overall' or r['S_DX_label']==group]
 x=np.log(np.array([float(r['ChPICV']) for r in rr]));y=np.log(np.array([float(r[bio]) for r in rr]));assert np.isfinite(x).all() and np.isfinite(y).all()
 X=np.column_stack([np.ones(len(x)),x]);coef=np.linalg.lstsq(X,y,rcond=None)[0];res=y-X@coef;s2=(res@res)/(len(x)-2);r2=1-(res@res)/sum((y-y.mean())**2)
 assert abs(r2-np.corrcoef(x,y)[0,1]**2)<1e-12
 xx=np.linspace(x.min(),x.max(),150);pred=coef[0]+coef[1]*xx;half=tcrit(len(x)-2)*np.sqrt(s2*(1/len(x)+(xx-x.mean())**2/sum((x-x.mean())**2)))
 st=next(s for s in stats if s['x']==bio and s['stage']==group);assert len(x)==int(st['n'])
 records.append({'biomarker':bio,'group':group,'n':len(x),'adjusted_beta_source':float(st['beta']),'adjusted_p_source':float(st['p_value']),'display_unadjusted_loglog_r2':r2,'display_loglog_slope':float(coef[1]),'t_critical_95':tcrit(len(x)-2)})
 return rr,x,y,xx,pred,half,r2,st
LIMITS={'MSD_STREM2CORRECTED':(400,20000,[500,1000,2000,5000,10000]),'TAU':(70,1100,[100,200,500,1000]),'PTAU':(6,120,[10,20,50,100]),'S_ABETA':(160,2200,[200,500,1000,2000])}
def panel(base,origin,bio,label,title,group='Overall',scale=1):
 # Draw at a common panel size, then uniformly scale for multi-panel supplements.
 im=Image.new('RGB',(1760,1120),'white');d=ImageDraw.Draw(im)
 rr,x,y,xx,pred,half,r2,st=fit(bio,group)
 d.text((200,30),title,font=font(60,True),fill='black')
 beta=f"β = {float(st['beta']):.3f}".replace('-','−');pv=float(st['p_value']);p='P < .001' if pv<.001 else f'P = {pv:.3f}'.replace('0.','.')
 right(d,(1570,25),beta,44)
 right(d,(1570,80),p,44)
 right(d,(1570,135),f'R² = {r2:.3f}',44)
 l,r,t,b=370,1570,200,850;xlo,xhi=-7.8,-6.;yl,yh,yt=LIMITS[bio];ylo,yhi=np.log([yl,yh])
 mx=lambda a:l+(np.asarray(a)-xlo)/(xhi-xlo)*(r-l)
 my=lambda a:b-(np.asarray(a)-ylo)/(yhi-ylo)*(b-t)
 assert x.min()>xlo and x.max()<xhi and y.min()>ylo and y.max()<yhi
 poly=list(zip(mx(xx),my(pred-half)))+list(zip(mx(xx)[::-1],my(pred+half)[::-1]))
 # Fit layer clipped to plotting rectangle to keep rare tail CI extrema in bounds.
 plot=Image.new('RGB',(r-l,b-t),'white');pd=ImageDraw.Draw(plot, 'RGBA')
 pd.polygon([(a-l,c-t) for a,c in poly],fill=(220,220,220))
 for row,px,py in zip(rr,mx(x),my(y)):
  color=(*COLORS[row['S_DX_label']],164);rad=7.5
  pd.ellipse((px-l-rad,py-t-rad,px-l+rad,py-t+rad),fill=color)
 pd.line([(a-l,c-t) for a,c in zip(mx(xx),my(pred))],fill='black',width=6)
 im.paste(plot,(l,t));d=ImageDraw.Draw(im)
 d.line((l,t,l,b,r,b),fill='black',width=4)
 for tick in [-7.5,-7.,-6.5,-6.]:
  px=float(mx(tick));d.line((px,b,px,b+14),fill='black',width=3);centered(d,(px,b+26),f'{tick:.1f}',50)
 for tick in yt:
  py=float(my(math.log(tick)));d.line((l-14,py,l,py),fill='black',width=3);right(d,(l-25,py-28),f'{tick:,}',50)
 centered(d,((l+r)/2,945),'ln(ChP/ICV)',56)
 lab=Image.new('RGBA',(850,90),(255,255,255,0));ld=ImageDraw.Draw(lab)
 centered(ld,(425,10),f'CSF {label.replace("Total Tau", "Tau")} (pg/mL)',54)
 lab=lab.rotate(90,expand=True);im.paste(lab,(85,100),lab)
 if scale!=1:im=im.resize((round(1760*scale),round(1120*scale)),Image.Resampling.LANCZOS)
 base.paste(im,origin)
def legend(im,y):
 d=ImageDraw.Draw(im);labels=[f"Cognitively normal (n = {counts['CN']})",f"MCI (n = {counts['MCI']})",f"Dementia (n = {counts['AD']})"];f=font(50)
 widths=[d.textbbox((0,0),s,font=f)[2]+100 for s in labels];gap=90;total=sum(widths)+2*gap;x=(im.width-total)//2
 for key,s,w in zip(COLORS,labels,widths):
  d.ellipse((x,y,x+36,y+36),fill=COLORS[key]);d.text((x+55,y-8),s,font=f,fill='black');x+=w+gap
def save(im,stem):
 im.save(FIG/f'{stem}.png',dpi=(600,600));im.save(FIG/f'{stem}.jpg',quality=97,subsampling=0,dpi=(600,600))
 prev=im.copy();prev.thumbnail((1400,1400));prev.save(FIG/f'{stem}_preview.jpg',quality=94)
main=Image.new('RGB',(3600,2440),'white');d=ImageDraw.Draw(main)
for i,(letter,label,bio) in enumerate(CONFIG):panel(main,(40+(i%2)*1760,45+(i//2)*1120),bio,label,f'{letter}. {label}')
legend(main,2320);save(main,'Figure_3')
for sn,configs in [(6,CONFIG[:2]),(7,CONFIG[2:])]:
 im=Image.new('RGB',(3600,4740),'white');d=ImageDraw.Draw(im)
 for j,(_,label,bio) in enumerate(configs):
  center_y=30+j*2290;centered(d,(1800,center_y),f'{chr(65+j)}. {label}',76,True)
  for k,group in enumerate(['Overall','CN','MCI','AD']):
   n=counts[group];gn={'AD':'Dementia'}.get(group,group)
   panel(im,(40+(k%2)*1760,120+j*2290+(k//2)*1080),bio,label,f'{gn} (n = {n})',group,scale=.98)
 legend(im,4600);save(im,f'Supplementary_Figure_S{sn}')
(FIG/'scatter_display_statistics.json').write_text(json.dumps({'panels':records},indent=2),encoding='utf8')
