import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools'))
from paths import arguments, font_path
from pathlib import Path
import csv,json,shutil,html
from PIL import Image,ImageDraw,ImageFont
args=arguments()
src=args.results/'incremental/incremental_results.csv'
rows=[r for r in csv.DictReader(src.open(encoding='utf-8-sig')) if r['specification']=='primary']
assert [r['biomarker'] for r in rows]==['sTREM2','Tau','P-Tau']
im=Image.new('RGB',(4800,2460),'white');d=ImageDraw.Draw(im)
svg=['<svg xmlns="http://www.w3.org/2000/svg" width="4800" height="2460" viewBox="0 0 4800 2460"><rect width="4800" height="2460" fill="white"/>']
font=lambda s,b=False:ImageFont.truetype(font_path(bold=b),int(s))
def text(x,y,s,size=60,bold=False,color='black',anchor='middle'):
 f=font(size,bold);bb=d.textbbox((0,0),s,font=f);xx=x-(bb[2]-bb[0])/2 if anchor=='middle' else x
 if s.startswith('P '):
  fi=ImageFont.truetype(font_path(italic=True),int(size));xx=x-(fi.getlength('P')+f.getlength(s[1:]))/2
  d.text((xx,y),'P',font=fi,fill=color);d.text((xx+fi.getlength('P'),y),s[1:],font=f,fill=color)
 else:d.text((xx,y),s,font=f,fill=color)
 svg.append(f'<text xml:space="preserve" x="{x}" y="{y+size*.92}" text-anchor="{anchor}" font-family="Arial" font-size="{size}" font-weight="{700 if bold else 400}" fill="{color}">{('<tspan font-style="italic">P</tspan>'+html.escape(s[1:])) if s.startswith('P ') else html.escape(s)}</text>')
def line(x1,y1,x2,y2,color='#333333',width=4):
 d.line((x1,y1,x2,y2),fill=color,width=width);svg.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{width}"/>')
def rect(x,y,w,h,fill,stroke=None):
 d.rectangle((x,y,x+w,y+h),fill=fill,outline=stroke,width=3)
 svg.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{fill}" stroke="{stroke or "none"}" stroke-width="3"/>')
ci={(r['biomarker'],r['model']):(float(r['lower']),float(r['upper'])) for r in csv.DictReader((args.results/'incremental/figure4_r2_bootstrap_ci.csv').open(encoding='utf8'))}
cols=['#3478B8','#E8751A','#CE4C62'];pales=['#E3EFF9','#FCEAD8','#F9E2E6']
top,bottom=500,2200;sy=lambda v:bottom-v/.40*(bottom-top)
for i,(r,col,pale) in enumerate(zip(rows,cols,pales)):
 left=400+i*1580;right=left+1080;cx=(left+right)/2
 rect(left-320,60,1480,2340,'white',pale);rect(left-320,60,1480,220,pale)
 text(left+420,125,f'{chr(65+i)}. '+('Total Tau' if i==1 else r['biomarker']),80,True)
 line(left,top,left,bottom);line(left,bottom,right,bottom)
 for j in range(9):
  val=j*.05;y=sy(val);line(left-15,y,left,y);text(left-100,y-28,f'{val:.2f}',55)
 xs=[left+270,left+810];vs=[float(r['base_r2']),float(r['extended_r2'])]
 for k,(x,v,c) in enumerate(zip(xs,vs,['#C9C9C9',col])):
  yy=sy(v);rect(x-180,yy,360,bottom-yy,c)
  text(x,yy+(bottom-yy)/2-32,f'{v:.4f}',65,True,'black')
 for k,x in enumerate(xs):
  low,high=ci[(r['biomarker'],['base','extended'][k])]
  line(x,sy(low),x,sy(high),width=6)
  line(x-28,sy(low),x+28,sy(low),width=6)
  line(x-28,sy(high),x+28,sy(high),width=6)
 by=sy(max(ci[(r['biomarker'],m)][1] for m in ['base','extended']))-90
 line(xs[0],by+42,xs[0],by);line(xs[0],by,xs[1],by);line(xs[1],by,xs[1],by+42)
 text(cx,by-310,f'ΔR² = {float(r["delta_r2"]):.4f}',72,True)
 text(cx,by-210,f'({100*float(r["delta_r2"]):.2f} percentage points)',60)
 text(xs[0],2240,'Base model',61);text(xs[0],2320,'(covariates)',53)
 text(xs[1],2240,'Extended model',61);text(xs[1],2320,'(+ ChP/ICV)',53)
 text(cx,by-110,'P = '+(f'{float(r["holm_p"]):.5f}' if float(r['holm_p'])<.001 else f'{float(r["holm_p"]):.4f}').lstrip('0'),57)
# Each panel has its own axis title, with separate gutters for title and ticks.
for i in range(3):
 left=400+i*1580;label='Explained variance (R²)';f=font(61)
 box=f.getbbox(label);lw=box[2]-box[0];lh=box[3]-box[1]
 lab=Image.new('RGBA',(lw+8,lh+8),(255,255,255,0))
 ImageDraw.Draw(lab).text((4-box[0],4-box[1]),label,font=f,fill='black')
 lab=lab.rotate(90,expand=True)
 xc=left-245;yc=(top+bottom)/2
 im.paste(lab,(round(xc-lab.width/2),round(yc-lab.height/2)),lab)
 svg.append(f'<text xml:space="preserve" x="{xc}" y="{yc}" text-anchor="middle" dominant-baseline="middle" transform="rotate(-90 {xc} {yc})" font-family="Arial" font-size="61">Explained variance (R²)</text>')

im.save(args.output/'Figure_4.png',dpi=(600,600))
im.save(args.output/'Figure_4.jpg',quality=97,subsampling=0,dpi=(600,600))
(args.output/'Figure_4.svg').write_text('\n'.join(svg+['</svg>']),encoding='utf8')
