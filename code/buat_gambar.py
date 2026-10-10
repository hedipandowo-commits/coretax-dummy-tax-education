import json, os, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__));ROOT=os.path.abspath(os.path.join(HERE,'..'))
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
plt.rcParams['font.family']='DejaVu Serif'
R=json.load(open(os.path.join(ROOT,'results','hasil_analisis.json')));B=R['boot']['B'];R2=R['pls']['R2'];M=R['boot']['M']
f=lambda x,dp=2:f"{x:.{dp}f}".replace('.',',').replace('-','−')
def sig(k): b=B[k]; return not(b['lo']<=0<=b['hi'])
# ---- Fig 2
fig=plt.figure(figsize=(7.2,4.6),dpi=220);ax=fig.add_axes([0,0,1,1]);ax.set_xlim(0,1584);ax.set_ylim(1012,0);ax.axis('off')
boxes={'TCSE':(80,481,356,597,'Tax Compliance\nSelf-Efficacy (TCSE)'),'PEOU':(469,271,745,386,'Perceived Ease\nof Use (PEOU)\nR² = %s'%f(R2['PEOU'])),
'PU':(469,692,745,807,'Perceived\nUsefulness (PU)\nR² = %s'%f(R2['PU'])),'ATU':(870,481,1146,597,'Attitude Toward\nUsing (ATU)\nR² = %s'%f(R2['ATU'])),'BI':(1235,481,1511,597,'Behavioral\nIntention (BI)\nR² = %s'%f(R2['BI']))}
for k,(x0,y0,x1,y1,t) in boxes.items():
    ax.add_patch(FancyBboxPatch((x0,y0),x1-x0,y1-y0,boxstyle='round,pad=0,rounding_size=12',fc='white',ec='black',lw=1.2))
    ax.text((x0+x1)/2,(y0+y1)/2,t,ha='center',va='center',fontsize=7.8 if k=='TCSE' else 8.6)
def arrow(p,q,lab,k,pos=0.5,off=(0,0),rad=0):
    s=sig(k);ls='-' if s else (0,(3,3))
    a=FancyArrowPatch(p,q,arrowstyle='-|>',mutation_scale=14,lw=1.6,ls=ls,color='black',connectionstyle=f'arc3,rad={rad}',shrinkA=0,shrinkB=0)
    ax.add_patch(a)
    return a
def lab(x,y,k):
    ax.text(x,y,f(B[k]['est'])+('*' if sig(k) else ' ns'),ha='center',va='center',fontsize=9.5,bbox=dict(fc='white',ec='none',pad=1.5))
arrow((337,474),(487,394),'', 'TCSE>PEOU');lab(412,432,'TCSE>PEOU')
arrow((337,603),(487,684),'', 'TCSE>PU');lab(412,642,'TCSE>PU')
arrow((607,402),(607,674),'', 'PEOU>PU');lab(607,538,'PEOU>PU')
arrow((730,393),(884,474),'', 'PEOU>ATU');lab(807,432,'PEOU>ATU')
arrow((730,686),(884,604),'', 'PU>ATU');lab(807,642,'PU>ATU')
arrow((748,765),(1342,598),'', 'PU>BI');lab(1046,682,'PU>BI')
arrow((1160,539),(1220,539),'', 'ATU>BI');lab(1190,498,'ATU>BI')
# dashed curves
for k,(x0,x1,h,lx,ly) in {'TCSE>ATU':(253,971,205,604,206),'TCSE>BI':(183,1374,90,786,90)}.items():
    t=np.linspace(0,1,100);xs=x0+(x1-x0)*t
    y0=474;y1=474
    ys=y0-4*h*(1-t)*t*1.0 if False else y0-(y0-h)*4*t*(1-t)
    ax.plot(xs[:-3],ys[:-3],ls=(0,(4,3)) if not sig(k) else '-',color='black',lw=1.4)
    ax.add_patch(FancyArrowPatch((xs[-4],ys[-4]),(xs[-1],y1),arrowstyle='-|>',mutation_scale=12,lw=1.4,color='black'))
    lab(lx,ly,k)
mi=M['inter'];ms='*' if not(mi['lo']<=0<=mi['hi']) else ' ns'
ax.text(792,890,'Moderasi TCSE × PEOU → ATU: β = %s%s'%(f(mi['est']),ms),ha='center',fontsize=9.5,style='italic')
ax.text(792,946,'Garis penuh = jalur signifikan; putus-putus = tidak signifikan (ns). * CI 95% bootstrap (5.000) tidak memuat 0.',ha='center',fontsize=7.6,style='italic')
fig.savefig(os.path.join(ROOT,'figures','gambar2_model_tax_tam.png'),dpi=220)
# ---- Fig 3
fig,ax=plt.subplots(figsize=(4.6,3.3),dpi=220)
x=np.linspace(-2,2,50)
for key,ls,lbl in (('s_lo','-','TCSE rendah (−1 SD)'),('s_mid',':','TCSE rata-rata'),('s_hi','--','TCSE tinggi (+1 SD)')):
    s=M[key]['est'];ax.plot(x,s*x,color='black',ls=ls,lw=1.6,label=f'{lbl}: kemiringan = {f(s)}')
ax.set_xlabel('PEOU (skor standar)',fontsize=9);ax.set_ylabel('ATU (skor standar)',fontsize=9)
ax.spines[['top','right']].set_visible(False);ax.tick_params(labelsize=8)
ax.legend(frameon=False,fontsize=6.8,loc='upper left')
fig.tight_layout();fig.savefig(os.path.join(ROOT,'figures','gambar3_simple_slopes.png'),dpi=220)
