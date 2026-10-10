import numpy as np, pandas as pd, json, sys
from scipy import stats
import os, argparse, warnings
warnings.filterwarnings('ignore')
ap=argparse.ArgumentParser(description='Analisis statistik artikel Sistem Dummy Coretax')
ap.add_argument('--data',default=None);ap.add_argument('--out',default=None);ap.add_argument('--boot',type=int,default=5000)
args,_=ap.parse_known_args()
HERE=os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
ROOT=os.path.abspath(os.path.join(HERE,'..')) if os.path.basename(HERE)=='code' else HERE
DATA=args.data or os.path.join(ROOT,'data','data_respons_individual.csv')
OUT=args.out or os.path.join(ROOT,'results')
os.makedirs(OUT,exist_ok=True)
d=pd.read_csv(DATA)
C=['TCSE','PEOU','PU','ATU','BI']
ind={c:[f'{c}{i}' for i in range(1,7)] for c in C}
paths=[('TCSE','PEOU'),('TCSE','PU'),('PEOU','PU'),('PEOU','ATU'),('PU','ATU'),('PU','BI'),('ATU','BI'),('TCSE','ATU'),('TCSE','BI')]
R={}
# ---------- H1
pre=d.PRE/30*100; post=d.POST/30*100
diff=post-pre; n=len(d)
sw=stats.shapiro(diff)
tt=stats.ttest_rel(post,pre,alternative='greater')
wl=stats.wilcoxon(post,pre,alternative='greater',method='approx',correction=False)
z=getattr(wl,'zstatistic',None)
dz=diff.mean()/diff.std(ddof=1)
se=diff.std(ddof=1)/np.sqrt(n); tc=stats.t.ppf(.975,n-1)
ng=(post-pre)/(100-pre)
rb=(np.sum(diff>0)-np.sum(diff<0))/n
R['H1']=dict(pre_m=pre.mean(),pre_sd=pre.std(ddof=1),post_m=post.mean(),post_sd=post.std(ddof=1),
 pct=(post.mean()-pre.mean())/pre.mean()*100,diff=diff.mean(),ci=[diff.mean()-tc*se,diff.mean()+tc*se],
 t=tt.statistic,p=tt.pvalue,dz=dz,sw_W=sw.statistic,sw_p=sw.pvalue,wil_z=z,wil_p=wl.pvalue,rb=rb,
 npos=int((diff>0).sum()),nneg=int((diff<0).sum()),nzero=int((diff==0).sum()),
 ng_m=ng.mean(),ng_sd=ng.std(ddof=1),ng_hi=int((ng>0.7).sum()),ng_mid=int(((ng>=0.3)&(ng<=0.7)).sum()),ng_lo=int((ng<0.3).sum()),
 raw_pre=d.PRE.mean(),raw_post=d.POST.mean())
# ---------- SUS
odd=[f'SUS{i}' for i in (1,3,5,7,9)];even=[f'SUS{i}' for i in (2,4,6,8,10)]
sus=((d[odd]-1).sum(1)+(5-d[even]).sum(1))*2.5
d['SUS']=sus
def one(x,ref):
    x=np.asarray(x);t=stats.ttest_1samp(x,ref,alternative='greater')
    w=stats.wilcoxon(x-ref,alternative='greater',method='approx',correction=False)
    return dict(m=x.mean(),sd=x.std(ddof=1),t=t.statistic,p=t.pvalue,z=getattr(w,'zstatistic',None),wp=w.pvalue,d=(x.mean()-ref)/x.std(ddof=1))
se=sus.std(ddof=1)/np.sqrt(n)
R['SUS']=dict(**one(sus,68),ci=[sus.mean()-tc*se,sus.mean()+tc*se],sw_p=stats.shapiro(sus).pvalue,sw_W=stats.shapiro(sus).statistic,vs70=one(sus,70),min=sus.min(),max=sus.max())
for c in ('PU','PEOU'):
    s=d[ind[c]].mean(1); d[c]=s
    R[c]=dict(**one(s,3.0),sw_p=stats.shapiro(s).pvalue)
for c in ('ATU','BI','TCSE'):
    d[c]=d[ind[c]].mean(1); R[c+'_mean']=(d[c].mean(),d[c].std(ddof=1))
R['corr_sus_peou']=dict(r=stats.pearsonr(sus,d.PEOU),rho=stats.spearmanr(sus,d.PEOU))
R['corr_sus_peou']={'r':R['corr_sus_peou']['r'][0],'rp':R['corr_sus_peou']['r'][1],'rho':R['corr_sus_peou']['rho'][0],'rhop':R['corr_sus_peou']['rho'][1]}
# ---------- PLS
def zs(a): return (a-a.mean(0))/a.std(0,ddof=1)
def pls(X,tol=1e-7,maxit=300):
    # X dict construct-> standardized matrix
    W={c:np.ones(X[c].shape[1]) for c in C}
    LV={c:None for c in C}
    def score(c,w):
        s=X[c]@w; return (s-s.mean())/s.std(ddof=1)
    for c in C: LV[c]=score(c,W[c])
    pred={c:[a for a,b in paths if b==c] for c in C}
    succ={c:[b for a,b in paths if a==c] for c in C}
    for it in range(maxit):
        Wold={c:W[c].copy() for c in C}
        inner={}
        for c in C:
            Z=np.zeros(len(LV[c]))
            if pred[c]:
                P=np.column_stack([LV[p] for p in pred[c]])
                b=np.linalg.lstsq(P,LV[c],rcond=None)[0]
                for p,bb in zip(pred[c],b): Z+=bb*LV[p]
            for s_ in succ[c]: Z+=np.corrcoef(LV[c],LV[s_])[0,1]*LV[s_]
            inner[c]=Z
        for c in C:
            w=X[c].T@inner[c]/(len(inner[c])-1)
            W[c]=w
            LV[c]=score(c,w)
            # rescale w so LV unit variance
            W[c]=w/ (X[c]@w).std(ddof=1)
        if max(np.abs(W[c]-Wold[c]).max() for c in C)<tol: break
    return W,LV
def paths_est(LV):
    B={}
    for c in C:
        pr=[a for a,b in paths if b==c]
        if pr:
            P=np.column_stack([LV[p] for p in pr]);b=np.linalg.lstsq(P,LV[c],rcond=None)[0]
            for p,bb in zip(pr,b): B[(p,c)]=bb
    return B
def r2s(LV):
    out={}
    for c in C:
        pr=[a for a,b in paths if b==c]
        if pr:
            P=np.column_stack([LV[p] for p in pr]);b=np.linalg.lstsq(P,LV[c],rcond=None)[0]
            e=LV[c]-P@b;out[c]=1-e.var(ddof=1)/LV[c].var(ddof=1)
    return out
def indirect(B):
    b=lambda a,c:B.get((a,c),0.0)
    ie={}
    ie['PEOU>PU>ATU']=b('PEOU','PU')*b('PU','ATU')
    ie['PU>ATU>BI']=b('PU','ATU')*b('ATU','BI')
    t={}
    t['TCSE>PEOU>PU>ATU>BI']=b('TCSE','PEOU')*b('PEOU','PU')*b('PU','ATU')*b('ATU','BI')
    t['TCSE>PEOU>PU>BI']=b('TCSE','PEOU')*b('PEOU','PU')*b('PU','BI')
    t['TCSE>PEOU>ATU>BI']=b('TCSE','PEOU')*b('PEOU','ATU')*b('ATU','BI')
    t['TCSE>PU>ATU>BI']=b('TCSE','PU')*b('PU','ATU')*b('ATU','BI')
    t['TCSE>PU>BI']=b('TCSE','PU')*b('PU','BI')
    t['TCSE>ATU>BI']=b('TCSE','ATU')*b('ATU','BI')
    ie['TCSE>BI total']=sum(t.values())
    ie.update(t)
    return ie
def mod(LV):
    tc=LV['TCSE'];pe=LV['PEOU'];pu=LV['PU'];at=LV['ATU']
    inter=tc*pe
    P=np.column_stack([pe,pu,tc,inter]);b=np.linalg.lstsq(P,at,rcond=None)[0]
    P0=P[:,:3];b0=np.linalg.lstsq(P0,at,rcond=None)[0]
    r2=lambda P,b:1-(at-P@b).var(ddof=1)/at.var(ddof=1)
    return dict(peou=b[0],pu=b[1],tcse=b[2],inter=b[3],s_lo=b[0]-b[3],s_mid=b[0],s_hi=b[0]+b[3],r2_0=r2(P0,b0),r2_1=r2(P,b))
Xraw={c:d[ind[c]].values.astype(float) for c in C}
X={c:zs(Xraw[c]) for c in C}
W,LV=pls(X)
B=paths_est(LV);R2=r2s(LV)
load={c:np.array([np.corrcoef(X[c][:,i],LV[c])[0,1] for i in range(6)]) for c in C}
# measurement
meas={}
for c in C:
    l=load[c];k=len(l)
    alpha=k/(k-1)*(1-k/(k+2*np.sum([np.corrcoef(X[c][:,i],X[c][:,j])[0,1] for i in range(k) for j in range(i+1,k)]))) 
    cr=l.sum()**2/(l.sum()**2+(1-l**2).sum());ave=(l**2).mean()
    w=W[c]/np.sqrt(W[c]@np.corrcoef(X[c].T)@W[c])
    S=np.corrcoef(X[c].T);wwT=np.outer(w,w)
    off=lambda M:M-np.diag(np.diag(M))
    rhoA=(w@w)**2*(w@off(S)@w)/(w@off(wwT)@w)
    meas[c]=dict(load_min=l.min(),load_max=l.max(),loads=l.round(3).tolist(),alpha=alpha,rhoA=rhoA,cr=cr,ave=ave)
LVdf=pd.DataFrame(LV)
LVc=LVdf.corr()
allX=np.column_stack([X[c] for c in C]);blk=np.repeat(range(5),6)
IC=np.corrcoef(allX.T)
def htmt(IC):
    H=np.zeros((5,5))
    for a in range(5):
        for b in range(5):
            if a==b: continue
            ma=IC[np.ix_(blk==a,blk==b)].mean()
            mono=lambda k:(IC[np.ix_(blk==k,blk==k)][np.triu_indices(6,1)]).mean()
            H[a,b]=ma/np.sqrt(mono(a)*mono(b))
    return H
H=htmt(IC)
# VIF
vif={}
for c in C:
    pr=[a for a,b in paths if b==c]
    if len(pr)>=2:
        for p in pr:
            others=[q for q in pr if q!=p]
            P=np.column_stack([LV[q] for q in others]);bb=np.linalg.lstsq(P,LV[p],rcond=None)[0]
            r2=1-(LV[p]-P@bb).var(ddof=1)/LV[p].var(ddof=1);vif[f'{p}->{c}']=1/(1-r2)
# f2
f2={}
for (a,c) in paths:
    pr=[p for p,q in paths if q==c]
    def r2f(prs):
        P=np.column_stack([LV[p] for p in prs]);b=np.linalg.lstsq(P,LV[c],rcond=None)[0]
        return 1-(LV[c]-P@b).var(ddof=1)/LV[c].var(ddof=1)
    full=r2f(pr);red=r2f([p for p in pr if p!=a]) if len(pr)>1 else 0
    f2[f'{a}->{c}']=(full-red)/(1-full)
# SRMR
Bm=np.zeros((5,5)) # row=target
idx={c:i for i,c in enumerate(C)}
for (a,c),v in B.items(): Bm[idx[c],idx[a]]=v
Psi=np.zeros((5,5))
for c in C:
    Psi[idx[c],idx[c]]=1-R2[c] if c in R2 else 1.0
Ainv=np.linalg.inv(np.eye(5)-Bm);SigLV=Ainv@Psi@Ainv.T
imp=np.zeros_like(IC)
lam=np.concatenate([load[c] for c in C])
for i in range(30):
    for j in range(30):
        if i==j: imp[i,j]=1
        elif blk[i]==blk[j]: imp[i,j]=lam[i]*lam[j]
        else: imp[i,j]=lam[i]*lam[j]*SigLV[blk[i],blk[j]]
res=IC-imp;tri=np.tril_indices(30)
srmr=np.sqrt((res[tri]**2).mean())
M0=mod(LV)
R['pls']=dict(B={f'{a}>{c}':v for (a,c),v in B.items()},R2=R2,meas=meas,fl=LVc.round(3).values.tolist(),sqrtAVE={c:np.sqrt(meas[c]['ave']) for c in C},
  htmt=H.round(3).tolist(),vif=vif,f2=f2,srmr=srmr,ind=indirect(B),mod=M0)
# ---------- bootstrap
rng=np.random.default_rng(20260110)  # seed tetap agar hasil bootstrap dapat direproduksi
NB=args.boot
bB={f'{a}>{c}':[] for a,c in paths};bI={k:[] for k in indirect(B)};bM={k:[] for k in M0};bH=[];bLoad={c:[] for c in C};bR2={c:[] for c in R2}
fails=0
for it in range(NB):
    ix=rng.integers(0,n,n)
    Xb={c:zs(Xraw[c][ix]) for c in C}
    try:
        Wb,LVb=pls(Xb)
    except Exception as e:
        fails+=1;continue
    # sign alignment with original weights
    for c in C:
        if np.corrcoef(LVb[c],LV[c][ix])[0,1]<0:
            Wb[c]=-Wb[c];LVb[c]=-LVb[c]
    Bb=paths_est(LVb)
    for (a,c),v in Bb.items(): bB[f'{a}>{c}'].append(v)
    for k,v in indirect(Bb).items(): bI[k].append(v)
    for k,v in mod(LVb).items(): bM[k].append(v)
    ICb=np.corrcoef(np.column_stack([Xb[c] for c in C]).T)
    bH.append(htmt(ICb))
    for c,v in r2s(LVb).items(): bR2[c].append(v)
print('boot fails',fails)
def summ(est,arr):
    arr=np.array(arr);se=arr.std(ddof=1);t=est/se
    return dict(est=est,se=se,t=t,p=2*stats.t.sf(abs(t),n-1) if False else 2*stats.t.sf(abs(t),95),lo=np.percentile(arr,2.5),hi=np.percentile(arr,97.5))
R['boot']=dict(B={k:summ(B[tuple(k.split('>'))],v) for k,v in bB.items()},
  I={k:summ(R['pls']['ind'][k],v) for k,v in bI.items()},
  M={k:summ(M0[k],v) for k,v in bM.items()},
  htmt_hi=np.percentile(np.array(bH),97.5,axis=0).round(3).tolist(),
  R2={k:summ(R2[k],v) for k,v in bR2.items()})
# direct effects & VAF
b=B
R['med']=dict(
 H11=dict(direct=b[('PEOU','ATU')],ind=R['pls']['ind']['PEOU>PU>ATU']),
 H12=dict(direct=b[('PU','BI')],ind=R['pls']['ind']['PU>ATU>BI']),
 H13=dict(direct=b[('TCSE','BI')],ind=R['pls']['ind']['TCSE>BI total']))
for k,v in R['med'].items(): v['vaf']=v['ind']/(v['ind']+v['direct'])*100
json.dump(R,open(os.path.join(OUT,'hasil_analisis.json'),'w'),indent=1,default=lambda o:o.tolist() if hasattr(o,'tolist') else float(o))
print('done')
