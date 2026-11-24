import sys
import numpy as np,re,glob,os
T0='/gpfs/hpc-home/projects/dtu_00075/people/sitali/FoamWaveTank_Training_test/'
def body(f):
    s=open(f).read();i=s.find('internalField');return s,i
def field(f,nc=1):
    s,i=body(f);h=s[i:min(i+120,len(s))]
    if 'nonuniform' not in h.split('\n')[0]:return None
    a=s.index('(',i);b=s.index('\n)',a)
    v=np.fromstring(s[a+1:b].replace('(',' ').replace(')',' '),sep=' ')
    return v.reshape(-1,nc) if nc>1 else v
def listbody(f):
    s=open(f).read();m=re.search(r'\n(\d+)\s*\n\(',s);return s[m.end():s.index('\n)',m.end())],int(m.group(1))
def pts(f):
    s,n=listbody(f);return np.fromstring(s.replace('(',' ').replace(')',' '),sep=' ').reshape(n,3)
def cells(p):
    q=p+'/constant/polyMesh/';F=[list(map(int,m.split())) for m in re.findall(r'\d+\(([^)]*)\)',listbody(q+'faces')[0])]
    C={}
    for n in ('owner','neighbour'):
        for f,c in enumerate(np.fromstring(listbody(q+n)[0],sep=' ',dtype=int)):C.setdefault(c,set()).update(F[f])
    return np.array([sorted(C[i]) for i in range(len(C))])
src=open(T0+'steps12_22Nov/build/plot_day3.py').read().split('\nRUNS = [')[0]
g={'__file__':T0+'steps12_22Nov/build/plot_day3.py','__name__':'c'};exec(src,g);N=g['NX']

import re as _re
RUNS=[('plain fine','meshes/plain_fine'),('plain coarse','d4quick_pCSu/plain_coarse'),('SW test 2','d4test2_fn11/case'),('SW M10','zonefix/M10')]
RUNS+=[(a.split('=')[0],a.split('=')[1]) for a in sys.argv[1:]]          # label=path relative to T0
RUNS=[r for r in RUNS if os.path.exists(T0+r[1]+'/23')]
per=5;a=23-per*g['T'];K=('surge','heave','pitch','tension')
Mo={n:g['numbers'](g['cfd'](T0+p),a,per) for n,p in RUNS}
print('MOTION amplitude, last %d periods (%.2f-23 s): error %% vs plain fine | vs experiment'%(per,a))
for n,_ in RUNS:
    print('%-13s'%n+''.join('%9.1f%% |%7.1f%%   '%(100*(Mo[n][k]['amp']/Mo['plain fine'][k]['amp']-1),100*(Mo[n][k]['amp']/N[k]['amp']-1)) for k in K))
X=np.array((-1.25,-1.,-.25,.25,1.,1.25))
tl=[round(a+.1*i,1) for i in range(0,int(per*g['T']*10)+2)];tl=[t for t in tl if a<=t<=23];L={}
for n,p in RUNS:
    P=T0+p;C=cells(P);P0=pts(P+'/constant/polyMesh/points')
    xc0=P0[C].mean(1)[:,0]; zc0=P0[C].mean(1)[:,2]
    key=np.round(xc0,6); ux=np.unique(key); cid=np.searchsorted(ux,key)
    nPer=np.bincount(cid); full=nPer==nPer.max()
    M={}
    for t in tl:
        d='%g'%t
        if not os.path.exists(f'{P}/{d}/alpha.water'):continue
        fp=f'{P}/{d}/polyMesh/points';Pt=(pts(fp) if os.path.exists(fp) else P0)[C]
        z=Pt[:,:,2];dz=z.max(1)-z.min(1);xt=Pt[:,:,0].mean(1)
        q=field(f'{P}/{d}/alpha.water')
        eta=np.bincount(cid,q*dz,len(ux))+np.array([z[cid==k].min() for k in range(len(ux))])
        xcol=np.bincount(cid,xt,len(ux))/nPer
        sel=full; o=np.argsort(xcol[sel])
        v=np.interp(X,xcol[sel][o],eta[sel][o])        # static part: linear in x
        for j,x in enumerate(X):
            if abs(x)<0.5:                                 # moving part: the two columns that straddle x at t=0,
                kr=np.searchsorted(ux,x); kl=kr-1          # linear weights from their initial x (columns move with the mesh)
                w=(x-ux[kl])/(ux[kr]-ux[kl]); v[j]=(1-w)*eta[kl]+w*eta[kr]
        M[t]=v
    L[n]=M
print('WAVE total eta vs plain fine, rms mm; |x|=1, 1.25 linear between static columns, |x|=0.25 linear between the straddling initial columns, moving with the mesh; x =',tuple(float(x) for x in X))
for n,_ in RUNS:
    if n=='plain fine':continue
    q=[t for t in L[n] if t in L['plain fine']];d=np.array([L[n][t]-L['plain fine'][t] for t in q])*1e3
    r=np.sqrt((d*d).mean(0));print('%-13s'%n,' '.join('%6.1f'%v for v in r),'| box +-0.25 %.1f, x+-1 %.1f, coarse +-1.25 %.1f  (n=%d)'%(r[[2,3]].mean(),r[[1,4]].mean(),r[[0,5]].mean(),len(q)))
