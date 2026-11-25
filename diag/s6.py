"""Harmonic analysis of existing runs (no solver run): waves vs Stokes II target, mean loads, RAOs.
Gauges as s5 (|x|<0.5 interpolated between the straddling initial columns, moving with the mesh).
Fit eta = m + Re(A1 e^{-iwt}) + Re(A2 e^{-2iwt}) by least squares; target A1 = H/2 e^{i(kx+3pi/2)} (ESI StokesII and
StokesIIIncident both: absolute x, wavePhase 3pi/2), A2 = kH^2/4 (3-s^2)/(4s^3) e^{2i(kx+3pi/2)}, mean 0."""
import numpy as np,os,re,pickle,sys
T0='/gpfs/hpc-home/projects/dtu_00075/people/sitali/FoamWaveTank_Training_test/'
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'s5.py')).read().split('import re as _re')[0].split('src=open')[0])
src=open(T0+'steps12_22Nov/build/plot_day3.py').read().split('\nX = experiment()')[0]
g={'__file__':T0+'steps12_22Nov/build/plot_day3.py','__name__':'c'};exec(src,g)
T=1.993;H=0.1195;h=0.5;w=2*np.pi/T;k=1.555392701892035;s=np.tanh(k*h);ph0=1.5*np.pi
A2c=k*H*H/4*(3-s*s)/(4*s**3)
RUNS=[('plain fine',os.environ.get('PFREF','meshes/plain_fine')),('H_P','zonefix/hybrid/H_P'),('H_S','zonefix/hybrid/H_S')]   # PFREF: reference run
RUNS+=[tuple(a.split('=')) for a in sys.argv[1:]]          # label=path relative to T0
X=np.array((-1.25,-1.,-.25,.25,1.,1.25))
tl=[round(8+.1*i,1) for i in range(151)]
cache=os.environ.get('S6CACHE','/tmp/s6_eta.pkl')
L=pickle.load(open(cache,'rb')) if os.path.exists(cache) else {}
if 1:
    for n,p in RUNS:
        if n in L:continue
        P=T0+p;C=cells(P);P0=pts(P+'/constant/polyMesh/points')
        y0=P0[:,1].min();Q=[[i for i in c if abs(P0[i,1]-y0)<1e-9] for c in C]   # front-face quad of each cell
        Q=np.array([sorted(q,key=lambda i,c=P0[q].mean(0):np.arctan2(P0[i,2]-c[2],P0[i,0]-c[0])) for q in Q])
        near=[np.where((P0[Q,0].min(1)<x+.15)&(P0[Q,0].max(1)>x-.15))[0] for x in X]
        M=[]
        for t in tl:
            d='%g'%t
            if not os.path.exists(f'{P}/{d}/alpha.water'):M.append([np.nan]*len(X));continue
            fp=f'{P}/{d}/polyMesh/points';Pt=pts(fp) if os.path.exists(fp) else P0;q=field(f'{P}/{d}/alpha.water');v=[]
            for x,nb in zip(X,near):            # vertical line x: exact chord length through each (convex) quad
                xa=Pt[Q[nb],0];za=Pt[Q[nb],2];xb=np.roll(xa,-1,1);zb=np.roll(za,-1,1)
                cr=((xa<=x)&(x<xb))|((xb<=x)&(x<xa));tt=np.where(cr,(x-xa)/np.where(cr,xb-xa,1),np.nan);zc=za+tt*(zb-za)
                ln=np.nan_to_num(np.nanmax(zc,1)-np.nanmin(zc,1));v.append(-h+np.sum(q[nb]*ln))
                if t in (tl[0],tl[-1]):print(n,t,'gauge',x,'column length %.6f'%ln.sum(),'top',Pt[:,2].max())
            M.append(v)
        L[n]=np.array(M)
    pickle.dump(L,open(cache,'wb'))
tl=np.array(tl)
def fit(t,y):
    B=np.c_[np.ones_like(t),np.cos(w*t),np.sin(w*t),np.cos(2*w*t),np.sin(2*w*t)];c=np.linalg.lstsq(B,y,rcond=None)[0]
    return c[0],c[1]+1j*c[2],c[3]+1j*c[4]          # Re(A e^{-iwt}) = a cos + b sin with A = a + i b
tgt1=lambda x:H/2*np.exp(1j*(k*x+ph0)); tgt2=lambda x:A2c*np.exp(2j*(k*x+ph0))
T1=float(os.environ.get('S6T1',23));sel=(tl>=13-1e-9)&(tl<=T1+1e-9)
print('target: A1 %.2f mm, A2 %.2f mm, mean 0'%(H/2e-3,A2c*1e3))
print('WAVE 13-%g s:'%T1+' A1 mm (ratio to target) | phase vs target deg | A2 mm | mean mm')
G=[0,1,'box',4,5];R={}
for n,_ in RUNS:
    E=L[n];out=[]
    for gi in G:
        if gi=='box':
            r=[fit(tl[sel],E[sel,j]) for j in (2,3)];m=np.mean([a[0] for a in r]);a1=np.mean([a[1] for a in r]);a2=np.mean([a[2] for a in r])
            t1=np.mean([tgt1(x) for x in (-.25,.25)]);t2=np.mean([tgt2(x) for x in (-.25,.25)]);lab='box'
        else:
            m,a1,a2=fit(tl[sel],E[sel,gi]);t1=tgt1(X[gi]);t2=tgt2(X[gi]);lab='%+.2f'%X[gi]
        out.append((lab,m,a1,a2,t1,t2))
    R[n]=out
    print('%-10s'%n+' '.join('%s: %.1f(%.3f) %+.1f %.2f %+.2f |'%(l,abs(a1)*1e3,abs(a1)/abs(t1),np.degrees(np.angle(a1/t1)),abs(a2)*1e3,m*1e3) for l,m,a1,a2,t1,t2 in out))
print('A2 target mm at gauges:',' '.join('%.2f'%abs(o[5]*1e3) for o in R['plain fine']))
print('PER PERIOD |A1|/target from 8 s (window start):')
nP=int((T1-8)/T)
for n,_ in RUNS:
    for gi,lab in ((0,'-1.25'),('box','box'),(5,'+1.25')):
        row=[]
        for p in range(nP):
            ss=(tl>=8+p*T-1e-9)&(tl<8+(p+1)*T-1e-9)
            if gi=='box':row.append(np.mean([abs(fit(tl[ss],L[n][ss,j])[1])/abs(tgt1(X[j])) for j in (2,3)]))
            else:row.append(abs(fit(tl[ss],L[n][ss,gi])[1])/abs(tgt1(X[gi])))
        print('%-10s %-6s'%(n,lab)+' '.join('%.3f'%v for v in row))
print('   windows start',' '.join('%.2f'%(8+p*T) for p in range(nP)))
print('MEANS 13-%g s:'%T1+' surge mm (x - x0) | tension N (T - T0) | pitch deg | heave mm')
Mo={};ref={}
for n,p in RUNS:
    S=g['cfd'](T0+p);Mo[n]=S
    o={};u=np.linspace(13,T1,5001)
    for kk,(t,y) in S.items():
        v=np.interp(u,t,y);y0=y[0] if kk=='surge' else 0;m,a1,a2=fit(u,v);o[kk]=(m-y0,a1,a2)
    ref[n]=o
    print('%-10s surge %+.2f  tension %+.4f  pitch %+.3f  heave %+.2f'%(n,o['surge'][0]*1e3,o['tension'][0],o['pitch'][0],o['heave'][0]*1e3))
print('RAO: motion A1 / box-gauge A1 (mm/mm, deg/m, N/m), phase motion - box gauge (deg); and vs plain fine')
for n,_ in RUNS:
    b=R[n][2][2];row=[]
    for kk in ('surge','heave','pitch','tension'):
        a=ref[n][kk][1];sc=1 if kk in('pitch','tension') else 1
        rao=a/b;rf=ref['plain fine'][kk][1]/R['plain fine'][2][2]
        row.append('%s %.3f %+6.1f (%+.1f%%, %+.1f deg)'%(kk,abs(rao),np.degrees(np.angle(np.conj(rao))),100*(abs(rao)/abs(rf)-1),np.degrees(np.angle(np.conj(rao/rf)))))
    print('%-10s '%n+' | '.join(row))
