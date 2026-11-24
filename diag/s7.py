"""etaI (alphaI), eta (alpha.water) and etaS = eta - etaI by s6's exact vertical line cut; 1st harmonic vs Stokes II, 13-23 s.
usage: s7.py label=path(rel. to T0) ...   (default: H_S, H_P, F3b)"""
import numpy as np,os,sys
T0='/gpfs/hpc-home/projects/dtu_00075/people/sitali/FoamWaveTank_Training_test/'
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'s5.py')).read().split('src=open')[0])
T=1.993;H=0.1195;h=0.5;w=2*np.pi/T;k=1.555392701892035;ph0=1.5*np.pi
X=np.array((-1.25,-1.,1.,1.25))
RUNS=[a.split('=') for a in sys.argv[1:]] or [('H_S','zonefix/hybrid/H_S'),('H_P','zonefix/hybrid/H_P'),('F3b','zonefix/F3b')]
t0=float(os.environ.get('S7T0',13));t1=float(os.environ.get('S7T1',23))
tl=[round(t0+.1*i,1) for i in range(int(round((t1-t0)*10))+1)]
def gauges(P,fields):
    C=cells(P);P0=pts(P+'/constant/polyMesh/points');y0=P0[:,1].min()
    Q=[[i for i in c if abs(P0[i,1]-y0)<1e-9] for c in C]
    Q=np.array([sorted(q,key=lambda i,c=P0[q].mean(0):np.arctan2(P0[i,2]-c[2],P0[i,0]-c[0])) for q in Q])
    near=[np.where((P0[Q,0].min(1)<x+.15)&(P0[Q,0].max(1)>x-.15))[0] for x in X]
    out={f:[] for f in fields};tt=[]
    for t in tl:
        d='%g'%t
        if not all(os.path.exists(f'{P}/{d}/{f}') for f in fields):continue
        fp=f'{P}/{d}/polyMesh/points';Pt=pts(fp) if os.path.exists(fp) else P0;tt.append(t)
        L=[]
        for x,nb in zip(X,near):
            xa=Pt[Q[nb],0];za=Pt[Q[nb],2];xb=np.roll(xa,-1,1);zb=np.roll(za,-1,1)
            cr=((xa<=x)&(x<xb))|((xb<=x)&(x<xa));s=np.where(cr,(x-xa)/np.where(cr,xb-xa,1),np.nan);zc=za+s*(zb-za)
            ln=np.nan_to_num(np.nanmax(zc,1)-np.nanmin(zc,1));assert abs(ln.sum()-1)<1e-9;L.append(ln)
        for f in fields:
            q=field(f'{P}/{d}/{f}');out[f].append([-h+np.sum(q[nb]*ln) for nb,ln in zip(near,L)])
    return np.array(tt),{f:np.array(v) for f,v in out.items()}
def fit(t,y):
    B=np.c_[np.ones_like(t),np.cos(w*t),np.sin(w*t),np.cos(2*w*t),np.sin(2*w*t)];c=np.linalg.lstsq(B,y,rcond=None)[0];return c[0],c[1]+1j*c[2]
tgt=H/2*np.exp(1j*(k*X+ph0))
print('A1/target (phase deg vs target) at x =',X,' window %g-%g s'%(t0,t1))
res={}
for n,p in RUNS:
    fl=['alpha.water']+(['alphaI'] if os.path.exists(T0+p+'/%g/alphaI'%tl[-1]) else [])
    tt,E=gauges(T0+p,fl);S={'eta':E['alpha.water']}
    if 'alphaI' in E:S['etaI']=E['alphaI'];S['etaS']=S['eta']-S['etaI']
    res[n]={}
    for kk,V in S.items():
        A=np.array([fit(tt,V[:,j])[1] for j in range(len(X))]);res[n][kk]=A
        print('%-6s %-5s'%(n,kk)+'  '.join('%.3f(%+5.1f)'%(abs(a)/abs(g),np.degrees(np.angle(a/g))) for a,g in zip(A,tgt))+'   n=%d'%len(tt))
if 'H_P' in res and 'H_S' in res:
    print('scattered amplitude, |eta - target A1| / |target|:')
    for n in res:print('%-6s'%n,'  '.join('%.3f'%(abs(a-g)/abs(g)) for a,g in zip(res[n]['eta'],tgt)))
