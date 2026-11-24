"""Checks of a run (no solver run): stray alpha above the incident surface, max air Uc, pitch harmonics vs plain fine.
usage: python3 s10.py T1 label=path [label=path ...]      (path relative to FoamWaveTank_Training_test/)
 1 per write 12.0-T1 s: max alphaI and max alpha1 over cells wholly above etaI + 2 dz (cell min z > max of the Stokes II
   incident elevation over the cell x extent + 2 x cell height); also restricted to cell centre z > +0.3.
 2 max|Uc| in air (log line 'SWENSE band', air) per 1 s and overall.
 3 pitch: fit m + sum_{n=1..3} Re(A_n e^{-inwt}) over 13-T1 s; |A_n| and phase vs plain fine over the same window."""
import numpy as np,os,re,sys
D=os.path.dirname(os.path.abspath(__file__))
T0='/gpfs/hpc-home/projects/dtu_00075/people/sitali/FoamWaveTank_Training_test/'
exec(open(os.path.join(D,'s5.py')).read().split('import re as _re')[0].split('src=open')[0])
src=open(T0+'steps12_22Nov/build/plot_day3.py').read().split('\nX = experiment()')[0]
g={'__file__':T0+'steps12_22Nov/build/plot_day3.py','__name__':'c'};exec(src,g)
T=1.993;H=0.1195;h=0.5;w=2*np.pi/T;k=1.555392701892035;s=np.tanh(k*h);ph0=1.5*np.pi
eta=lambda x,t:H/2*np.cos(k*x-w*t+ph0)+k*H*H/4*(3-s*s)/(4*s**3)*np.cos(2*(k*x-w*t+ph0))
T1=float(sys.argv[1]);RUNS=[tuple(a.split('=')) for a in sys.argv[2:]]
def fit3(t,y):
    B=np.c_[np.ones_like(t)]
    for n in (1,2,3):B=np.c_[B,np.cos(n*w*t),np.sin(n*w*t)]
    c=np.linalg.lstsq(B,y,rcond=None)[0];return [c[2*n-1]+1j*c[2*n] for n in (1,2,3)]
for n,p in RUNS:
    P=T0+p;print('=== %s (%s)'%(n,p))
    if os.environ.get('S10NOALPHA')!='1':
        C=cells(P);hasI=os.path.exists(P+'/12/alphaI');rows=[]
        for i in range(120,int(round(T1*10))+1):
            t=i/10;d='%g'%t
            if not os.path.exists(f'{P}/{d}/alpha.water'):continue
            fp=f'{P}/{d}/polyMesh/points';X=pts(fp if os.path.exists(fp) else P+'/constant/polyMesh/points')[C]
            x0,x1=X[:,:,0].min(1),X[:,:,0].max(1);z0,z1=X[:,:,2].min(1),X[:,:,2].max(1);xc=X.mean(1)
            e=eta(np.linspace(0,1,9)[:,None]*(x1-x0)+x0,t).max(0);m=np.where(z0>e+2*(z1-z0))[0];u=m[xc[m,2]>0.3]
            a1=field(f'{P}/{d}/alpha.water');aI=field(f'{P}/{d}/alphaI') if hasI else a1*np.nan
            rows.append((t,np.nanmax(aI[m]) if hasI else np.nan,xc[m][np.nanargmax(aI[m])] if hasI else (np.nan,)*3,a1[m].max(),xc[m][a1[m].argmax()],np.nanmax(aI[u]) if hasI else np.nan,a1[u].max()))
        R=np.array([(r[0],r[1],r[3],r[5],r[6]) for r in rows])
        print('stray alpha, writes %.1f-%.1f (%d): max alphaI %.3e  max alpha1 %.3e | z>0.3: alphaI %.3e alpha1 %.3e | writes <= 1e-6: alphaI %d alpha1 %d'%(R[0,0],R[-1,0],len(R),np.nanmax(R[:,1]),R[:,2].max(),np.nanmax(R[:,3]),R[:,4].max(),(R[:,1]<=1e-6).sum(),(R[:,2]<=1e-6).sum()))
        for r in rows[::10]+[rows[-1]]:print('  %.1f  alphaI %.3e (%+.3f,%+.3f)  alpha1 %.3e (%+.3f,%+.3f)  z>0.3: %.3e %.3e'%(r[0],r[1],r[2][0],r[2][2],r[3],r[4][0],r[4][2],r[5],r[6]))
    lg=P+'/log.waveIsoFoamSWENSE_dyn'
    if os.path.exists(lg):
        tt=None;A=[]
        for l in open(lg):
            if l.startswith('Time = '):tt=float(l.split()[2])
            elif l.startswith('SWENSE band:') and tt is not None:A.append((tt,float(re.search(r'air = ([-+0-9.eE]+)',l).group(1))))
        A=np.array(A);print('max air Uc: overall %.2f m/s at %.3f s (last t %.3f); per second:'%(A[:,1].max(),A[A[:,1].argmax(),0],A[-1,0]),' '.join('%d:%.1f'%(q,A[(A[:,0]>=q)&(A[:,0]<q+1),1].max()) for q in range(int(A[0,0]),int(A[-1,0])+1) if ((A[:,0]>=q)&(A[:,0]<q+1)).any()))
pf=g['cfd'](T0+'meshes/plain_fine')['pitch'];u=np.linspace(13,T1,5001);Pf=fit3(u,np.interp(u,*pf))
print('pitch harmonics 13-%g s, deg (vs plain fine: amplitude %%, phase deg):'%T1)
print('  %-10s h1 %.3f  h2 %.3f  h3 %.3f'%('plain fine',*[abs(a) for a in Pf]))
for n,p in RUNS:
    S=g['cfd'](T0+p)['pitch']
    if S[0][-1]<T1-0.05:print('  %-10s ends at %.2f s < %g'%(n,S[0][-1],T1));continue
    A=fit3(u,np.interp(u,*S))
    print('  %-10s '%n+'  '.join('h%d %.3f (%+.1f%%, %+.1f deg)'%(i+1,abs(a),100*(abs(a)/abs(b)-1),np.degrees(np.angle(np.conj(a/b)))) for i,(a,b) in enumerate(zip(A,Pf))))
