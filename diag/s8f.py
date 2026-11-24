"""Field maxima at the 0.1 s writes: max |U| water/air, max |Uc|, max Co ~ dt(|Ux-Umx|/dx + |Uz-Umz|/dz) with location."""
import numpy as np,os,sys
exec(open('s5.py').read().split('src=open')[0])
T0='/gpfs/hpc-home/projects/dtu_00075/people/sitali/FoamWaveTank_Training_test/zonefix/hybrid/'
import pickle
times=[float(a) for a in sys.argv[1].split(',')]
for run in sys.argv[2:]:
    P=T0+run;C=cells(P);D=pickle.load(open(f'/tmp/claude-7201142/s8_{run}.pkl','rb'))
    Pprev=None
    for t in times:
        d='%g'%t;Pt=pts(f'{P}/{d}/polyMesh/points') if os.path.exists(f'{P}/{d}/polyMesh/points') else pts(P+'/constant/polyMesh/points')
        X=Pt[C];xc=X.mean(1);dx=np.ptp(X[:,:,0],1);dz=np.ptp(X[:,:,2],1)
        U=field(f'{P}/{d}/U',3);Uc=field(f'{P}/{d}/Uc',3);a=field(f'{P}/{d}/alpha.water')
        dt=D['dt'][np.argmin(abs(D['t']-t))]
        co=dt*(abs(U[:,0])/dx+abs(U[:,2])/dz);w=a>0.5;mU=np.linalg.norm(U,axis=1);mUc=np.linalg.norm(Uc,axis=1)
        f=lambda m,s:'%.2f @(%+.3f,%+.3f) a=%.2f'%(m[s].max(),*xc[s][np.argmax(m[s])][[0,2]],a[s][np.argmax(m[s])])
        iw=np.where(w)[0];ia=np.where(~w)[0]
        print('%-5s %5.2f |U|w %s | |U|a %s | |Uc| %s | Co %s'%(run,t,f(mU,iw),f(mU,ia),f(mUc,np.arange(len(a))),f(co,np.arange(len(a)))))
