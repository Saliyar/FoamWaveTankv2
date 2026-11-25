"""Theory matrix of SWENSE (no solver run): empty coarse tank, one case per theory x wave (zonefix/th/<case>).
usage: python3 s14.py [t0 t1] case [case ...]        (window default 9-15 s; cases relative to zonefix/th/)
Per case, from log.waveIsoFoamSWENSE_dyn and postProcessing:
  start-up residual line (kinematic, mass flux of omega H/2); FatalError text if the run stopped
  first 200 steps and whole run: max of 'SWENSE check: max|alpha1 - alphaI|' and of max|Uc|
  gauges x = -1, +1: eta (waveGauges, alpha.water) and etaI (waveGaugesI, alphaI) against the theory eta as SWENSE uses
  it (postProcessing/incidentEta/eta.dat, written at start-up from the forcing, ramp included). Fit m + sum_{n=1..3}
  Re(A_n e^{-inwt}) over the window at T 1.993 s; A1 % and deg, A2 mm (gauge / theory) and % when theory A2 > 0.5 mm,
  mean mm. Irregular (name ir*): rms(etaI - theory)/rms(theory) %, correlation, mean mm instead."""
import numpy as np,re,sys,os
Z='/gpfs/hpc-home/projects/dtu_00075/people/sitali/FoamWaveTank_Training_test/zonefix/th/'
a=sys.argv[1:];t0,t1=9.0,15.0
if len(a)>1 and re.match(r'^[0-9.]+$',a[0]):t0,t1=float(a[0]),float(a[1]);a=a[2:]
w=2*np.pi/1.993;h=0.5
def fit(t,y):
    B=np.c_[np.ones_like(t)]
    for n in (1,2,3):B=np.c_[B,np.cos(n*w*t),np.sin(n*w*t)]
    c=np.linalg.lstsq(B,y,rcond=None)[0];return c[0],[c[2*n-1]+1j*c[2*n] for n in (1,2,3)]
def gauges(f):
    X=[float(re.search(r'\(([-+0-9.eE]+)',l).group(1)) for l in open(f) if l.startswith('# Location ') and ':' in l]
    return X,np.loadtxt(f,comments='#')
for c in a:
    P=Z+c;L=open(P+'/log.waveIsoFoamSWENSE_dyn',errors='replace').read() if os.path.exists(P+'/log.waveIsoFoamSWENSE_dyn') else ''
    print('=== %s'%c)
    r=re.search(r'SWENSE incident residual \(([^)]*)\): kinematic (\S+), mass flux (\S+)',L)
    print('  residual: '+('%s  kinematic %.3g  mass flux %.3g  -> default 0.005 %s'%(r.group(1),float(r.group(2)),float(r.group(3)),'STOPS' if max(float(r.group(2)),float(r.group(3)).__float__())>0.005 else 'passes') if r else 'none printed'))
    f=re.search(r'FOAM FATAL (?:IO )?ERROR.*?\n\n(.*?)\n\n',L,re.S)
    if f:print('  FATAL: '+' '.join(f.group(1).split())[:400]);continue
    ck=[(float(m.group(1)),float(m.group(2))) for m in re.finditer(r'SWENSE check: max\|alpha1 - alphaI\| = (\S+)\s+max\|Uc\| = (\S+)',L)]
    if ck:
        A=np.array(ck);print('  steps %d: first 200 max|alpha1-alphaI| %.3g max|Uc| %.3g | all max|alpha1-alphaI| %.3g max|Uc| %.3g'%(len(A),A[:200,0].max(),A[:200,1].max(),A[:,0].max(),A[:,1].max()))
    tf=P+'/postProcessing/incidentEta/eta.dat'
    if not os.path.exists(tf):print('  no theory eta file');continue
    th=np.loadtxt(tf,comments='#')
    for name,g in (('eta','waveGauges'),('etaI','waveGaugesI')):
        fg=P+'/postProcessing/%s/0/height.dat'%g
        if not os.path.exists(fg):continue
        X,d=gauges(fg);s=(d[:,0]>=t0-1e-9)&(d[:,0]<=t1+1e-9);t=d[s,0]
        if len(t)<10:print('  %s: window not reached (last t %.2f)'%(name,d[-1,0]));continue
        out=[]
        for k,x in enumerate((-1.0,1.0)):
            y=d[s,1+2*X.index(x)]-h;yt=np.interp(t,th[:,0],th[:,1+k])
            if c.startswith('ir'):
                e=y-yt;out.append('x %+.0f: rms err %.1f %%, corr %.4f, mean %+.2f mm'%(x,100*np.sqrt(np.mean((e-e.mean())**2))/np.sqrt(np.mean((yt-yt.mean())**2)),np.corrcoef(y,yt)[0,1],1e3*(y.mean()-yt.mean())))
                continue
            m,F=fit(t,y);mT,FT=fit(t,yt)
            a1=100*(abs(F[0])/abs(FT[0])-1);p1=np.degrees(np.angle(np.conj(F[0]/FT[0])))
            a2='A2 %.2f/%.2f mm'%(1e3*abs(F[1]),1e3*abs(FT[1]))
            if abs(FT[1])>5e-4:a2+=' (%+.1f %%, %+.1f deg)'%(100*(abs(F[1])/abs(FT[1])-1),np.degrees(np.angle(np.conj(F[1]/FT[1]))))
            out.append('x %+.0f: A1 %+.2f %% %+.2f deg, %s, mean %+.2f mm'%(x,a1,p1,a2,1e3*(m-mT)))
        print('  %-4s '%name+' | '.join(out))
