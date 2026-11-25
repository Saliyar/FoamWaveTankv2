"""Gate c of the height ramp (no solver run): coarse empty tank, streamFunction N = 10, rampMode height.
usage: python3 s13.py <case> [t0 t1]   (default window 12-15 s, after the ramp)
eta (alpha.water, waveGauges) and etaI (alphaI, waveGaugesI) at x = -1, +1 from interfaceHeight hB - h; fit
m + sum_{n=1..3} Re(A_n e^{-inwt}) over the window; compare A1, A2 (amplitude %, phase deg) and m with the stream
function theory of <zonefix>/hr/sf/fenton4Foam (as streamFunctionProperties rescales it, phi = 1.5 pi)."""
import numpy as np,re,sys
case=sys.argv[1];t0,t1=(float(sys.argv[2]),float(sys.argv[3])) if len(sys.argv)>3 else (12.0,15.0)
g,h=9.81,0.5;phi=1.5*np.pi
s=open('/gpfs/hpc-home/projects/dtu_00075/people/sitali/FoamWaveTank_Training_test/zonefix/hr/sf/fenton4Foam').read()
num=lambda key:float(re.search(key+r'\s+'+key+r'\s+\[[^\]]*\]\s+([-+0-9.eE]+);',s).group(1))
lst=lambda key:np.array([float(v) for v in re.search(key+r'[^(]*\(([^)]*)\)',s).group(1).split()])
k=num('depth')/h;T=num('period')/np.sqrt(g*k);w=2*np.pi/T;A=lst('aCoeffs')/k
etaT=lambda x,t:sum(A[j-1]*np.cos(j*(k*x-w*t+phi)) for j in range(1,11))
def fit(t,y):
    B=np.c_[np.ones_like(t)]
    for n in (1,2,3):B=np.c_[B,np.cos(n*w*t),np.sin(n*w*t)]
    c=np.linalg.lstsq(B,y,rcond=None)[0];return c[0],[c[2*n-1]+1j*c[2*n] for n in (1,2,3)]
def gauges(f):
    hdr=[l for l in open(f) if l.startswith('# Location ') and ':' in l]
    X=[float(re.search(r'\(([-+0-9.eE]+)',l).group(1)) for l in hdr]
    d=np.loadtxt(f,comments='#');return X,d
print('window %.1f-%.1f s, theory k %.6f T %.6f'%(t0,t1,k,T))
print('%-6s %5s %9s %9s %9s %9s %9s'%('field','x','A1 %','A1 deg','A2 %','A2 deg','mean mm'))
for name,f in (('eta',case+'/postProcessing/waveGauges/0/height.dat'),('etaI',case+'/postProcessing/waveGaugesI/0/height.dat')):
    X,d=gauges(f);s_=(d[:,0]>=t0-1e-9)&(d[:,0]<=t1+1e-9);t=d[s_,0]
    for x in (-1.0,1.0):
        j=X.index(x);y=d[s_,1+2*j]-h
        m,F=fit(t,y);mT,FT=fit(t,etaT(x,t))
        r=[(100*(abs(F[n])/abs(FT[n])-1),np.degrees(np.angle(np.conj(F[n]/FT[n])))) for n in (0,1)]
        print('%-6s %+5.1f %+9.2f %+9.2f %+9.2f %+9.2f %+9.2f'%(name,x,r[0][0],r[0][1],r[1][0],r[1][1],1e3*(m-mT)))
