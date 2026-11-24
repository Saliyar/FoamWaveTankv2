import numpy as np,os,sys
exec(open('s7.py').read().split("tgt=H/2")[0].replace("X=np.array((-1.25,-1.,1.,1.25))","X=np.array(%s)"%os.environ['XS']))
P=T0+sys.argv[1];tt,E=gauges(P,['alphaI'])
tgt=H/2*np.exp(1j*(k*X+ph0))
for j,x in enumerate(X):
    a=fit(tt,E['alphaI'][:,j]);print('x %+5.2f  etaI A1/tgt %.3f  phase %+5.1f  mean %+.2f mm'%(x,abs(a[1])/abs(tgt[j]),np.degrees(np.angle(a[1]/tgt[j])),a[0]*1e3))
