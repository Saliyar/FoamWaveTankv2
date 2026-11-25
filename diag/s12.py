"""Offline free-surface residuals of the regular wave theories for H 0.1195, T 1.993, h 0.5 (no solver run).
Formulas as coded in src/waves2Foam/waveTheories/regular (ramp factor 1, seaLevel 0, wave along +x):
  stokesSecond   (k linear; return current -gH^2/(8ch) in u)
  stokesFifth    (k from stokesFifthProperties::eval with stokesDrift Q = 0)
  streamFunction N = 10 (A, B, uBar, k, T from fenton4Foam, Stokes velocity 0, read from <sf10>/constant/fenton4Foam
                 as streamFunctionProperties::set rescales them)
  StokesIIIncident, fluxConsistentU off and on (the H_S3 incident field)
On one wavelength at t = 0, theta = kx - wt, nTh points, spectral x derivatives, eta_t = -c eta_x:
  mass-flux residual  R_m = eta_t + dQ/dx,  Q = int_{-h}^{eta} u dz (Gauss-Legendre)
  kinematic residual  R_k = eta_t + u_s eta_x - w_s at z = eta
  both scaled by w H/2 (max |.|); dynamic: range over theta of 1/2((u-c)^2 + w^2) + g eta at the surface, scaled by g H.
usage: python3 s12.py <dir holding constant/fenton4Foam for N = 10>"""
import numpy as np,re,sys
H,T,h,g=0.1195,1.993,0.5,9.81;w0=2*np.pi/T;n=2048
th=2*np.pi*np.arange(n)/n
def klin(om):
    k=om*om/g
    for _ in range(200):k=om*om/(g*np.tanh(k*h))
    return k
def ddx(f,k):F=np.fft.rfft(f);return np.fft.irfft(1j*np.arange(len(F))*k*F,n)
xg,wg=np.polynomial.legendre.leggauss(64)
def residuals(name,k,om,eta,uw):
    c=om/k;e=eta(th);ex=ddx(e,k);et=-c*ex
    us,ws=uw(th,e)
    Z=(e[None,:]+h)*(xg[:,None]+1)/2-h;u,_=uw(th[None,:],Z);Q=((e+h)/2)*(wg[:,None]*u).sum(0)
    Rm=et+ddx(Q,k);Rk=et+us*ex-ws;B=0.5*((us-c)**2+ws**2)+g*e
    s=om*H/2
    return name,k,om,np.abs(Rm).max()/s,np.abs(Rk).max()/s,(B.max()-B.min())/(g*H),Q.mean(),e.max()-e.min()
R=[]
# stokesSecond (stokesSecond.C:107-222), arg = wt - kx + phi -> use -theta
k=klin(w0);c=w0/k
R.append(residuals('stokesSecond',k,w0,
    lambda t:H/2*np.cos(t)+k*H*H/16*(3/np.tanh(k*h)**3-1/np.tanh(k*h))*np.cos(2*t),
    lambda t,z:(np.pi*H/T*np.cosh(k*(z+h))/np.sinh(k*h)*np.cos(t)+3/16*c*(k*H)**2*np.cosh(2*k*(z+h))/np.sinh(k*h)**4*np.cos(2*t)-g*H*H/(8*c*h)+0*z,
                np.pi*H/T*np.sinh(k*(z+h))/np.sinh(k*h)*np.sin(t)+3/16*c*(k*H)**2*np.sinh(2*k*(z+h))/np.sinh(k*h)**4*np.sin(2*t))))
# StokesIIIncident (StokesIIIncident.C:118-163, 185-263), fU off / on
for fUon in (False,True):
    k=klin(w0);s=np.tanh(k*h);fU=1.0
    if fUon:
        t_=2*np.pi*(np.arange(3600)+0.5)/3600;e_=H/2*np.cos(t_)+k*H*H/4*(3-s*s)/(4*s**3)*np.cos(2*t_)
        A_=H/2*w0*np.cos(t_)*np.sinh(k*(h+e_))/(k*np.sinh(k*h));B_=3/4*H*H/4*w0*k*np.sinh(2*k*(h+e_))/(2*k*np.sinh(k*h)**4)*np.cos(2*t_)
        A1=(2*A_*np.cos(t_)).mean();B1=(2*B_*np.cos(t_)).mean();fU=(w0/k*H/2-B1)/A1
    R.append(residuals('StokesIIIncident fU %.4f'%fU,k,w0,
        lambda t:H/2*np.cos(t)+k*H*H/4*(3-s*s)/(4*s**3)*np.cos(2*t),
        lambda t,z,fU=fU:(fU*H/2*w0*np.cos(t)*np.cosh(k*(z+h))/np.sinh(k*h)+3/4*H*H/4*w0*k*np.cosh(2*k*(z+h))/np.sinh(k*h)**4*np.cos(2*t),
                         fU*H/2*w0*np.sin(t)*np.sinh(k*(z+h))/np.sinh(k*h)+3/4*H*H/4*w0*k*np.sinh(2*k*(z+h))/np.sinh(k*h)**4*np.sin(2*t))))
# stokesFifth (stokesFifthProperties.C:48-109, stokesFifth.C:73-329), Q = 0
def S5(k):
    S=1/np.cosh(2*k*h);sh=np.sinh(k*h);th_=np.tanh(k*h)
    C=dict(A11=1/sh,A22=3*S**2/(2*(1-S)**2),A31=(-4-20*S+10*S**2-13*S**3)/(8*sh*(1-S)**3),A33=(-2*S**2+11*S**3)/(8*sh*(1-S)**3),
      A42=(12*S-14*S**2-264*S**3-45*S**4-13*S**5)/(24*(1-S)**5),A44=(10*S**3-174*S**4+291*S**5+278*S**6)/(48*(3+2*S)*(1-S)**5),
      A51=(-1184+32*S+13232*S**2+21712*S**3+20940*S**4+12554*S**5-500*S**6-3341*S**7-670*S**8)/(64*sh*(3+2*S)*(4+S)*(1-S)**6),
      A53=(4*S+105*S**2+198*S**3-1376*S**4-1302*S**5-117*S**6+58*S**7)/(32*sh*(3+2*S)*(1-S)**6),
      A55=(-6*S**3+272*S**4-1552*S**5+852*S**6+2029*S**7+430*S**8)/(64*sh*(3+2*S)*(4+S)*(1-S)**6),
      B22=(1/th_)*(1+2*S)/(2*(1-S)),B31=-3*(1+3*S+3*S**2+2*S**3)/(8*(1-S)**3),
      B42=(1/th_)*(6-26*S-182*S**2-204*S**3-25*S**4+26*S**5)/(6*(3+2*S)*(1-S)**4),
      B44=(1/th_)*(24+92*S+122*S**2+66*S**3+67*S**4+34*S**5)/(24*(3+2*S)*(1-S)**4),
      B53=9*(132+17*S-2216*S**2-5897*S**3-6292*S**4-2687*S**5+194*S**6+467*S**7+82*S**8)/(128*(3+2*S)*(4+S)*(1-S)**6),
      B55=5*(300+1579*S+3176*S**2+2949*S**3+1188*S**4+675*S**5+1326*S**6+827*S**7+130*S**8)/(384*(3+2*S)*(4+S)*(1-S)**6),
      C0=np.sqrt(th_),C2=np.sqrt(th_)*(2+7*S**2)/(4*(1-S)**2),C4=np.sqrt(th_)*(4+32*S-116*S**2-400*S**3-71*S**4+146*S**5)/(32*(1-S)**5),
      D2=-np.sqrt(1/th_)/2,D4=np.sqrt(1/th_)*(2+4*S+S**2+2*S**3)/(8*(1-S)**3))
    return C
def ev5(k,Q=0.0):
    C=S5(k);return np.sqrt(k/g)*Q-2*np.pi/(T*np.sqrt(g*k))+C['C0']+(k*H/2)**2*(C['C2']+C['D2']/(k*h))+(k*H/2)**4*(C['C4']+C['D4']/(k*h))
lo,hi=1e-7,max(4*np.pi/(T*np.sqrt(g*h)),2*np.pi/T**2)
for _ in range(200):
    m=0.5*(lo+hi)
    if np.sign(ev5(lo))==np.sign(ev5(m)):lo=m
    else:hi=m
k=0.5*(lo+hi);C=S5(k);ep=k*H/2;c=w0/k;co=C['C0']*np.sqrt(g/k**3);ub=(C['C0']+ep**2*C['C2']+ep**4*C['C4'])/np.sqrt(k/g)
def eta5(t):
    return (ep*np.cos(t)+ep**2*C['B22']*np.cos(2*t)+ep**3*C['B31']*(np.cos(t)-np.cos(3*t))+ep**4*(C['B42']*np.cos(2*t)+C['B44']*np.cos(4*t))
            +ep**5*(-(C['B53']+C['B55'])*np.cos(t)+C['B53']*np.cos(3*t)+C['B55']*np.cos(5*t)))/k
TERMS=[(1,1,'A11'),(2,2,'A22'),(3,1,'A31'),(3,3,'A33'),(4,2,'A42'),(4,4,'A44'),(5,1,'A51'),(5,3,'A53'),(5,5,'A55')]
def uw5(t,z):
    Zb=z+h;u=c-ub+0*z;v=0*z
    for p,j,a in TERMS:u=u+ep**p*j*co*k*C[a]*np.cos(j*t)*np.cosh(j*k*Zb);v=v+ep**p*j*co*k*C[a]*np.sin(j*t)*np.sinh(j*k*Zb)
    return u,v
R.append(residuals('stokesFifth Q=0',k,w0,eta5,uw5))
# streamFunction N = 10 (streamFunctionProperties.C:145-165, streamFunction.C:112-223)
s=open(sys.argv[1]+'/constant/fenton4Foam').read()
num=lambda key:float(re.search(key+r'\s+'+key+r'\s+\[[^\]]*\]\s+([-+0-9.eE]+);',s).group(1))
lst=lambda key:np.array([float(v) for v in re.search(key+r'[^(]*\(([^)]*)\)',s).group(1).split()])
k=num('depth')/h;Tsf=num('period')/np.sqrt(g*k);om=2*np.pi/Tsf;ubar=num('uBar')*np.sqrt(g/k)
A=lst('aCoeffs')/k;Bc=lst('bCoeffs')*np.arange(1,11)*np.sqrt(g/k);J=np.arange(1,11)
def etaS(t):return sum(A[j-1]*np.cos(j*t) for j in J)
def uwS(t,z):
    u=om/k-ubar+0*z;v=0*z
    for j in J:u=u+Bc[j-1]*np.cosh(j*k*(z+h))/np.cosh(j*k*h)*np.cos(j*t);v=v+Bc[j-1]*np.sinh(j*k*(z+h))/np.cosh(j*k*h)*np.sin(j*t)
    return u,v
R.append(residuals('streamFunction N=10 (T %.5f)'%Tsf,k,om,etaS,uwS))
print('H %g T %g h %g; residuals scaled by w H/2 (R_m, R_k, max over the wavelength) and g H (dynamic range)'%(H,T,h))
print('%-34s %9s %9s %9s %9s %9s %10s %8s'%('theory','k','omega','R_m','R_k','dyn','mean Q','H eta'))
for r in R:print('%-34s %9.5f %9.5f %9.2e %9.2e %9.2e %+10.2e %8.5f'%r)
