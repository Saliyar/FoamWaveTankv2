"""2nd-harmonic pitch deficit, H_S3 vs H_P2, 13-23 s, from the 0.1 s writes (no solver run).
usage: python3 s11.py [S=zonefix/hybrid/H_S3] [P=zonefix/hybrid/H_P2]
 1 pitch moment about the centre of rotation (uniform/sixDoFRigidBodyMotionState), per box face
   (front x-, back x+, bottom, top; the y planes are empty in 2D) and, for S, per pressure part
   pInc, pc, Drho gh (= p - pInc - pc, checked against (rho - rhoI) g.Cf) and viscous
   (mu (U_P - U_b)_t / d at the owner cell, checked against box_forces viscous at the writes).
   Harmonics: fit m + sum_{n=1..3} Re(A_n e^{-inwt}) over the writes 13-23 s, T = 1.993 s.
 2 2nd harmonic of eta (alpha columns next to the front and back faces and at x = -1.25) and of p at
   3 body-fixed faces per face; S total and incident (alphaI, pInc) vs P total.
 3 max air |U| (alpha.water < 0.5) per second: band A = the 3 cells above the top wet cell of each column,
   band B = top 3 rows of each column; S U (total) and Uc, P U."""
import numpy as np,os,re,sys
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'s5.py')).read().split('src=open')[0])
A=dict(a.split('=') for a in sys.argv[1:]);S=T0+A.get('S','zonefix/hybrid/H_S3');PP=T0+A.get('P','zonefix/hybrid/H_P2')
T=1.993;w=2*np.pi/T;G=9.81;RW,RA,MUW,MUA=1000.,1.,1e-3,1.48e-5
TL=[round(13+0.1*i,1) for i in range(101)]
def fit(t,y):
    B=np.c_[np.ones_like(t)]
    for n in (1,2,3):B=np.c_[B,np.cos(n*w*t),np.sin(n*w*t)]
    c=np.linalg.lstsq(B,y,rcond=None)[0];return np.array([c[0]]+[c[2*n-1]+1j*c[2*n] for n in (1,2,3)])
def patch(f,name='box',nc=1):
    s=open(f).read();i=s.index('boundaryField');i=s.index('\n    '+name+'\n',i);j=s.index('\n    }',i);b=s[i:j]
    k=b.find('value');
    if k<0:return None
    m=re.search(r'value\s+uniform\s+([^;]+);',b[k:])
    if m and 'nonuniform' not in b[k:k+30]:
        v=np.fromstring(m.group(1).replace('(',' ').replace(')',' '),sep=' ');return np.tile(v,(52,1)) if nc>1 else np.full(52,v[0])
    a=b.index('(',k);e=b.index('\n)',a);v=np.fromstring(b[a+1:e].replace('(',' ').replace(')',' '),sep=' ')
    return v.reshape(-1,nc) if nc>1 else v
def mesh(P):
    q=P+'/constant/polyMesh/';s=open(q+'boundary').read()
    sf=int(re.search(r'box\s*\{[^}]*startFace\s+(\d+)',s).group(1))
    F=[list(map(int,m.split())) for m in re.findall(r'\d+\(([^)]*)\)',listbody(q+'faces')[0])]
    own=np.fromstring(listbody(q+'owner')[0],sep=' ',dtype=int)
    return [F[f] for f in range(sf,sf+52)],own[sf:sf+52],cells(P)
def cor(P,d):
    s=open(f'{P}/{d}/uniform/sixDoFRigidBodyMotionState').read()
    return np.fromstring(re.search(r'centreOfRotation\s*\(([^)]*)\)',s).group(1),sep=' ')
FACE={'front':range(32,42),'back':range(42,52),'bottom':range(0,16),'top':range(16,32)}
PTS={'front':(32,34,36),'back':(42,44,46),'bottom':(0,8,15),'top':(16,23,31)}
def col_eta(P,C,P0,Xs):
    """columns: initial x of the cell centres; eta = lowest z + sum alpha dz; per x in Xs the initial column nearest"""
    xc0=P0[C].mean(1)[:,0];key=np.round(xc0,6);ux=np.unique(key);cid=np.searchsorted(ux,key)
    return ux,cid,[int(np.argmin(abs(ux-x))) for x in Xs]
def run(P,split):
    Fb,own,C=mesh(P);P0=pts(P+'/constant/polyMesh/points')
    Xs=(-0.1-0.006,0.1+0.006,-1.25);ux,cid,kc=col_eta(P,C,P0,Xs)
    print('  eta columns at initial x',[round(ux[k],4) for k in kc])
    nCol=len(ux);order=[np.where(cid==k)[0] for k in range(nCol)]
    out={k:[] for k in ('t','M','eta','etaI','pp','ppI','UA','UB','UcA','UcB','chk','vchk')}
    for t in TL:
        d='%g'%t;fp=f'{P}/{d}/polyMesh/points';X=pts(fp) if os.path.exists(fp) else P0
        Xf=np.array([X[f] for f in Fb]);Cf=Xf.mean(1)
        Sf=np.array([0.5*np.cross(q[2]-q[0],q[3]-q[1]) for q in Xf])     # OpenFOAM orientation: out of the fluid
        R=Cf-cor(P,d);mom=lambda F:R[:,2]*F[:,0]-R[:,0]*F[:,2]
        p=patch(f'{P}/{d}/p');rho=patch(f'{P}/{d}/rho') if split else None
        Ub=patch(f'{P}/{d}/U',nc=3);U=field(f'{P}/{d}/U',3);a1=field(f'{P}/{d}/alpha.water')
        Xc=X[C].mean(1);n=Sf/np.linalg.norm(Sf,axis=1)[:,None]
        dU=U[own]-Ub;dU-=(dU*n).sum(1)[:,None]*n;dist=abs(((Xc[own]-Cf)*n).sum(1))
        mu=a1[own]*MUW+(1-a1[own])*MUA;Fv=mu[:,None]*dU/dist[:,None]*np.linalg.norm(Sf,axis=1)[:,None]
        parts={'total':mom(p[:,None]*Sf),'visc':mom(Fv)}
        if split:
            pI=patch(f'{P}/{d}/pInc');pc=patch(f'{P}/{d}/pc');aIb=patch(f'{P}/{d}/alphaI')
            rI=aIb*RW+(1-aIb)*RA;dg=(rho-rI)*(-G*Cf[:,2]);rs=p-pI-pc-dg;out['chk'].append(np.abs(rs).max())
            parts.update(pInc=mom(pI[:,None]*Sf),pc=mom(pc[:,None]*Sf),Drhogh=mom(dg[:,None]*Sf),resid=mom(rs[:,None]*Sf))
            out['ppI'].append([pI[i] for f in PTS for i in PTS[f]])
        out['M'].append({k:np.array([v[list(FACE[f])].sum() for f in FACE]) for k,v in parts.items()})
        out['pp'].append([p[i] for f in PTS for i in PTS[f]])
        # viscous check vs box_forces (moment about origin, viscous_y)
        Fp=p[:,None]*Sf;out['vchk'].append(((Cf[:,2]*Fp[:,0]-Cf[:,0]*Fp[:,2]).sum(),(Cf[:,2]*Fv[:,0]-Cf[:,0]*Fv[:,2]).sum()))
        z=X[C][:,:,2];dz=z.max(1)-z.min(1);zc=z.mean(1)
        def eta(q):return [q[cid==k].dot(dz[cid==k])+z[cid==k].min() for k in kc]
        out['eta'].append(eta(a1))
        if split:out['etaI'].append(eta(field(f'{P}/{d}/alphaI')))
        Um=np.linalg.norm(U,axis=1);Ucm=np.linalg.norm(field(f'{P}/{d}/Uc',3),axis=1) if split else None
        bA=[];bB=[]
        for k in range(nCol):
            c=order[k][np.argsort(zc[order[k]])];wet=np.where(a1[c]>=0.5)[0]
            top=wet.max() if len(wet) else -1;bA+=list(c[top+1:top+4]);bB+=list(c[-3:])
        bA=np.array(bA);bB=np.array(bB);bA=bA[a1[bA]<0.5];bB=bB[a1[bB]<0.5]
        out['UA'].append((Um[bA].max(),*Xc[bA[Um[bA].argmax()]][[0,2]]));out['UB'].append((Um[bB].max(),*Xc[bB[Um[bB].argmax()]][[0,2]]))
        if split:out['UcA'].append(Ucm[bA].max());out['UcB'].append(Ucm[bB].max())
        out['t'].append(t)
    return out
O={}
for n,P,sp in (('H_S3',S,True),('H_P2',PP,False)):
    print('=== %s (%s)'%(n,P));O[n]=run(P,sp)
t=np.array(TL)
# viscous check against box_forces at the writes
for n,P in (('H_S3',S),('H_P2',PP)):
    m=np.loadtxt(P+'/postProcessing/box_forces/0/moment.dat',comments='#');v=np.array(O[n]['vchk'])
    i=np.searchsorted(m[:,0],t-1e-9);rms=lambda x:np.sqrt((x**2).mean())
    print('%s moment about origin vs box_forces at the writes (time %.4f..%.4f): pressure rms %.3e, rms diff %.3e; viscous rms %.3e, box_forces %.3e, rms diff %.3e'%(n,m[i[0],0],m[i[-1],0],rms(m[i,5]),rms(v[:,0]-m[i,5]),rms(v[:,1]),rms(m[i,8]),rms(v[:,1]-m[i,8])))
print('H_S3: Drho gh = (rho_b - rhoI_b)(-g z_f) (hRef 0); resid = p - pInc - pc - Drho gh, max %.2f Pa'%max(O['H_S3']['chk']))
F=list(FACE)
def H(n,k):return np.array([fit(t,np.array([M[k][j] for M in O[n]['M']])) for j in range(4)])   # face x harmonic
print('\n1 PITCH MOMENT about CoR, 2nd harmonic |A2| (N m, per cell depth) and phase (deg, + = lead as in s10 conj convention)')
Hs={k:H('H_S3',k) for k in ('total','pInc','pc','Drhogh','resid','visc')};Hp={k:H('H_P2',k) for k in ('total','visc')}
ph=lambda a:np.degrees(np.angle(np.conj(a)))
tot=lambda h:h[:,2].sum()
print('   %-8s %-18s %-18s %-18s | H_S3 parts: %-16s %-16s %-16s %-16s %-16s'%('face','H_P2 p+visc','H_S3 p+visc','diff S-P','pInc','pc','Drho gh','resid','visc'))
for j,f in enumerate(F+['sum']):
    g=lambda h:(h[j,2] if j<4 else h[:,2].sum())
    a=g(Hp['total'])+g(Hp['visc']);b=g(Hs['total'])+g(Hs['visc']);c=b-a
    print('   %-8s %.3e %+6.1f   %.3e %+6.1f   %.3e %+6.1f   | '%(f,abs(a),ph(a),abs(b),ph(b),abs(c),ph(c))+'  '.join('%.3e %+6.1f'%(abs(g(Hs[k])),ph(g(Hs[k]))) for k in ('pInc','pc','Drhogh','resid','visc')))
ref=Hp['total'][:,2].sum()+Hp['visc'][:,2].sum();u=ref/abs(ref)
print('   projection of S-P on the H_P2 total 2nd-harmonic phasor (in-phase / quadrature, N m), per face:')
for j,f in enumerate(F):
    c=(Hs['total'][j,2]+Hs['visc'][j,2])-(Hp['total'][j,2]+Hp['visc'][j,2]);r=c/u;print('     %-7s %+.3e  %+.3e'%(f,r.real,r.imag))
c=(tot(Hs['total'])+tot(Hs['visc']))-ref;r=c/u;print('     %-7s %+.3e  %+.3e   (H_P2 |A2| %.3e)'%('sum',r.real,r.imag,abs(ref)))
print('   1st harmonic, sum over faces: H_P2 %.3e %+.1f  H_S3 %.3e %+.1f'%(abs(Hp['total'][:,1].sum()+Hp['visc'][:,1].sum()),ph(Hp['total'][:,1].sum()),abs(Hs['total'][:,1].sum()+Hs['visc'][:,1].sum()),ph(Hs['total'][:,1].sum())))
print('   H_S3 parts, per face, projection on the H_P2 total 2nd-harmonic phasor (in-phase, quadrature):')
for j,f in enumerate(F):print('     %-7s '%f+'  '.join('%s %+.2e %+.2e'%(k,(Hs[k][j,2]/u).real,(Hs[k][j,2]/u).imag) for k in ('pInc','pc','Drhogh','resid','visc')))
print('\n2 SECOND HARMONIC of eta (mm) and p at body-fixed faces (Pa); |A2| and phase; H_S3 total, H_S3 incident, H_P2 total; also |A1|')
E=lambda n,k:np.array(O[n][k]);lab=('front col','back col','x=-1.25')
for j,l in enumerate(lab):
    a=fit(t,E('H_S3','eta')[:,j]);b=fit(t,E('H_S3','etaI')[:,j]);c=fit(t,E('H_P2','eta')[:,j])
    print('   eta %-9s A2: S %.2f %+6.1f  SI %.2f %+6.1f  P %.2f %+6.1f   A1: S %.1f SI %.1f P %.1f'%(l,1e3*abs(a[2]),ph(a[2]),1e3*abs(b[2]),ph(b[2]),1e3*abs(c[2]),ph(c[2]),1e3*abs(a[1]),1e3*abs(b[1]),1e3*abs(c[1])))
i=0
for f in PTS:
    for q in PTS[f]:
        a=fit(t,E('H_S3','pp')[:,i]);b=fit(t,E('H_S3','ppI')[:,i]);c=fit(t,E('H_P2','pp')[:,i])
        print('   p %-6s face %2d A2: S %6.1f %+6.1f  SI %6.1f %+6.1f  P %6.1f %+6.1f   A1: S %6.1f SI %6.1f P %6.1f'%(f,q,abs(a[2]),ph(a[2]),abs(b[2]),ph(b[2]),abs(c[2]),ph(c[2]),abs(a[1]),abs(b[1]),abs(c[1])));i+=1
print('\n3 MAX AIR |U| per second (m/s), writes only; A = 3 cells above the top wet cell, B = top 3 rows')
for q in range(13,23):
    s=(t>=q-1e-9)&(t<q+1-1e-9);r=lambda n,k:np.array(O[n][k])[s]
    print('   %d-%d  A: S %5.2f (Uc %5.2f) P %5.2f   B: S %5.2f (Uc %5.2f) P %5.2f'%(q,q+1,r('H_S3','UA')[:,0].max(),r('H_S3','UcA').max(),r('H_P2','UA')[:,0].max(),r('H_S3','UB')[:,0].max(),r('H_S3','UcB').max(),r('H_P2','UB')[:,0].max()))
for n in O:
    for k in ('UA','UB'):a=np.array(O[n][k]);j=a[:,0].argmax();print('   %s %s overall %.2f at t %.1f x %+.3f z %+.3f'%(n,k,a[j,0],t[j],a[j,1],a[j,2]))
