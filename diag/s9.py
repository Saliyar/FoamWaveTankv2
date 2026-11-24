"""H_S2 restart 21.4-21.8 s (H_S2r, controlDict ucCellLog): reproduction check against H_S2 and per-cell terms of the Uc equation
in the top two rows at x 0.9-1.4 (postProcessing/ucCellLog/<first step>/ucCellLog.dat, swenseUcCellLog.H).
Growth: first step at which max over cells of |term| exceeds K x its max over the reference window 21.40-21.60 s and stays above it to the end of
a 0.01 s window."""
import numpy as np,re,sys,glob
T0='/gpfs/hpc-home/projects/dtu_00075/people/sitali/FoamWaveTank_Training_test/zonefix/hybrid/'
ta,tb=21.4,21.8;K=float(sys.argv[1]) if len(sys.argv)>1 else 2.0
num=r'([-+0-9.eE]+)'
def parse(run):
    R=[];cur=None
    for l in open(T0+run+'/log.waveIsoFoamSWENSE_dyn'):
        if l.startswith('Time = '):
            t=float(l.split()[2])
            if cur and 'dt' in cur:R.append(cur)
            cur={'t':t} if ta-1e-9<t<=tb+1e-9 else None
            if t>tb+1e-6:break
            continue
        if cur is None:continue
        if l.startswith('    Orientation:'):q=[float(v) for v in l[l.index('(')+1:l.index(')')].split()];cur['pitch']=np.degrees(np.arctan2(q[2],q[0]))
        elif l.startswith('    Angular velocity:'):cur['prate']=np.degrees(float(l[l.index('(')+1:l.index(')')].split()[1]))
        elif l.startswith('SWENSE check:'):v=re.findall(r'= '+num,l);cur['dA']=float(v[0]);cur['Uc']=float(v[1]);cur['pc']=float(v[2])
        elif l.startswith('SWENSE band:'):m=re.search(r'max\|Uc\| water = '+num+'  air = '+num,l);cur['UcW']=float(m.group(1));cur['UcA']=float(m.group(2))
        elif l.startswith('deltaT ='):cur['dt']=float(l.split()[2])
    if cur and 'dt' in cur and 'Uc' in cur:R.append(cur)
    keys=sorted(set().union(*[r.keys() for r in R]))
    return {k:np.array([r.get(k,np.nan) for r in R]) for k in keys}
def moor(run,off=0.0):
    out={}
    for i in range(1,5):
        t=[];F=[]
        with open(T0+run+f'/Mooring/lines_v2_bodyC_Line{i}.out') as f:
            next(f);next(f)
            for l in f:
                tt=float(l[:12])+off
                if tt<ta-1e-9:continue
                if tt>tb:break
                t.append(tt);F.append(float(l.split()[-1]))
        out[i]=(np.array(t),np.array(F))
    return out
A,B=parse('H_S2'),parse('H_S2r')
print('steps 21.4-21.8: H_S2 %d  H_S2r %d   first steps H_S2 %s  H_S2r %s'%(len(A['t']),len(B['t']),A['t'][:3],B['t'][:3]))
print('reproduction, max |H_S2r - H_S2| (H_S2 interpolated to H_S2r steps) and the H_S2 range in the window:')
for k in ('pitch','prate','Uc','UcW','UcA','pc','dt','dA'):
    for t1 in (21.5,21.6,21.7,21.8):
        s=B['t']<=t1+1e-9;d=np.abs(B[k][s]-np.interp(B['t'][s],A['t'],A[k]))
        sA=A['t']<=t1+1e-9
        print('  %-6s to %.1f: %.3g  (H_S2 %.4g..%.4g)'%(k,t1,np.nanmax(d),np.nanmin(A[k][sA]),np.nanmax(A[k][sA])),end='' if t1<21.8 else '\n')
try:
    MA,MB=moor('H_S2'),moor('H_S2r',21.4)  # MoorDyn restarts its clock at 0
    for i in range(1,5):
        for t1 in (21.5,21.6,21.7,21.8):
            s=MB[i][0]<=t1;d=np.abs(MB[i][1][s]-np.interp(MB[i][0][s],MA[i][0],MA[i][1]))
            print('  T%d    to %.1f: %.3g N'%(i,t1,d.max()),end='' if t1<21.8 else '\n')
except Exception as e:print('mooring:',e)
# per-cell terms
LOG=glob.glob(T0+'H_S2r/postProcessing/ucCellLog/*/ucCellLog.dat')[0]
L=[l.split() for l in open(LOG) if l[0]=='C']
M=np.array([[float(v) for v in l[1:]] for l in L])
cols='t dt cell x z a1 aI rho rhoI lam w Ucx Ucz pc'.split()
names='ddt conv visc inc clo face gpc grav R'.split()
for n in names:cols+=[n+'x',n+'z']
ix={c:i for i,c in enumerate(cols)}
ts=np.unique(M[:,0]);cells=np.unique(M[:,2]).astype(int);nc=len(cells)
M=M[np.lexsort((M[:,2],M[:,0]))];M=M[:len(ts)*nc].reshape(len(ts),nc,-1)
g=lambda n:np.hypot(M[:,:,ix[n+'x']],M[:,:,ix[n+'z']])
print('cells',cells.tolist(),' steps',len(ts))
ref=ts<=21.6+1e-9
print('first growth (max over cells > %g x its 21.4-21.6 max, sustained 0.01 s):'%K)
res=[]
for n in names+['Uc','pc','rho']:
    G=g(n) if n in names else (np.hypot(M[:,:,ix['Ucx']],M[:,:,ix['Ucz']]) if n=='Uc' else np.abs(M[:,:,ix[n]]))
    m=G.max(1);r=m[ref].max();ex=m>K*r;first=None
    for j in np.where(ex&~ref)[0]:
        w=(ts>=ts[j])&(ts<=ts[j]+0.01)
        if ex[w].all():first=j;break
    res.append((ts[first] if first is not None else 99,n,r,m.max(),cells[G[first].argmax()] if first is not None else -1))
res.sort()
for f,n,r,mx,c in res:print('  %-5s %s  ref max %.3g  run max %.3g  cell %s'%(n,'%.4f'%f if f<99 else 'none  ',r,mx,c if c>=0 else '-'))
# time series at the cell with the largest final |Uc|
UcM=np.hypot(M[:,:,ix['Ucx']],M[:,:,ix['Ucz']]);kc=int(UcM[-1].argmax())
print('cell %d (x %.3f z %.3f): per step, z components (N/m3)'%(cells[kc],M[0,kc,ix['x']],M[0,kc,ix['z']]))
print('    t        dt    a1   aI  rho/rhoI lam   w     Ucx    Ucz     pc   '+' '.join('%8s'%n for n in names))
step=max(1,len(ts)//40)
for j in list(range(0,len(ts),step))+[len(ts)-1]:
    r=M[j,kc];print('%8.4f %.2e %.2f %.2f %6.2f %.2f %.2f %6.2f %6.2f %7.1f '%(r[0],r[1],r[ix['a1']],r[ix['aI']],r[ix['rho']]/r[ix['rhoI']],r[ix['lam']],r[ix['w']],r[ix['Ucx']],r[ix['Ucz']],r[ix['pc']])+' '.join('%8.1f'%r[ix[n+'z']] for n in names))
# atmosphere boundary
Bl=np.array([[float(v) for v in l.split()[1:]] for l in open(LOG) if l[0]=='B'])
print('atmosphere faces x 0.9-1.4: max|pc_b| %.3g, per 0.05 s: max phiComp/|Sf| (out +), min, max|snGrad pc|, max|Uc_b|'%np.abs(Bl[:,4]).max())
for t0 in np.arange(21.4,21.8,0.05):
    s=(Bl[:,0]>t0)&(Bl[:,0]<=t0+0.05)
    if s.any():print('  %.2f-%.2f  %+.3f %+.3f  %.3g  %.3g'%(t0,t0+0.05,Bl[s,7].max(),Bl[s,7].min(),np.abs(Bl[s,8]).max(),np.hypot(Bl[s,5],Bl[s,6]).max()))
np.save('/tmp/claude-7201142/s9_M.npy',M)
# Second test: growth against the first 0.01 s (21.40-21.41), so a term already rising in 21.4-21.6 is not hidden.
# Quantities: each term (max over cells of |.|), Drho = rho - rhoI, alphaI, |Uc|, and the atmosphere inflow max(-phiComp/|Sf|).
print('growth against 21.40-21.41 (first step with max over cells >= 2x and >= 1.5x, sustained 0.01 s); terms with 21.40 max < 1 N/m3 skipped:')
r0=ts<=21.41+1e-9
Q={n:g(n) for n in names}
Q['Drho']=np.abs(M[:,:,ix['rho']]-M[:,:,ix['rhoI']]);Q['alphaI']=M[:,:,ix['aI']];Q['Uc']=np.hypot(M[:,:,ix['Ucx']],M[:,:,ix['Ucz']])
tb_=np.unique(Bl[:,0]);inf=np.array([-Bl[Bl[:,0]==t,7].min() for t in tb_])
def first(tt,m,f):
    r=m[tt<=21.41+1e-9].max();ex=m>=f*r
    for j in np.where(ex&(tt>21.41))[0]:
        if ex[(tt>=tt[j])&(tt<=tt[j]+0.01)].all():return tt[j],r
    return 99,r
out=[]
for n,G in Q.items():
    m=G.max(1)
    if n in names and m[r0].max()<1:out.append((99,99,n,m[r0].max(),m.max(),'-'));continue
    a,r=first(ts,m,1.5);b,_=first(ts,m,2.0);out.append((a,b,n,r,m.max(),cells[G[ts==a].argmax()] if a<99 else (cells[G[ts==b].argmax()] if b<99 else '-')))
a,r=first(tb_,inf,1.5);b,_=first(tb_,inf,2.0);out.append((a,b,'inflow',r,inf.max(),'atm'))
out.sort()
for a,b,n,r,mx,c in out:print('  %-7s x1.5 at %s  x2 at %s  (21.40 %.3g, max %.3g, cell %s)'%(n,'%.4f'%a if a<99 else ' none ','%.4f'%b if b<99 else ' none ',r,mx,c))
