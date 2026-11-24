"""H_S2 blow-up (22.1 s): per-step log quantities, 4 MoorDyn lines, field maxima at 0.1 s writes; H_S2 against H_S.
Departure: first time |H_S2 - H_S| exceeds 5x its rms over the baseline window 20.0-21.5 s and stays above it for 0.02 s."""
import numpy as np,re,os,sys,pickle
T0='/gpfs/hpc-home/projects/dtu_00075/people/sitali/FoamWaveTank_Training_test/zonefix/hybrid/'
ta,tb=20.0,22.7
num=r'([-+0-9.eE]+)'
def parse(run):
    c=f'/tmp/claude-7201142/s8_{run}.pkl'
    if os.path.exists(c):return pickle.load(open(c,'rb'))
    R=[];cur=None
    for l in open(T0+run+'/log.waveIsoFoamSWENSE_dyn'):
        if l.startswith('Time = '):
            t=float(l.split()[2])
            if cur and 'dt' in cur:R.append(cur)
            cur={'t':t} if ta-0.01<=t<=tb else None
            if t>tb:break
            continue
        if cur is None:continue
        if l.startswith('    Orientation:'):q=[float(v) for v in l[l.index('(')+1:l.index(')')].split()];cur['pitch']=np.degrees(np.arctan2(q[2],q[0]))
        elif l.startswith('    Angular velocity:'):cur['prate']=np.degrees(float(l[l.index('(')+1:l.index(')')].split()[1]))
        elif l.startswith('SWENSE check:'):
            v=re.findall(r'= '+num,l);cur['dA']=float(v[0]);cur['Uc']=float(v[1])
            m=re.search(r'alpha1 \['+num+', '+num+r'\]\s+alphaI \['+num+', '+num,l);cur.update(a1min=float(m.group(1)),a1max=float(m.group(2)),aImin=float(m.group(3)),aImax=float(m.group(4)))
        elif l.startswith('SWENSE band:'):
            m=re.search(r'max\|Uc\| water = '+num+'  air = '+num,l);cur['UcW']=float(m.group(1));cur['UcA']=float(m.group(2))
        elif l.startswith('Courant Number mean'):cur['Co']=float(l.split()[-1])
        elif l.startswith('Interface Courant'):cur['CoI']=float(l.split()[-1])
        elif l.startswith('deltaT ='):cur['dt']=float(l.split()[2])
    keys=sorted(set().union(*[r.keys() for r in R]))
    D={k:np.array([r.get(k,np.nan) for r in R]) for k in keys};pickle.dump(D,open(c,'wb'));return D
def moor(run):
    c=f'/tmp/claude-7201142/s8m_{run}.pkl'
    if os.path.exists(c):return pickle.load(open(c,'rb'))
    out={}
    for i in range(1,5):
        t=[];F=[];mn=[]
        with open(T0+run+f'/Mooring/lines_v2_bodyC_Line{i}.out') as f:
            next(f);next(f)
            for l in f:
                tt=float(l[:12])
                if tt<ta-0.01:continue
                if tt>tb:break
                v=np.array(l.split()[-40:],float);t.append(tt);F.append(v[-1]);mn.append(v.min())
        out[i]=(np.array(t),np.array(F),np.array(mn))
    pickle.dump(out,open(c,'wb'));return out
A,B=parse('H_S'),parse('H_S2');MA,MB=moor('H_S'),moor('H_S2')
S={k:(B['t'],B[k],A['t'],A[k]) for k in ('pitch','prate','Uc','UcW','UcA','Co','CoI','dt','dA','a1min','a1max','aImin','aImax')}
for i in range(1,5):
    S['T%d'%i]=(MB[i][0],MB[i][1],MA[i][0],MA[i][1]);S['Tmin%d'%i]=(MB[i][0],MB[i][2],MA[i][0],MA[i][2])
print('steps H_S %d, H_S2 %d (20.0-22.7 s)'%(len(A['t']),len(B['t'])))
print('slack (min segment tension <= 0): H_S2',{i:(MB[i][0][MB[i][2]<=0][:1].tolist() or '-') for i in MB},' H_S',{i:(MA[i][0][MA[i][2]<=0][:1].tolist() or '-') for i in MA})
dep=[]
for k,(tB,yB,tA,yA) in S.items():
    ok=np.isfinite(yB);tB,yB=tB[ok],yB[ok];d=yB-np.interp(tB,tA,yA)
    base=(tB>=20)&(tB<=21.5);r=np.sqrt(np.mean(d[base]**2))+1e-300;ex=np.abs(d)>5*r
    first=None
    for j in np.where(ex&(tB>21.5))[0]:
        w=(tB>=tB[j])&(tB<=tB[j]+0.02)
        if ex[w].all() and w.sum()>2:first=tB[j];break
    dep.append((first if first is not None else 99,k,r))
dep.sort()
print('first departure |H_S2 - H_S| > 5 x baseline rms, sustained 0.02 s:')
for f,k,r in dep:print('  %-6s %s  (baseline rms %.3g)'%(k,'%.4f'%f if f<99 else 'none',r))
print('time   pitch(S2/S)    rate deg/s(S2/S)  Uc max(S2/S)  UcW band  UcA band   Co(S2/S)    dt(S2/S)        a1min a1max aImin aImax  T1..T4 fairlead N (S2)  [H_S]')
for t in np.arange(20.0,22.701,0.1 if True else 0.05):
    j=np.argmin(abs(B['t']-t));i=np.argmin(abs(A['t']-t))
    if abs(B['t'][j]-t)>0.01:continue
    T=[np.interp(t,MB[q][0],MB[q][1]) for q in range(1,5)];TA=[np.interp(t,MA[q][0],MA[q][1]) for q in range(1,5)]
    print('%5.2f %+6.2f/%+6.2f %+7.1f/%+7.1f %6.2f/%6.2f %5.3f/%5.3f %5.1f/%5.1f %.3f/%.3f %.2e/%.2e %g %g %g %g  %s [%s]'%(t,B['pitch'][j],A['pitch'][i],B['prate'][j],A['prate'][i],B['Uc'][j],A['Uc'][i],B['UcW'][j],A['UcW'][i],B['UcA'][j],A['UcA'][i],B['Co'][j],A['Co'][i],B['dt'][j],A['dt'][i],B['a1min'][j],B['a1max'][j],B['aImin'][j],B['aImax'][j],' '.join('%.3f'%v for v in T),' '.join('%.3f'%v for v in TA)))
