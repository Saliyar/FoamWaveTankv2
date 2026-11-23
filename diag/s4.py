import sys
import numpy as np,re,glob,os
T0='/gpfs/hpc-home/projects/dtu_00075/people/sitali/FoamWaveTank_Training_test/'
def body(f):
    s=open(f).read();i=s.find('internalField');return s,i
def field(f,nc=1):
    s,i=body(f);h=s[i:min(i+120,len(s))]
    if 'nonuniform' not in h.split('\n')[0]:return None
    a=s.index('(',i);b=s.index('\n)',a)
    v=np.fromstring(s[a+1:b].replace('(',' ').replace(')',' '),sep=' ')
    return v.reshape(-1,nc) if nc>1 else v
def listbody(f):
    s=open(f).read();m=re.search(r'\n(\d+)\s*\n\(',s);return s[m.end():s.index('\n)',m.end())],int(m.group(1))
def pts(f):
    s,n=listbody(f);return np.fromstring(s.replace('(',' ').replace(')',' '),sep=' ').reshape(n,3)
def cells(p):
    q=p+'/constant/polyMesh/';F=[list(map(int,m.split())) for m in re.findall(r'\d+\(([^)]*)\)',listbody(q+'faces')[0])]
    C={}
    for n in ('owner','neighbour'):
        for f,c in enumerate(np.fromstring(listbody(q+n)[0],sep=' ',dtype=int)):C.setdefault(c,set()).update(F[f])
    return np.array([sorted(C[i]) for i in range(len(C))])
src=open(T0+'steps12_22Nov/build/plot_day3.py').read().split('\nRUNS = [')[0]
g={'__file__':T0+'steps12_22Nov/build/plot_day3.py','__name__':'c'};exec(src,g);N=g['NX']
R=[('plain fine','meshes/plain_fine'),('plain coarse','d4quick_pCSu/plain_coarse'),('plain medium','meshes/plain_medium'),('SW test 2','d4test2_fn11/case'),('SW B05','zonefix/B05')]
R+=[('SW '+c,'zonefix/'+c) for c in sys.argv[1:]]
R=[r for r in R if os.path.exists(T0+r[1]+'/23')]
per=5;a=23-per*g['T'];K=('surge','heave','pitch','tension')
Mo={n:g['numbers'](g['cfd'](T0+p),a,per) for n,p in R}
print('MOTION amplitude, last %d periods (%.2f-23 s): error %% vs plain fine | vs experiment'%(per,a))
print('%-13s'%'run'+''.join('%22s'%k for k in K))
for n,_ in R:
    print('%-13s'%n+''.join('%9.1f%% |%7.1f%%   '%(100*(Mo[n][k]['amp']/Mo['plain fine'][k]['amp']-1),100*(Mo[n][k]['amp']/N[k]['amp']-1)) for k in K))
X=(-1.,-.25,.25,1.);tl=[round(a+.1*i,1) for i in range(0,int(per*g['T']*10)+2)];tl=[t for t in tl if a<=t<=23];L={}
for n,p in R:
    P=T0+p;C=cells(P);P0=pts(P+'/constant/polyMesh/points');xc=P0[C].mean(1)[:,0];M={}
    col={x:abs(xc-xc[np.argmin(abs(xc-x))])<1e-5 for x in X}
    for t in tl:
        d='%g'%t
        if not os.path.exists(f'{P}/{d}/alpha.water'):continue
        fp=f'{P}/{d}/polyMesh/points';z=(pts(fp) if os.path.exists(fp) else P0)[C][:,:,2];dz=z.max(1)-z.min(1)
        q=field(f'{P}/{d}/alpha.water');M[t]=[z[col[x]].min()+(q[col[x]]*dz[col[x]]).sum() for x in X]
    L[n]=M
print('WAVE total eta vs plain fine, rms mm, x =',X)
for n,_ in R:
    if n=='plain fine':continue
    q=[t for t in L[n] if t in L['plain fine']];d=np.array([np.subtract(L[n][t],L['plain fine'][t]) for t in q])*1e3
    r=np.sqrt((d*d).mean(0));print('%-13s'%n,' '.join('%6.1f'%v for v in r),'| x+-1 mean %.1f, box +-0.25 mean %.1f  (n=%d)'%(r[[0,3]].mean(),r[[1,2]].mean(),len(q)))