import math, os
F="/gpfs/hpc-home/projects/dtu_00075/apps/FoamWaveTank/FoamWaveTankv2/tutorials/SWENSE_bc_template/constant/modes_HOS_SWENSE_clean.dat"
n1=512; reclen=18*n1
nt=os.path.getsize(F)//(reclen*8)-1
h=5.0; dtd=0.05
eta=[]
with open(F,'rb') as f:
    for it in range(nt):
        f.seek(((it+1)*8+4)*reclen); b=f.read(18*n1)
        s=0.0
        for i in range(n1):
            s+=float(b[18*i:18*i+18].replace(b'D',b'E').replace(b'd',b'E'))
        eta.append(s*h)
t=[i*dtd for i in range(nt)]
m=sum(eta)/nt
x=[e-m for e in eta]
# period from the first upward zero crossings of the mean-removed signal
zc=[t[i]+(t[i+1]-t[i])*(-x[i])/(x[i+1]-x[i]) for i in range(nt-1) if x[i]<0<=x[i+1]]
T=(zc[-1]-zc[0])/(len(zc)-1) if len(zc)>1 else float('nan')
print("nt",nt,"span",t[-1],"eta(0) max,min",max(eta),min(eta))
print("zero-crossings",len(zc),"mean period T (s)",T)
# least-squares first harmonic at period T over the whole record
w=2*math.pi/T
S=[[0.0]*3 for _ in range(3)]; r=[0.0]*3
for i in range(nt):
    b=[math.cos(w*t[i]),math.sin(w*t[i]),1.0]
    for p in range(3):
        r[p]+=b[p]*eta[i]
        for q in range(3): S[p][q]+=b[p]*b[q]
# solve 3x3
import copy
A=copy.deepcopy(S); y=r[:]
for p in range(3):
    piv=A[p][p]
    for q in range(p+1,3):
        fct=A[q][p]/piv
        for k in range(3): A[q][k]-=fct*A[p][k]
        y[q]-=fct*y[p]
c=[0,0,0]
for p in (2,1,0):
    c[p]=(y[p]-sum(A[p][k]*c[k] for k in range(p+1,3)))/A[p][p]
a1=math.hypot(c[0],c[1])
print("first harmonic amplitude a1 (m)",a1,"H = 2 a1 (m)",2*a1,"mean level",c[2])
