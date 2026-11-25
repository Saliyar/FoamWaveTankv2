"""Builds the theory-matrix cases zonefix/th/<case> from hr/C_sf and th/gen/constant/waveProperties (setWaveParameters
output). After generation the streamFunction cases need Niter 40, direction (1 0 0) and period 1.993 added back (the
setWaveProperties output drops them), the irregular cases the farfieldCoeffs Air2Air block removed."""
import re,sys,os,shutil,subprocess
Z='/gpfs/hpc-home/projects/dtu_00075/people/sitali/FoamWaveTank_Training_test/zonefix'
src=Z+'/hr/C_sf'
gen=open(Z+'/th/gen/constant/waveProperties').read()
def block(s,name):
    m=re.search(r'\n'+name+r'\s*\n\{',s);i=m.end();d=1;j=i
    while d: d+=(s[j]=='{')-(s[j]=='}');j+=1
    return m.start()+1,j,s[i:j-1]
def make(case,coeff=None,ext=None,set_={},end=15,dt=None,mesh=None):
    dst=Z+'/th/'+case
    if os.path.exists(dst): shutil.rmtree(dst)
    os.makedirs(dst)
    for f in ('0.orig','constant','system','Allrun'):
        (shutil.copytree if os.path.isdir(src+'/'+f) else shutil.copy)(src+'/'+f,dst+'/'+f)
    wp=open(dst+'/constant/waveProperties').read()
    a,b,_=block(wp,'waveTheoryForcingCoeffs')
    if coeff:
        _,_,body=block(gen,coeff+'Coeffs')
        if 'phi ' not in body and 'cn_' not in coeff and 'ir_' not in coeff: body+='    phi 4.71238898038469;\n'
        wp=wp[:a]+'waveTheoryForcingCoeffs\n{'+body+'}'+wp[b:]
    if ext:   # irregularFast
        _,_,body=block(gen,ext+'Coeffs')
        wp=wp[:a]+'externalForcingCoeffs\n{'+body+'}'+wp[b:]
        set_=dict(set_,externalForcing='irregularFast',ignoreMeshMotion='no',searchTolerance='1e-6')
    for k,v in set_.items():
        if re.search(r'\n'+k+r'\s',wp): wp=re.sub(r'\n'+k+r'\s[^;]*;','\n%s %s;'%(k,v),wp,count=1)
        else: wp=wp.replace('\nseaLevel','\n%s %s;\nseaLevel'%(k,v),1)
    open(dst+'/constant/waveProperties','w').write(wp)
    cd=open(dst+'/system/controlDict').read()
    cd=re.sub(r'\nendTime\s[^;]*;','\nendTime         %g;'%end,cd)
    if dt: cd=re.sub(r'\ndeltaT\s[^;]*;','\ndeltaT          %g;'%dt,cd)
    open(dst+'/system/controlDict','w').write(cd)
    al=open(dst+'/Allrun').read().replace('fwt_zone_fix/private_hr','fwt_zone_fix/private_th')
    if mesh: al=al.replace('MESH=${MESH:-coarse}','MESH=%s'%mesh)
    open(dst+'/Allrun','w').write(al)
    print(case)
common={'incidentResidualMax':'1','incidentEtaProbes':'((-1 0 0) (1 0 0))','rampTsoft':'6'}
for w in ('mild','ours'):
    for t in ('s1','s2','s5','sf'): make(f'{t}_{w}',coeff=f'{t}_{w}',set_=common)
    make(f'cn_{w}',coeff=f'cn_{w}',set_=common,end=0.01)
    make(f'ir_{w}',coeff=f'ir_{w}',set_=dict(common,rampMode='factor',airFill='water'))
make('irF_mild',ext='ir_mild',set_=dict(common,rampMode='factor',airFill='water'),end=0.01)
g={'incidentResidualMax':'1','incidentEtaProbes':'((-1 0 0) (1 0 0))','rampTsoft':'0'}
make('G1_coarse_noramp',coeff='sf_ours',set_=g,end=10)
make('G2_medium_noramp',coeff='sf_ours',set_=g,end=10,mesh='medium')
make('G3_coarse_noramp_dt',coeff='sf_ours',set_=g,end=10,dt=0.0005)
