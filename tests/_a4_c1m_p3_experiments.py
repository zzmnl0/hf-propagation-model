"""P3 immutable production experiments, using the existing C1M MATLAB exporter."""
from __future__ import annotations
import os
for _k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_k]='1'
import sys
sys.dont_write_bytecode=True
import argparse, hashlib, json, shutil, subprocess, time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from scipy.io import savemat
import _run_a4_c1m_attribution as common

ROOT=Path(__file__).resolve().parents[1]
PARENT=ROOT/'output2/a4_pharlap/c1/c1m/a2_20260928_v1/a2'
COEFFICIENT=80.61638494312905
EXPORTER=ROOT/'tests/matlab/a4_c1m_uniform_export.m'
MEX=Path('D:/pharlap/pharlap_4.7.4/mex/raytrace_3d_sp.mexw64')

def read(p): return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def write(p,d):
    Path(p).parent.mkdir(parents=True,exist_ok=True)
    common.write_json(Path(p),d)
def digest(d): return hashlib.sha256(json.dumps(d,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()

def prepare(run):
    directory=run/'experiments'
    if directory.exists(): raise FileExistsError(directory)
    directory.mkdir()
    source=directory/'source';source.mkdir()
    for p in (Path(__file__),EXPORTER,ROOT/'tests/_run_a4_c1m_attribution.py',ROOT/'tests/_a4_c1m_p2_analysis_core_v2.py'):
        shutil.copy2(p,source/p.name)
    exposures=[];checked=0
    def scan(d,p):
        if isinstance(d,dict):
            if d.get('frequency_MHz',d.get('freq_MHz'))==12 and d.get('elevation_deg',d.get('elev_deg')) in (25,35): exposures.append(str(p))
            for v in d.values():scan(v,p)
        elif isinstance(d,list):
            for v in d:scan(v,p)
    for p in (ROOT/'output2/a4_pharlap').rglob('*.json'):
        if run in p.parents:continue
        try: d=read(p)
        except (ValueError,UnicodeError):continue
        checked+=1;scan(d,p)
    if exposures: raise ValueError('Requested holdout previously exposed: register a new holdout before proceeding')
    write(directory/'holdout_exposure_audit.json',{'files_checked':checked,'requested_holdout':[{'frequency_MHz':12,'elevation_deg':e} for e in (25,35)],'prior_matches':exposures,'unexposed':True})
    base=read(PARENT/'input_manifest.json')['jobs']
    batches=[('R1',read(PARENT/'identified/input_manifest.json')['jobs']),
             ('R2',[j for j in base if j['stage']=='smooth']),
             ('R3',read(PARENT/'identified_holdout/input_manifest.json')['jobs']),
             ('R5_grid',[j for j in base if j['stage']=='affine' and j['frequency_MHz']==10])]
    hold=[]
    template=next(j for j in base if j['stage']=='smooth')
    for e in (25.,35.):
        for level in (1,2,3):
            j=dict(template,frequency_MHz=12.,elevation_deg=e,level=level,tol=common.TOL[level-1].tolist(),job_id=f'holdout_f12_e{e:g}_t{level}')
            hold.append(j)
    batches.insert(3,('R4',hold))
    jobs=[];grids=directory/'grids';grids.mkdir()
    for role,batch in batches:
        for old in batch:
            job=dict(old);original=Path(job['grid_mat']);target=grids/original.name
            if common.sha(original)!=job['grid_sha256']:raise ValueError('Historical grid changed')
            if not target.exists():shutil.copy2(original,target)
            job.update(grid_mat=str(target.resolve()),launch_k2=COEFFICIENT,role=role,
                       scope='natural_full_path' if role in ('R2','R4') else 'first_100_km',
                       spatial_fd_base_km=.125,compatibility_mode='pharlap_radial_b0',
                       exporter_sha256=common.sha(EXPORTER),p1_case='P3',parent_job_sha256=old['job_sha256'])
            job.pop('job_sha256',None);job['job_sha256']=digest(job);jobs.append(job)
    assert len(jobs)==168 and sum(j['role']!='R5_grid' for j in jobs)==72
    fd=[]
    probes=[j for j in jobs if j['role']=='R1']+[j for j in jobs if j['role']=='R2' and j['level']==3 and j['elevation_deg'] in (21.94290326504532,39.58701138438118)]
    for j in probes:
        for h in (.5,.25):
            if j['role']=='R1' and j['level']!=3 and h==.25:continue
            q=dict(j,job_id=j['job_id']+f'_fd{h:g}'.replace('.','p'),role='R5_fd',ph_job_id=j['job_id'],spatial_fd_base_km=h)
            q.pop('job_sha256');q['job_sha256']=digest(q);fd.append(q)
    alljobs=jobs+fd
    payload=np.empty((1,len(jobs)),dtype=object)
    for i,j in enumerate(jobs):payload[0,i]=j
    savemat(directory/'input.mat',{'a2_jobs':payload},long_field_names=True)
    for name,subset in [('training',[j for j in jobs if j['role']!='R4']),('holdout',[j for j in jobs if j['role']=='R4'])]:
        part=np.empty((1,len(subset)),dtype=object)
        for i,j in enumerate(subset):part[0,i]=j
        savemat(directory/(name+'.mat'),{'a2_jobs':part},long_field_names=True)
    write(directory/'manifest.json',{'schema':'c1m-p3-production-v1','coefficient':COEFFICIENT,'jobs':alljobs,'ph_jobs':len(jobs),'ah_jobs':len(alljobs),'core_ph_jobs':72,'ah_profiles':common.AH_PROFILES,'source_sha256':{str(p):common.sha(p) for p in source.iterdir()},'input_sha256':common.sha(directory/'input.mat'),'formal_method':'physical_derivative_hermite','reconstruction_scope':'B0 collisionless only','ph_reused':False})
    print(json.dumps({'prepared':len(alljobs),'ph':len(jobs)}),flush=True)

def run_ph(run,holdout=False):
    d=run/'experiments';m=read(d/'manifest.json')
    if common.sha(d/'input.mat')!=m['input_sha256']:raise ValueError('Input changed')
    exporter=d/'source/a4_c1m_uniform_export.m'
    if common.sha(exporter)!=m['jobs'][0]['exporter_sha256']:raise ValueError('Exporter changed')
    q=lambda p:str(p.resolve()).replace("'","''").replace('\\','/')
    batch='holdout' if holdout else 'training'
    if holdout and not (d/'production_parameter_freeze.json').exists():raise ValueError('Freeze production source and parameters before holdout')
    command=[str(common.MATLAB_EXE),'-batch',f"addpath('{q(exporter.parent)}');a4_c1m_uniform_export('{q(d/(batch+'.mat'))}','{q(d/'pharlap_raw')}','{q(d/(batch+'_pharlap_summary.json'))}')"]
    with (d/'matlab.log').open('a',encoding='utf-8') as log:
        result=subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
    if result.returncode:raise RuntimeError(f'MATLAB exit {result.returncode}; see {d}/matlab.log')

def ah_one(job,directory,source_hashes):
    from models.ne_provider_geocentric import GeocentricGridNeProvider3D
    from models.ray_tracer_ah_3d_ecef import RefractiveIndexAH3DECEF,integrate_fixed_path_ah3d_ecef,trace_ray_ah3d_ecef
    from config import UPDATE2
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    dest=directory/(job['job_id']+'.json');raw=dest.with_suffix('.npz')
    identity={'job_sha256':job['job_sha256'],'source_sha256':source_hashes}
    if dest.exists():
        saved=read(dest)
        if any(saved.get(k)!=v for k,v in identity.items()):raise ValueError('AH resume identity changed; use new run')
        if saved['status']=='completed' and common.sha(raw)==saved['raw_sha256']:return saved
        raise ValueError('Existing failed raw result immutable; use new attempt directory')
    start=time.perf_counter();result=dict(identity,job_id=job['job_id'])
    try:
        if common.sha(Path(job['grid_mat']))!=job['grid_sha256']:raise ValueError('Grid changed')
        g=common.mat(Path(job['grid_mat']))
        provider=GeocentricGridNeProvider3D(g['Ne_m3'],g['latitude_deg'],g['longitude_deg'],g['height_km'],earth_radius_km=6371.,interpolation='radial_hermite')
        name='p3_'+str(job['level']);UPDATE2['AH3D_INTEGRATION'][name]=dict(common.AH_PROFILES[job['level']-1])
        model=RefractiveIndexAH3DECEF(provider,earth_radius_km=6371.,freq_MHz=job['frequency_MHz'],wave_mode='O',gyro_MHz=0.,collision_s=0.,integration_profile=name,spatial_fd_base_km=job['spatial_fd_base_km'],compatibility_mode='pharlap_radial_b0',plasma_coefficient=COEFFICIENT)
        args=dict(tx_lat_deg=30.,tx_lon_deg=120.,tx_height_km=job['origin_height_km'],elevation_deg=job['elevation_deg'],bearing_deg=0.,profile=name)
        if job['scope']=='first_100_km':trace=integrate_fixed_path_ah3d_ecef(model,**args,path_length_km=100.,stop_on_physical_events=False,apply_dispersion_projection=True)
        else:trace=trace_ray_ah3d_ecef(model,**args,h_stop_km=460.,max_path_km=5000.)
        xyz=np.asarray(trace['trajectory_ecef_km']);p=np.asarray(trace['wavevector_ecef_km^-1'])/(2*np.pi*job['frequency_MHz']*1e6/299792.458)
        s=np.asarray(trace['trajectory_geom_path_km']);group=np.asarray(trace['trajectory_group_path_km'])
        llh=np.asarray(trace['trajectory_lat_lon_height']);ne=np.array([provider.get_Ne_geocentric(*row) for row in llh])
        np.savez_compressed(raw,s=s,xyz=xyz,p=p,group=group,tangent=p/np.linalg.norm(p,axis=1)[:,None],llh=llh,ne=ne,chord_tangent=np.asarray(trace['trajectory_group_tangent_ecef']))
        result.update(status='completed',raw_sha256=common.sha(raw),termination=trace.get('termination_status'),failure_reason=trace.get('failure_reason'),group_path_km=trace.get('group_path_km'),geom_path_km=float(s[-1]),point_count=len(s),initial_residual=float(p[0]@p[0]-(1-COEFFICIENT*ne[0]/(job['frequency_MHz']*1e6)**2)),metadata=model.metadata(),group_path_semantics=trace.get('trajectory_group_path_semantics'))
    except Exception as e:result.update(status='failed',error=repr(e))
    result['elapsed_s']=time.perf_counter()-start;write(dest,result);return result

def run_ah(run,workers,roles):
    d=run/'experiments';m=read(d/'manifest.json');files=list((ROOT/'models').glob('*.py'))+[ROOT/'config.py',ROOT/'utils.py',Path(__file__),ROOT/'tests/_run_a4_c1m_attribution.py']
    hashes={str(p):common.sha(p) for p in files};freeze=d/'ah_source'
    if not freeze.exists():
        freeze.mkdir()
        for p in files:shutil.copy2(p,freeze/p.name)
        write(d/'ah_source_sha256.json',hashes)
    elif read(d/'ah_source_sha256.json')!=hashes:raise ValueError('AH source changed after freeze')
    jobs=[j for j in m['jobs'] if not roles or j['role'] in roles]
    if any(j['role']=='R4' for j in jobs) and not (d/'production_parameter_freeze.json').exists():raise ValueError('Holdout requires production parameter freeze')
    with ProcessPoolExecutor(max_workers=min(2,workers)) as pool:
        futures={pool.submit(ah_one,j,str(d/'ah'),hashes):j['job_id'] for j in jobs}
        for f in as_completed(futures):
            r=f.result();print(json.dumps({'job':r['job_id'],'status':r['status'],'elapsed_s':r['elapsed_s'],'error':r.get('error')}),flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('stage',choices=('prepare','ph','ah'));p.add_argument('--run',type=Path,required=True);p.add_argument('--workers',type=int,default=2);p.add_argument('--roles',nargs='*');p.add_argument('--holdout',action='store_true');a=p.parse_args()
    if a.stage=='prepare':prepare(a.run.resolve())
    elif a.stage=='ph':run_ph(a.run.resolve(),a.holdout)
    else:run_ah(a.run.resolve(),a.workers,a.roles)

if __name__=='__main__':main()
