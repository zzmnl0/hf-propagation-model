"""P3.1 immutable entrance controls; no production compatibility changes."""
from __future__ import annotations
import argparse, csv, hashlib, json, os, shutil, subprocess
from pathlib import Path
import numpy as np
from scipy.io import loadmat, savemat
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator
from scipy.integrate import solve_ivp
import _run_a4_c1m_attribution as common

ROOT=Path(__file__).resolve().parents[1]
C=80.61638494312905
EXPORT=ROOT/'tests/matlab/a4_c1m_uniform_export.m'
def read(p): return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def write(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2,allow_nan=False),encoding='utf-8')
def digest(x): return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def medium_spline(m):
    h=m['height_km'];v=m['ne_profile_m3'];return CubicHermiteSpline(h,v,np.gradient(v,h,edge_order=2))
def prepare(run):
    d=run/'boundary'
    if d.exists(): raise FileExistsError(d)
    (d/'source').mkdir(parents=True);(d/'grids').mkdir()
    for p in (Path(__file__),EXPORT,ROOT/'tests/_run_a4_c1m_attribution.py'):shutil.copy2(p,d/'source'/p.name)
    original=ROOT/'output2/a4_pharlap/c1/c1m/p3_20261009_v1/experiments/grids/smooth_s5.mat'
    m={k:v for k,v in loadmat(original,simplify_cells=True).items() if not k.startswith('__')}
    jobs=[];up=np.array([-.4330127018922193,.75,.5]);north=np.array([.25,-.4330127018922193,.8660254037844386])
    for amplitude in (0,1,10,100):
        mm=dict(m);mm['Ne_m3']=m['Ne_m3']*amplitude;mm['ne_profile_m3']=m['ne_profile_m3']*amplitude
        gp=d/'grids'/f'smooth_s5_x{amplitude}.mat';savemat(gp,mm);sp=medium_spline(mm)
        for f in ((10,8,15) if amplitude==1 else (10,)):
            for h,e in ((16,22),(16,39.5),(24,-22),(24,-39.5)):
                n2=1-C*float(sp(h))/(f*1e6)**2;p=np.sqrt(n2)*(np.cos(np.deg2rad(e))*north+np.sin(np.deg2rad(e))*up)
                for level,tol in enumerate(common.TOL,1):
                    name=f'boundary_x{amplitude}_f{f}_h{h}_e{e:g}_t{level}'.replace('.','p').replace('-','m')
                    j=dict(job_id=name,case_id=name.rsplit('_t',1)[0],stage='p31_boundary',p1_case='P31',grid_mat=str(gp.resolve()),grid_sha256=common.sha(gp),frequency_MHz=float(f),elevation_deg=float(e),origin_height_km=float(h),origin_lat_deg=30.,origin_lon_deg=120.,initial_xyz_km=((6371+h)*up).tolist(),initial_p=p.tolist(),explicit_state=True,launch_k2=C,level=level,tol=tol.tolist(),amplitude=amplitude,scope_km=20.,role='training' if f==10 else 'holdout',exporter_sha256=common.sha(EXPORT))
                    j['job_sha256']=digest(j);jobs.append(j)
    for role in ('training','holdout'):
        part=[j for j in jobs if j['role']==role];a=np.empty((1,len(part)),object)
        for i,j in enumerate(part):a[0,i]=j
        savemat(d/f'{role}.mat',{'a2_jobs':a},long_field_names=True)
    write(d/'manifest.json',dict(schema='p31-boundary-v1',jobs=jobs,parent_grid_sha256=common.sha(original),source_sha256={p.name:common.sha(p) for p in (d/'source').iterdir()},input_sha256={r:common.sha(d/f'{r}.mat') for r in ('training','holdout')},coefficient=C))
def ph(run,role):
    d=run/'boundary';m=read(d/'manifest.json')
    for name,hash_ in m['source_sha256'].items():
        if common.sha(d/'source'/name)!=hash_:raise ValueError('Frozen source mismatch')
    if role=='holdout' and not (d/'training_analysis/candidate_freeze.json').exists():raise ValueError('Freeze training candidate status first')
    if common.sha(d/f'{role}.mat')!=m['input_sha256'][role]:raise ValueError('Input mismatch')
    q=lambda p:str(p.resolve()).replace('\\','/').replace("'","''")
    expr=f"addpath('{q(d/'source')}');a4_c1m_uniform_export('{q(d/f'{role}.mat')}','{q(d/'pharlap_raw')}','{q(d/f'{role}_summary.json')}')"
    with (d/f'{role}_matlab.log').open('a',encoding='utf-8') as log:
        x=subprocess.run([str(common.MATLAB_EXE),'-batch',expr],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
    if x.returncode:raise RuntimeError(f'MATLAB returned {x.returncode}')
def candidate(sp,job,kind,s):
    a=C/(job['frequency_MHz']*1e6)**2;der=sp.derivative()
    def rhs(_,y):
        r=y[:3];p=y[3:6];radius=np.linalg.norm(r);h=radius-6371.;n=np.linalg.norm(p);slope=float(der(h))
        if kind=='clip' and float(sp(h))<=0:slope=0.
        return np.r_[p/n,-a*slope*r/radius/(2*n),1/n]
    y=np.r_[job['initial_xyz_km'],job['initial_p'],0.]
    sol=solve_ivp(rhs,(0,20),y,method='DOP853',rtol=2.3e-14,atol=1e-13,max_step=.02,dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    return sol.sol(s).T
def analyze(run,role):
    d=run/'boundary';out=d/f'{role}_analysis'
    if out.exists():raise FileExistsError(out)
    out.mkdir();records=[];rows=[]
    for j in read(d/'manifest.json')['jobs']:
        if j['role']!=role:continue
        path=d/'pharlap_raw'/f"{j['job_id']}.mat";raw=loadmat(path,simplify_cells=True);st=raw['ray_state_vec'];rp=raw['ray_path_data']
        v=lambda k:common.vector(st,k)
        s=v('geometrical_path')/1000.;mask=(s>=0)&(s<=20+1e-9);s=s[mask]
        xyz=np.column_stack([v(k) for k in ('pos_x','pos_y','pos_z')])[mask]/1000.;p=np.column_stack([v(k) for k in ('dir_x','dir_y','dir_z')])[mask];h=np.linalg.norm(xyz,axis=1)-6371.;p2=np.sum(p*p,axis=1);g=v('group_path')[mask]/1000.
        m=loadmat(j['grid_mat'],simplify_cells=True);sp=medium_spline(m);ne=sp(h);neout=common.vector(rp,'electron_density')[mask]*1e6;nout=common.vector(rp,'refractive_index')[mask];a=C/(j['frequency_MHz']*1e6)**2
        pred={k:candidate(sp,j,k,s) for k in ('raw','clip')};scores={}
        for k,y in pred.items():scores[k]=dict(max_position_m=float(np.max(np.linalg.norm(y[:,:3]-xyz,axis=1))*1000),max_p_error=float(np.max(np.linalg.norm(y[:,3:6]-p,axis=1))),max_group_m=float(np.max(abs(y[:,6]-g))*1000))
        neg=ne<0
        rec=dict(job_id=j['job_id'],case_id=j['case_id'],level=j['level'],amplitude=j['amplitude'],frequency_MHz=j['frequency_MHz'],elevation_deg=j['elevation_deg'],native_points=len(s),native_negative_points=int(sum(neg)),max_output_clip_ne_error=float(np.max(abs(neout-np.maximum(ne,0)))),max_raw_residual=float(np.max(abs(p2-(1-a*ne)))),max_clip_residual=float(np.max(abs(p2-(1-a*np.maximum(ne,0))))),initial_p_error=float(np.linalg.norm(p[0]-j['initial_p'])),raw_sha256=common.sha(path),candidate_scores=scores)
        records.append(rec)
        for i in range(len(s)):rows.append(dict(job_id=j['job_id'],s_km=s[i],h_km=h[i],raw_ne_m3=ne[i],output_ne_m3=neout[i],p2=p2[i],output_n2=nout[i]**2,raw_residual=p2[i]-1+a*ne[i],clip_residual=p2[i]-1+a*max(ne[i],0),group_km=g[i]))
    with (out/'pointwise.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    write(out/'metrics.json',records)
    write(out/'candidate_freeze.json',dict(output_rule='compare_centered_Hermite_clipped_zero',propagation_rule='unresolved_pending_state_and_convergence_review',boundary_convention_identified=False,production_patch_authorized_by_evidence=False,candidates=['clipped_value_and_gradient','raw_propagation_clipped_output','mixed_or_state_reset'],note='Raw and mixed-gradient laws can coincide in p dynamics; density output alone cannot distinguish them.'))
    print(json.dumps({'role':role,'records':len(records),'analysis':str(out)}))
def diagnose(run):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    d=run/'boundary';out=d/'diagnostics'
    if out.exists():raise FileExistsError(out)
    out.mkdir();shutil.copy2(__file__,out/Path(__file__).name)
    data=[];fd=[];resets=[]
    for j in read(d/'manifest.json')['jobs']:
        if j['level']!=3:continue
        raw=loadmat(d/'pharlap_raw'/f"{j['job_id']}.mat",simplify_cells=True);st=raw['ray_state_vec'];rp=raw['ray_path_data'];s=np.asarray(st['geometrical_path'])/1000
        xyz=np.column_stack([st[k] for k in ('pos_x','pos_y','pos_z')])/1000;p=np.column_stack([st[k] for k in ('dir_x','dir_y','dir_z')]);h=np.linalg.norm(xyz,axis=1)-6371;p2=np.sum(p*p,axis=1);m=loadmat(j['grid_mat'],simplify_cells=True);sp=medium_spline(m);ne=sp(h);a=C/(j['frequency_MHz']*1e6)**2;phase=np.asarray(rp['refractive_index'])**2
        mask=(s<=20);negative=mask&(ne<0);nodes=np.flatnonzero(negative)
        mismatches=[i for i in nodes if abs(p2[i]-phase[i])>1e-13]
        for i in mismatches:resets.append(dict(job_id=j['job_id'],index=int(i),h_km=float(h[i]),s_km=float(s[i]),raw_ne=float(ne[i]),p2_minus1=float(p2[i]-1),output_n2_minus1=float(phase[i]-1),next_p2_minus1=float(p2[i+1]-1),next_h_km=float(h[i+1])))
        # This tests a directional derivative of reconstructed states, not a full gradient.
        ip=PchipInterpolator(s,p,axis=0);ix=PchipInterpolator(s,xyz,axis=0)
        q=np.arange(1.,19.01,.25);rr=ix(q);hh=np.linalg.norm(rr,axis=1)-6371;pp=ip(q);norm=np.linalg.norm(pp,axis=1);grad=sp.derivative()(hh)[:,None]*rr/np.linalg.norm(rr,axis=1)[:,None];expected=-a*grad/(2*norm[:,None]);expected_clip=expected.copy();expected_clip[sp(hh)<=0]=0;previous=None
        for step in (.5,.25,.125):
            observed=(ip(q+step)-ip(q-step))/(2*step)
            fd.append(dict(job_id=j['job_id'],difference_halfwidth_km=step,raw_directional_force_max_error=float(np.max(np.linalg.norm(observed-expected,axis=1))),clip_directional_force_max_error=float(np.max(np.linalg.norm(observed-expected_clip,axis=1))),successive_derivative_change=None if previous is None else float(np.max(np.linalg.norm(observed-previous,axis=1))),full_gradient_claim=False,independently_validated_reconstruction=False));previous=observed
        data.append(dict(job_id=j['job_id'],negative_native_count=len(nodes),mismatched_state_output_count=len(mismatches),max_negative_p2_minus1=float(max(abs(p2[negative]-1),default=0)),max_negative_raw_prediction=float(max(abs(a*ne[negative]),default=0))))
        if j['amplitude']==100 or (j['amplitude']==1 and j['frequency_MHz']!=10):
            fig,axs=plt.subplots(2,1,figsize=(8,7),sharex=True)
            axs[0].plot(s[mask],p2[mask]-1,'o-',label='PH p2 - 1');axs[0].plot(s[mask],phase[mask]-1,'.-',label='PH output n2 - 1');axs[0].plot(s[mask],-a*ne[mask],label='signed Hermite prediction');axs[0].plot(s[mask],-a*np.maximum(ne[mask],0),label='clipped prediction');axs[0].legend(fontsize=8);axs[0].set_ylabel('Index squared residual')
            axs[1].plot(s[mask],ne[mask],label='raw Hermite Ne');axs[1].plot(s[mask],np.asarray(rp['electron_density'])[mask]*1e6,'.',label='PH Ne');axs[1].legend();axs[1].set_xlabel('Native geometric path (km)');axs[1].set_ylabel('Ne (m^-3)');fig.suptitle(j['job_id']);fig.tight_layout();fig.savefig(out/(j['job_id']+'.png'),dpi=140);plt.close(fig)
    write(out/'state_reset_evidence.json',resets);write(out/'directional_difference_diagnostics.json',fd);write(out/'negative_support_summary.json',data)
    metrics=read(d/'training_analysis/metrics.json')+read(d/'holdout_analysis/metrics.json')
    report=dict(schema='p31-boundary-conclusion-v1',jobs=len(metrics),all_custom_initial_states_preserved=bool(all(x['initial_p_error']==0 for x in metrics)),boundary_convention_identified=False,output_rule_identified=True,propagation_rule_identified=False,production_repair_permitted=False,causes=[dict(candidate='raw propagation / raw-gradient mixed observable class',status='excluded_as_complete_mechanism',evidence='In negative Hermite support most PH states have p2=1, contradicting signed-density invariant.'),dict(candidate='pure clipped Hamiltonian with no output/state reset',status='unresolved_not_sufficient',evidence='PH output n2=1 can precede state p2 reset by one native sample; finite-resolution crossing errors are not isolated.'),dict(candidate='vacuum classification and state normalization / delayed output semantics',status='confirmed_observation_unresolved_internal_mechanism',evidence='Native state-output mismatches recur in amplified training and held-out data. Normalization location in the closed MEX remains unknown.')],state_output_mismatch_count=len(resets),gradient_diagnostics_only=True,scope_km=20.,recommend_p4=False)
    write(out/'conclusion.json',report)
    (out/'boundary_report.txt').write_text('P3.1-A入口归因\n72次PH调用，48训练之后冻结未解决状态，再运行24留出。\n全部共同初态逐值保持：'+str(report['all_custom_initial_states_preserved'])+'\nNe输出符合截零中心斜率Hermite，但这不足以确定内部梯度。\n负支撑区大量状态p²=1，排除未重置的纯有符号传播作为完整机制。\n负区首节点存在输出n²=1而p²非1、后续节点恢复1的现象；状态投影位置与输出语义尚未唯一归因。\n已执行三个方向差分宽度；基于PCHIP状态重建，仅作为诊断，不宣称完整梯度或独立梯度收敛。\n不实施生产裁剪或经验重置；P4准入保持false。\n',encoding='utf-8')
    write(out/'sha256.json',{str(p.relative_to(d)):common.sha(p) for p in d.rglob('*') if p.is_file() and p.name!='sha256.json'})
def verify(run):
    d=run/'boundary';version=1
    while (d/f'verification_v{version}').exists():version+=1
    out=d/f'verification_v{version}'
    out.mkdir();shutil.copy2(__file__,out/Path(__file__).name);m=read(d/'manifest.json');checks=[]
    def check(name,passed):checks.append(dict(name=name,passed=bool(passed)))
    check('matrix_72',len(m['jobs'])==72)
    check('training_48',sum(j['role']=='training' for j in m['jobs'])==48)
    for name,h in m['source_sha256'].items():check('frozen_source:'+name,common.sha(d/'source'/name)==h)
    for role,h in m['input_sha256'].items():check('input:'+role,common.sha(d/f'{role}.mat')==h)
    for j in m['jobs']:
        q=dict(j);ident=q.pop('job_sha256');check(j['job_id']+':identity',digest(q)==ident)
        check(j['job_id']+':grid_hash',common.sha(Path(j['grid_mat']))==j['grid_sha256'])
        raw=loadmat(d/'pharlap_raw'/f"{j['job_id']}.mat",simplify_cells=True);st=raw['ray_state_vec'];check(j['job_id']+':raw_identity',raw['job_sha256']==ident and raw['grid_sha256']==j['grid_sha256'] and raw['exporter_sha256']==j['exporter_sha256'])
        s=np.asarray(st['geometrical_path'])/1000;check(j['job_id']+':full20_covered',s[0]==0 and s[-1]>=20 and np.all(np.diff(s)>0))
        p=np.array([np.asarray(st[k])[0] for k in ('dir_x','dir_y','dir_z')]);xyz=np.array([np.asarray(st[k])[0] for k in ('pos_x','pos_y','pos_z')])/1000;check(j['job_id']+':initial_p_exact',np.array_equal(p,j['initial_p']));check(j['job_id']+':initial_position',np.max(abs(xyz-j['initial_xyz_km']))<1e-11)
    # The candidate implementation must reproduce straight vacuum before it is used diagnostically.
    j=next(j for j in m['jobs'] if j['amplitude']==0);sp=medium_spline(loadmat(j['grid_mat'],simplify_cells=True));s=np.array([0.,5.,20.]);y=candidate(sp,j,'clip',s);direction=np.array(j['initial_p'])/np.linalg.norm(j['initial_p']);expected=np.array(j['initial_xyz_km'])+s[:,None]*direction
    check('candidate_vacuum_position',np.max(abs(y[:,:3]-expected))<1e-8);check('candidate_vacuum_group',np.max(abs(y[:,6]-s))<1e-10)
    check('training_before_holdout', (d/'training_analysis/candidate_freeze.json').stat().st_mtime_ns<(d/'holdout_summary.json').stat().st_mtime_ns)
    summary=dict(checks=checks,evidence_checks_passed=all(c['passed'] for c in checks),boundary_convention_identified=False,formal_20km_endpoint_acceptance=False,native_node_diagnostics_cover_requested_interval=True,scope_note='All native outputs bracket 20 km. Candidate scores are diagnostic at native nodes; no strict trajectory acceptance is claimed.',source_sha256=common.sha(Path(__file__)))
    write(out/'verification.json',summary);write(out/'sha256.json',{str(p.relative_to(d)):common.sha(p) for p in d.rglob('*') if p.is_file() and p.name!='sha256.json'})
    print(json.dumps({'passed':summary['evidence_checks_passed'],'checks':len(checks),'failures':[x for x in checks if not x['passed']]}))
def main():
    p=argparse.ArgumentParser();p.add_argument('stage',choices=('prepare','ph','analyze','diagnose','verify'));p.add_argument('--run',type=Path,required=True);p.add_argument('--role',choices=('training','holdout'),default='training');a=p.parse_args()
    if a.stage=='prepare':prepare(a.run)
    elif a.stage=='ph':ph(a.run,a.role)
    elif a.stage=='analyze':analyze(a.run,a.role)
    elif a.stage=='diagnose':diagnose(a.run)
    else:verify(a.run)
if __name__=='__main__':main()
