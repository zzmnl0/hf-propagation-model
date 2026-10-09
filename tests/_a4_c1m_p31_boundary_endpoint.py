"""Exact 20 km supplemental entrance diagnostics, without state-tangent reconstruction."""
from pathlib import Path
import argparse, csv, json, shutil
import numpy as np
from scipy.io import loadmat
from scipy.interpolate import PchipInterpolator
from scipy.integrate import solve_ivp
import _a4_c1m_p31_boundary as base

def arrays(path):
    raw=loadmat(path,simplify_cells=True);st=raw['ray_state_vec'];s=np.asarray(st['geometrical_path'])/1000
    y=np.column_stack([st[k] for k in ('pos_x','pos_y','pos_z')]+[st['group_path']])/1000
    p=np.column_stack([st[k] for k in ('dir_x','dir_y','dir_z')]);n2=np.asarray(raw['ray_path_data']['refractive_index'])**2
    return s,y,p,n2

def stencil(s):
    right=int(np.searchsorted(s,20,side='left'))
    if right==0 or right+1>=len(s):raise ValueError('Real endpoint bracket and following real neighbor required')
    return slice(max(0,right-2),right+2),right

def reconstruct(s,y):
    sl,right=stencil(s)
    pchip=PchipInterpolator(s[sl],y[sl],axis=0)(20.)
    frac=(20-s[right-1])/(s[right]-s[right-1]);chord=y[right-1]+frac*(y[right]-y[right-1])
    return pchip,chord

def solution(sp,j,end):
    a=base.C/(j['frequency_MHz']*1e6)**2;der=sp.derivative()
    def build(kind):
        def rhs(_,y):
            r=y[:3];p=y[3:6];rad=np.linalg.norm(r);h=rad-6371;n=np.linalg.norm(p);grad=float(der(h))
            if kind=='clip' and float(sp(h))<=0:grad=0.
            return np.r_[p/n,-a*grad*r/(2*n*rad),1/n]
        result=solve_ivp(rhs,(0,end),np.r_[j['initial_xyz_km'],j['initial_p'],0.],method='DOP853',rtol=2.3e-14,atol=1e-13,max_step=.02,dense_output=True)
        if not result.success:raise RuntimeError(result.message)
        return result
    return {kind:build(kind) for kind in ('raw','clip')}

def run(run):
    d=run/'boundary';out=d/'endpoint_supplement_v1'
    if out.exists():raise FileExistsError(out)
    out.mkdir();shutil.copy2(__file__,out/Path(__file__).name);shutil.copy2(base.__file__,out/Path(base.__file__).name)
    jobs=base.read(d/'manifest.json')['jobs'];cache={};records=[];controls=[]
    for j in jobs:
        s,y,p,n2=arrays(d/'pharlap_raw'/f"{j['job_id']}.mat");sl,right=stencil(s);ss=s[sl];ph,chord=reconstruct(s,y)
        if j['case_id'] not in cache:
            siblings=[q for q in jobs if q['case_id']==j['case_id']];end=max(arrays(d/'pharlap_raw'/f"{q['job_id']}.mat")[0][stencil(arrays(d/'pharlap_raw'/f"{q['job_id']}.mat")[0])[1]+1] for q in siblings)
            sp=base.medium_spline(loadmat(j['grid_mat'],simplify_cells=True));cache[j['case_id']]=(sp,solution(sp,j,float(end)))
        sp,sols=cache[j['case_id']];hhh=np.linalg.norm(y[sl,:3],axis=1)-6371
        # Ne roots and state/output mismatch are diagnostic indicators only.
        roots=sp.roots(extrapolate=False);cross=bool(np.any((roots>=min(hhh))&(roots<=max(hhh))))
        mismatch=bool(np.any(abs(np.sum(p[sl]**2,axis=1)-n2[sl])>1e-13))
        item=dict(job_id=j['job_id'],level=j['level'],endpoint_s_km=20.,ph_endpoint_xyz_km=ph[:3].tolist(),ph_endpoint_group_km=float(ph[3]),real_bracket_km=[float(s[right-1]),float(s[right])],real_following_neighbor_km=float(s[right+1]),reconstruction='PCHIP_position_and_native_group_only',state_derivatives_used=False,stencil_contains_density_root=cross,stencil_state_output_mismatch=mismatch,pchip_chord_position_gap_m=float(np.linalg.norm(ph[:3]-chord[:3])*1000),pchip_chord_group_gap_m=float(abs(ph[3]-chord[3])*1000),method_gap_is_uncertainty_bound=False,reconstruction_independently_certified=False,formal_endpoint_passed=False,candidates={})
        for kind,sol in sols.items():
            exact=sol.sol(20)[[0,1,2,6]];nodes=sol.sol(ss)[[0,1,2,6],:].T;recon=PchipInterpolator(ss,nodes,axis=0)(20)
            item['candidates'][kind]=dict(endpoint_xyz_km=exact[:3].tolist(),endpoint_group_km=float(exact[3]),ph_endpoint_difference_m=float(np.linalg.norm(ph[:3]-exact[:3])*1000),ph_group_difference_m=float(abs(ph[3]-exact[3])*1000),ph_group_relative_difference=float(abs(ph[3]-exact[3])/abs(exact[3])),candidate_only_reconstruction_position_error_m=float(np.linalg.norm(recon[:3]-exact[:3])*1000),candidate_only_reconstruction_group_error_m=float(abs(recon[3]-exact[3])*1000),candidate_reconstruction_error_bounds_PH=False)
        records.append(item)
        # Analytic stencil controls verify the endpoint and retained-neighbor operation.
        straight=np.column_stack((ss,np.zeros_like(ss),np.zeros_like(ss),ss));v=PchipInterpolator(ss,straight,axis=0)(20)
        radius=6371.;curve=np.column_stack((radius*np.cos(ss/radius),radius*np.sin(ss/radius),np.zeros_like(ss),ss+ss*ss*1e-7));truth=np.array([radius*np.cos(20/radius),radius*np.sin(20/radius),0,20+400e-7]);cv=PchipInterpolator(ss,curve,axis=0)(20)
        controls.append(dict(job_id=j['job_id'],straight_error_m=float(np.max(abs(v-[20,0,0,20]))*1000),curved_position_error_m=float(np.linalg.norm(cv[:3]-truth[:3])*1000),curved_group_error_m=float(abs(cv[3]-truth[3])*1000)))
    base.write(out/'endpoint_metrics.json',records);base.write(out/'analytic_controls.json',controls)
    flat=[]
    for x in records:
        for kind,c in x['candidates'].items():flat.append(dict(job_id=x['job_id'],level=x['level'],candidate=kind,endpoint_difference_m=c['ph_endpoint_difference_m'],group_difference_m=c['ph_group_difference_m'],group_relative_difference=c['ph_group_relative_difference'],pchip_chord_position_gap_m=x['pchip_chord_position_gap_m'],pchip_chord_group_gap_m=x['pchip_chord_group_gap_m'],stencil_contains_density_root=x['stencil_contains_density_root'],formal_endpoint_passed=False))
    with (out/'endpoint_metrics.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(flat[0]));w.writeheader();w.writerows(flat)
    base.write(out/'conclusion.json',dict(exact20_endpoint_diagnostics_complete=True,case_levels=72,physical_conditions=24,all_real_brackets=True,all_real_following_neighbors=True,state_tangent_used=False,analytic_straight_checks_passed=bool(max(x['straight_error_m'] for x in controls)<1e-7),formal_20km_endpoint_acceptance=False,boundary_convention_identified=False,recommend_p4=False,reason='PH endpoint uses positional/native-group PCHIP. Its method gap and candidate-only reconstruction errors are diagnostics, not independent PH reconstruction bounds. No state derivative is applied across reset.'))
    (out/'endpoint_report.txt').write_text('P3.1-A精确20km端点补充\n72例均使用真实包围节点及下一真实邻居，以PCHIP分别重建位置和原生积分群程；不使用p或输出群折射率作为插值导数。\n两候选独立ODE直接查询s=20km，与PH端点和群程比较；同时列出线性弦/PCHIP差值及候选自身重建误差。\n方法差不是误差上界；候选自身误差不能移用为PH误差上界。存在跨零点或状态重置的模板单独标记。\n完成精确20km诊断，但PH重建尚无独立认证，正式端点验收false，入口机制未解决。\n',encoding='utf-8')
    base.write(out/'sha256.json',{str(p.relative_to(out)):base.common.sha(p) for p in out.iterdir() if p.is_file()})
    print(json.dumps({'records':len(records),'output':str(out)}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);run(p.parse_args().run)
