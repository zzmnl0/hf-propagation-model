"""Independent physical gates and derivative-aware R3 certification from MP records."""
from __future__ import annotations
import argparse,shutil
from pathlib import Path
import numpy as np
import mpmath as mp
from scipy.interpolate import CubicHermiteSpline
import _a4_c1m_p31_certify as cert
from _a4_c1m_p31_reference import AffineReference
comp=cert.comp


def run_case(parent,previous,out,job):
    case=job['job_id'].rsplit('_t',1)[0];src=previous/case;dest=out/case;dest.mkdir()
    original=comp.exp.read(src/'result.json');old=comp.exp.read(parent/'analysis/p3_refined_v1'/(case+'.json'))
    with np.load(src/'frozen_field.npz') as raw:field={k:np.array(raw[k]) for k in raw.files}
    with np.load(src/'reference_dps160_n32.npz') as raw:truth={k:np.array(raw[k]) for k in raw.files}
    d=parent/'experiments';ph=comp.read_trace(d,job,'ph');ap=Path(old['levels']['ah'][-1]['raw_path']);meta=comp.exp.read(ap.with_suffix('.json'))
    if comp.common.sha(ap)!=meta['raw_sha256']:raise ValueError('AH checksum mismatch')
    with np.load(ap) as raw:ah={k:np.array(raw[k]) for k in ('s','xyz','p','group','tangent','ne','llh')}
    traces={s:comp.p2.prefix(t) for s,t in {'ah':ah,'ph':ph}.items()}
    extra=np.unique(np.concatenate([t['s'][t['s']>100.] for t in traces.values()]))
    supplements={};proof=None
    for dps,order in ((80,32),(120,32),(160,16),(160,24),(160,32)):
        r=AffineReference(field,ph['xyz'][0],ph['p'][0],comp.exp.COEFFICIENT,job['frequency_MHz'],dps,order);proof=r.proof
        v=r.query(extra);supplements[dps,order]=v
        comp.common.write_json(dest/f'support_dps{dps}_n{order}.json',{'s':extra.tolist(),'decimal':v['decimal']})
    changes=[cert.decimal_difference(supplements[a]['decimal'],supplements[b]['decimal']) for a,b in (((80,32),(120,32)),((120,32),(160,32)),((160,16),(160,24)),((160,24),(160,32)))]
    supplement_stable=all(float(x['position_max_m'])<=1e-12 and float(x['group_max_km'])<=1e-15 for x in changes)
    truth={k:np.concatenate([truth[k],supplements[160,32][k]],axis=0) for k in ('xyz','p','group')} | {'s':np.r_[truth['s'],extra]}
    order=np.argsort(truth['s']);truth={k:v[order] for k,v in truth.items()}
    q=comp.query(.0625);idx=np.searchsorted(truth['s'],q);formal={k:truth[k][idx] for k in ('xyz','p','group')}
    alpha=comp.exp.COEFFICIENT/(job['frequency_MHz']*1e6)**2;gradient=float(abs(field['coefficients'][2,0]));n2min=float(min(np.sum(truth['p']**2,axis=1)))
    M2=alpha*gradient/(2*n2min);G2=alpha*gradient/(2*n2min**1.5)
    with mp.workdps(170):
        v=supplements[160,32]
        converted=[[mp.nstr(cert.exact(x),170) for x in (*pos,g)] for pos,g in zip(v['xyz'],v['group'])]
    support_rounding=cert.decimal_difference(v['decimal'],converted)
    support_pos=max(float(x['position_max_m']) for x in changes[-1:])+float(support_rounding['position_max_m'])
    support_group=max(float(x['group_max_km']) for x in changes[-1:])+float(support_rounding['group_max_km'])
    refbudget=dict(original['reference_budget'])
    for key in ('rms_m','max_m','landing_m'):refbudget[key]=max(refbudget[key],support_pos)
    refbudget['group_rel']=max(refbudget['group_rel'],support_group/formal['group'][-1])
    bounds={};envelopes={};native_rows=[]
    for side,t in traces.items():
        ix=np.searchsorted(truth['s'],t['s']);x=truth['xyz'][ix];p=truth['p'][ix];g=truth['group'][ix];tangent=p/np.linalg.norm(p,axis=1)[:,None]
        spline=CubicHermiteSpline(t['s'],x,tangent);gs=CubicHermiteSpline(t['s'],g,1/np.linalg.norm(p,axis=1))
        error=np.linalg.norm(spline(q)-formal['xyz'],axis=1);deriv=np.linalg.norm(spline.derivative()(q)-formal['p']/np.linalg.norm(formal['p'],axis=1)[:,None],axis=1)
        fill=.0625/2*max(deriv)+.0625**2/8*(comp.p2.polynomial_second_max(spline)+M2)
        ge=np.abs(gs(q)-formal['group']);gd=np.abs(gs.derivative()(q)-1/np.linalg.norm(formal['p'],axis=1));gf=.0625/2*max(gd)+.0625**2/8*(comp.p2.polynomial_second_max(gs)+G2)
        bounds[side]={'rms_m':float((max(error)+fill)*1000+refbudget['rms_m']),'max_m':float((max(error)+fill)*1000+refbudget['max_m']),'landing_m':float(error[-1]*1000+refbudget['landing_m']),'group_rel':float((ge[-1]+refbudget['group_rel']*formal['group'][-1])/formal['group'][-1])}
        dx=np.linalg.norm(t['xyz']-x,axis=1);dt=np.linalg.norm(t['p']/np.linalg.norm(t['p'],axis=1)[:,None]-tangent,axis=1);dg=np.abs(t['group']-g);dgd=np.abs(1/np.linalg.norm(t['p'],axis=1)-1/np.linalg.norm(p,axis=1))
        eligible=np.flatnonzero(t['s'][:-1]<100.);h=np.diff(t['s'])[eligible]
        pos=float(max(np.maximum(dx[eligible],dx[eligible+1])+h/4*np.maximum(dt[eligible],dt[eligible+1]))*1000)
        group=float(max(np.maximum(dg[eligible],dg[eligible+1])+h/4*np.maximum(dgd[eligible],dgd[eligible+1]))/formal['group'][-1])
        envelopes[side]={'rms_m':pos+refbudget['rms_m'],'max_m':pos+refbudget['max_m'],'landing_m':pos+refbudget['landing_m'],'group_rel':group+refbudget['group_rel']}
        native_rows.extend({'side':side,'s_km':float(s),'position_m':float(e*1000),'direction_error':float(a),'group_error_km':float(b),'group_derivative_error':float(c),'support_only':bool(s>100)} for s,e,a,b,c in zip(t['s'],dx,dt,dg,dgd))
        bounds[side]['continuous_group_max_rel']=float((max(ge)+gf)/formal['group'][-1])
    grid=comp.common.mat(Path(job['grid_mat']));phraw=comp.common.mat(Path(ph['raw_path']));path=phraw['ray_path_data']
    initialxyz=6471.*np.array([np.cos(np.pi/6)*np.cos(2*np.pi/3),np.cos(np.pi/6)*np.sin(2*np.pi/3),.5]);north=np.array([-.5*np.cos(2*np.pi/3),-.5*np.sin(2*np.pi/3),np.cos(np.pi/6)])
    ne100=field['coefficients'][3,0]+field['coefficients'][2,0]*(100-field['h'][0]);e=np.deg2rad(job['elevation_deg']);expectedp=np.sqrt(1-alpha*ne100)*(np.sin(e)*initialxyz/6471.+np.cos(e)*north)
    ah_expected=field['coefficients'][3,0]+field['coefficients'][2,0]*(ah['llh'][:,2]-field['h'][0])
    ph_h=comp.common.vector(path,'height');ph_expected=field['coefficients'][3,0]+field['coefficients'][2,0]*(ph_h-field['h'][0]);ph_ne=comp.common.vector(path,'electron_density')*1e6
    local=ph['s']<=100.;phmask=ph_h<=max(np.linalg.norm(ph['xyz'][local],axis=1)-6371.)+1e-6
    m=meta['metadata'];physical={'same_frequency':bool(m['frequency_MHz']==job['frequency_MHz'] and np.all(np.asarray(path['frequency'])==job['frequency_MHz'])),'same_coefficient':bool(m['plasma_coefficient']==comp.exp.COEFFICIENT and job['launch_k2']==comp.exp.COEFFICIENT),'sphere_6371':bool(float(grid['earth_radius_m'])==6371000 and m['earth_radius_km']==6371.),'input_B0_collisionless':bool(all(np.all(np.asarray(grid[k])==0) for k in ('Bx_T','By_T','Bz_T','collision_freq_Hz'))),'AH_B0_collisionless':bool(m['gyrofrequency_MHz']==0 and m['collision_frequency_s^-1']==0),'PH_B0_collisionless':bool(all(np.all(comp.common.vector(path,k)==0) for k in ('geomag_x','geomag_y','geomag_z','collision_frequency'))),'AH_density':bool(np.max(np.abs(ah['ne']-ah_expected))<=.1),'PH_density_local':bool(np.max(np.abs(ph_ne[phmask]-ph_expected[phmask]))<=.1),'initial_position_each':bool(all(np.linalg.norm(t['xyz'][0]-initialxyz)*1000<=1e-5 for t in traces.values())),'initial_momentum_each':bool(all(np.linalg.norm(t['p'][0]-expectedp)<=1e-10 for t in traces.values())),'coverage_each':bool(all(t['s'][-1]>=100 and t['s'][0]==0 for t in traces.values())),'AH_termination':meta['termination']=='fixed_geometric_path_completed','PH_local_scope':bool(ph['s'][-1]>=100 and old['scope']=='first_100_km')}
    budgets={k:3*(sum(envelopes[s][k]+bounds[s][k] for s in traces)+original['sampling'][k]) for k in cert.LIMITS}
    checks=dict(original['checks']);checks.update(common_physics=all(physical.values()),reference_support=supplement_stable,reconstruction=all(cert.gate(b) and b['continuous_group_max_rel']<=2e-6 for b in bounds.values()),budget=all(original['actual'][k]<=budgets[k] for k in cert.LIMITS),native_envelope=all(cert.gate(v) for v in envelopes.values()))
    result=dict(original,reference_certified_passed=bool(all(checks.values())),checks=checks,physical_checks=physical,reconstruction=bounds,side_envelopes=envelopes,budget=budgets,supplement_changes=changes,support_rounding=support_rounding,reference_budget=refbudget,reference_affine_proof=proof,propagation_bound='Each Hermite segment: max endpoint delta position + h/4 max endpoint delta tangent; group analogous. No propagation term in reconstruction.',reconstruction_proof={'validation_step_km':.0625,'position_M2':M2,'group_M2':G2,'formula':'sample maximum + step/2 max derivative error + step^2/8(max spline second derivative + reference M2)'})
    comp.p2.write_csv(dest/'native_errors.csv',native_rows);comp.common.write_json(dest/'result.json',result);return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--parent',type=Path,required=True);p.add_argument('--run',type=Path,required=True);a=p.parse_args();out=a.run/'reference_r3_v4';out.mkdir(exist_ok=True)
    source=out/'source';source.mkdir(exist_ok=True)
    for name in ('_a4_c1m_p31_reference.py','_a4_c1m_p31_certify.py','_a4_c1m_p31_certify_finalize.py'):
        src=Path(__file__).with_name(name);dest=source/name
        if dest.exists() and comp.common.sha(src)!=comp.common.sha(dest):raise ValueError('Source mutation requires new analysis version')
        if not dest.exists():shutil.copy2(src,dest)
    rows=[]
    for job in comp.exp.read(a.parent/'experiments/manifest.json')['jobs']:
        if job['role']!='R3' or job['level']!=3:continue
        case=job['job_id'].rsplit('_t',1)[0]
        if not (a.run/'reference_r3_v2'/case/'result.json').exists():continue
        path=out/case/'result.json';r=comp.exp.read(path) if path.exists() else run_case(a.parent,a.run/'reference_r3_v2',out,job);rows.append(r);print(case,r['reference_certified_passed'],flush=True)
    comp.common.write_json(out/'summary.json',{'case_count':len(rows),'reference_certified_count':sum(r['reference_certified_passed'] for r in rows),'all_16_certified':len(rows)==16 and all(r['reference_certified_passed'] for r in rows),'cases':rows})
if __name__=='__main__':main()
