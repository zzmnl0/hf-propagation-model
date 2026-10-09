"""Read-only P3 comparison with the P2-v2 physical Hermite contract."""
from __future__ import annotations
import os
for _k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[_k]='1'
import sys
sys.dont_write_bytecode=True
import argparse,json,shutil
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicHermiteSpline,PchipInterpolator
import _a4_c1m_p3_experiments as exp
import _a4_c1m_p2_analysis_core_v2 as p2
import _run_a4_c1m_attribution as common

LIMITS={'rms_m':10.,'max_m':50.,'landing_m':50.,'group_rel':1e-5}
STEPS=(1.,.5,.25)

def read_trace(d,job,side):
    if side=='ah':
        meta=exp.read(d/'ah'/(job['job_id']+'.json'))
        if meta['status']!='completed':raise ValueError(meta.get('error','AH failed'))
        path=d/'ah'/(job['job_id']+'.npz')
        if common.sha(path)!=meta['raw_sha256']:raise ValueError('AH hash changed')
        with np.load(path) as raw:r={k:np.array(raw[k]) for k in ('s','xyz','p','group','tangent')}
        r['termination']=meta.get('termination');r['meta']=meta
    else:
        path=d/'pharlap_raw'/(job.get('ph_job_id',job['job_id'])+'.mat');raw=common.mat(path)
        identity=job
        if 'ph_job_id' in job:identity=next(j for j in exp.read(d/'manifest.json')['jobs'] if j['job_id']==job['ph_job_id'])
        if raw['job_sha256']!=identity['job_sha256'] or raw['grid_sha256']!=job['grid_sha256']:raise ValueError('PH identity mismatch')
        st=raw['ray_state_vec'];p=np.column_stack([common.vector(st,k) for k in ('dir_x','dir_y','dir_z')])
        r={'s':common.vector(st,'geometrical_path')/1000.,'xyz':np.column_stack([common.vector(st,k) for k in ('pos_x','pos_y','pos_z')])/1000.,'p':p,'group':common.vector(st,'group_path')/1000.,'tangent':p/np.linalg.norm(p,axis=1)[:,None],'termination':float(np.asarray(raw['ray_data']['ray_label']).reshape(-1)[-1]),'meta':{'mex_sha256':raw['mex_sha256'],'exporter_sha256':raw['exporter_sha256']}}
    if len(r['s'])<2 or np.any(np.diff(r['s'])<=0) or not all(np.isfinite(r[k]).all() for k in ('s','xyz','p','group')):raise ValueError('Nonfinite/nonmonotonic native state')
    if abs(r['s'][0])>1e-8 or abs(r['group'][0])>1e-8:raise ValueError('Nonzero origin path')
    r['raw_path']=str(path.resolve());r['raw_sha256']=common.sha(path)
    return r

def at(trace,q,method='derivative_hermite'):
    if q[-1]>trace['s'][-1]+1e-8:raise ValueError('Early termination: requested path not covered')
    if method=='derivative_hermite':
        xyz=CubicHermiteSpline(trace['s'],trace['xyz'],trace['tangent'],axis=0)(q)
        group=CubicHermiteSpline(trace['s'],trace['group'],1/np.linalg.norm(trace['p'],axis=1))(q)
    else:
        xyz=PchipInterpolator(trace['s'],trace['xyz'],axis=0)(q);group=PchipInterpolator(trace['s'],trace['group'])(q)
    return {'xyz':xyz,'g100':float(group[-1])}

def query(step,end=100.):return np.linspace(0.,end,int(np.ceil(end/step))+1)
def metric(a,b):return p2.errors(a,b)
def controls(values):return p2.controls(values)

def equivalent_reference_field(field):
    coefficients=field['coefficients'];heights=field['h']
    exact=bool(np.all(coefficients[:2]==0) and np.ptp(coefficients[2])==0 and np.array_equal(coefficients[3],coefficients[3,0]+coefficients[2,0]*(heights[:-1]-heights[0])))
    proof={'exact_global_affine':exact,'original_intervals':len(heights)-1,'zero_higher_coefficients':bool(np.all(coefficients[:2]==0)),'slope_spread':float(np.ptp(coefficients[2])),'intercept_max_difference':float(np.max(abs(coefficients[3]-(coefficients[3,0]+coefficients[2,0]*(heights[:-1]-heights[0]))))),'only_reference_knots_changed':True}
    if not exact:return field,proof
    reduced=dict(field,h=heights[[0,-1]],coefficients=coefficients[:,:1].copy())
    if 'ne' in reduced:reduced['ne']=np.asarray(reduced['ne'])[[0,-1]]
    return reduced,proof

def compare_case(d,jobs,out):
    jobs=sorted(jobs,key=lambda j:j['level']);first=jobs[0]
    row={'case_id':first['job_id'].rsplit('_t',1)[0],'role':first['role'],'scope':first['scope'],'jobs':[j['job_id'] for j in jobs],'strict_passed':False,'limits':LIMITS}
    if len(jobs)!=3:raise ValueError('Missing integration level')
    traces={side:[read_trace(d,j,side) for j in jobs] for side in ('ah','ph')}
    row['levels']={side:[{'level':j['level'],'raw_path':t['raw_path'],'raw_sha256':t['raw_sha256'],'termination':t['termination'],'native_s_end_km':float(t['s'][-1]),'native_group_end_km':float(t['group'][-1])} for j,t in zip(jobs,traces[side])] for side in traces}
    if first['scope']!='first_100_km':
        raise ValueError('Full-path reference unavailable for rejected signed-density Hermite medium; no truncated acceptance')
    q=query(.25);v={side:[at(t,q) for t in ts] for side,ts in traces.items()}
    actual=metric(v['ah'][-1],v['ph'][-1]);integration={side:controls(values) for side,values in v.items()}
    initial=max(float(np.linalg.norm(a['xyz'][0]-b['xyz'][0])*1000.) for a,b in zip(traces['ah'],traces['ph']))
    pdelta=max(float(np.linalg.norm(a['p'][0]-b['p'][0])) for a,b in zip(traces['ah'],traces['ph']))
    sampled=[metric(at(traces['ah'][-1],query(h)),at(traces['ph'][-1],query(h))) for h in STEPS]
    sampling={k:abs(sampled[-1][k]-sampled[-2][k]) for k in LIMITS}
    sampling12={k:abs(sampled[1][k]-sampled[0][k]) for k in LIMITS}
    sampling_monotonic={k:sampling[k]<=sampling12[k] for k in LIMITS}
    ode=p2.load_ode(exp.COEFFICIENT);field,field_reduction=equivalent_reference_field(ode.medium(first,False))
    max_s=max(float(p2.prefix(t)['s'][-1]) for side in traces for t in traces[side])
    ref=p2.RadialReference(traces['ph'][-1],field,ode,exp.COEFFICIENT,first['frequency_MHz'],5,max_s)
    reference_levels=[]
    for level in (3,4):
        candidate=p2.RadialReference(traces['ph'][-1],field,ode,exp.COEFFICIENT,first['frequency_MHz'],level,max_s)
        reference_levels.append(candidate.query(q))
        common.write_json(out/(row['case_id']+f'_quadrature_l{level}.json'),{'estimates':list(candidate.uncertainty_records.values()),'warnings':candidate.warnings,'interval_certified':False})
    reference_levels.append(ref.query(q))
    reference_convergence=controls(reference_levels)
    proof=p2.field_bounds(field,ref.alpha,ref.E,max_s)
    reconstruction={side:p2.reconstruction_bound(traces[side][-1],'derivative_hermite',ref,proof,side) for side in traces}
    bound={k:sum(reconstruction[side]['bounds'][k] for side in traces) for k in LIMITS}
    budget={k:3*(integration['ah']['delta23'][k]+integration['ph']['delta23'][k]+sampling[k]+bound[k]) for k in LIMITS}
    topology={'matched':True,'eligible':True,'expected':'first_100_km','coverage_complete':all(t['s'][-1]>=100.-1e-8 for ts in traces.values() for t in ts) and all(abs(t['s'][-1]-100.)<=1e-8 and t['termination']=='fixed_geometric_path_completed' for t in traces['ah'])}
    reference_valid=reference_convergence['passed'] and not ref.warnings
    checks={'thresholds':all(actual[k]<=LIMITS[k] for k in LIMITS),'budget':all(actual[k]<=budget[k] for k in LIMITS),'integration':all(x['passed'] for x in integration.values()),'sampling':all(sampling[k]<=.2*LIMITS[k] and sampling_monotonic[k] for k in LIMITS),'reconstruction':reference_valid and all(bound[k]<=.2*LIMITS[k] for k in LIMITS),'initial_state':initial<=1e-3 and pdelta<=1e-10,'scope':topology['coverage_complete']}
    reference=ref.query(q);delta=v['ah'][-1]['xyz']-v['ph'][-1]['xyz']
    common.write_json(out/(row['case_id']+'_quadrature.json'),{'estimates':list(ref.uncertainty_records.values()),'warnings':ref.warnings,'interval_certified':False})
    p2.write_csv(out/(row['case_id']+'_points.csv'),[{'s_km':float(s),'dx_m':float(e[0]*1000),'dy_m':float(e[1]*1000),'dz_m':float(e[2]*1000),'error_m':float(np.linalg.norm(e)*1000)} for s,e in zip(q,delta)])
    row.update(actual=actual,integration=integration,resampling={'steps_km':STEPS,'metrics':sampled,'delta12':sampling12,'monotonic':sampling_monotonic,'final_bound':sampling},reconstruction={'validated':reference_valid,'sides':reconstruction,'bound':bound},budget=budget,topology=topology,checks=checks,initial_position_m=initial,initial_p_difference=pdelta,strict_passed=all(checks.values()),independent_reference={'AH':metric(v['ah'][-1],reference),'PH':metric(v['ph'][-1],reference),'convergence':reference_convergence,'field_reduction':field_reduction,'quadrature_level':5,'interval_certified':False})
    return row

def analyze(run):
    d=run/'experiments';out=run/'analysis/p3_v2'
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True);shutil.copy2(Path(__file__),out/'analysis_source.py')
    manifest=exp.read(d/'manifest.json');cases={}
    for j in manifest['jobs']:
        if j['role'].startswith('R5_'):continue
        cases.setdefault(j['job_id'].rsplit('_t',1)[0],[]).append(j)
    results=[]
    for key,jobs in cases.items():
        try:row=compare_case(d,jobs,out)
        except Exception as e:row={'case_id':key,'role':jobs[0]['role'],'scope':jobs[0]['scope'],'jobs':[j['job_id'] for j in jobs],'strict_passed':False,'error':repr(e)}
        exp.write(out/(key+'.json'),row);results.append(row)
        print(json.dumps({'case':key,'passed':row['strict_passed'],'error':row.get('error')}),flush=True)
    formal=[r for r in results if r['role']!='R5_grid']
    exp.write(out/'summary.json',{'schema':'c1m-p3-comparison-v1','case_count':len(formal),'expected_case_count':24,'cases':results,'formal_passed_count':sum(r['strict_passed'] for r in formal),'p3_scope_passed':len(formal)==24 and all(r['strict_passed'] for r in formal),'global_passed':False,'recommend_p4':False,'analysis_sha256':common.sha(Path(__file__))})

def diagnostic_scans(run):
    d=run/'experiments';out=run/'analysis/p3_r5_v1'
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True);manifest=exp.read(d/'manifest.json');rows=[];cache={}
    for j in manifest['jobs']:
        if not j['role'].startswith('R5_'):continue
        try:
            a=read_trace(d,j,'ah');b=read_trace(d,j,'ph');cache[j['job_id']]=a
            if j['scope']!='first_100_km':raise ValueError('Full trajectory input rejected; no local substitute')
            rows.append({'job_id':j['job_id'],'role':j['role'],'scope':j['scope'],'cross':metric(at(a,query(.25)),at(b,query(.25))),'status':'computed','strict_acceptance':False})
        except Exception as e:rows.append({'job_id':j['job_id'],'role':j['role'],'status':'failed','error':repr(e),'strict_acceptance':False})
    fd=[]
    for base in [j for j in manifest['jobs'] if j['role'] in ('R1','R2') and j['level']==3]:
        probes=[j for j in manifest['jobs'] if j.get('ph_job_id')==base['job_id']]
        if len(probes)!=2:continue
        try:
            seq=sorted(probes+[base],key=lambda j:-j['spatial_fd_base_km'])
            traces=[read_trace(d,j,'ah') for j in seq]
            c=controls([at(t,query(.25)) for t in traces])
            fd.append({'base_job_id':base['job_id'],'FD_km':[j['spatial_fd_base_km'] for j in seq],'controls':c,'converged':c['passed']})
        except Exception as e:fd.append({'base_job_id':base['job_id'],'converged':False,'error':repr(e)})
    grids=[];groups={}
    for j in manifest['jobs']:
        if j['role']=='R5_grid' and j['level']==3:
            key=(j['case_id'].split('_d')[0],j['elevation_deg']);groups.setdefault(key,[]).append(j)
    for key,jobs in groups.items():
        try:
            seq=sorted(jobs,key=lambda j:-float(np.diff(common.mat(Path(j['grid_mat']))['height_km'])[0]))
            values=[at(read_trace(d,j,'ah'),query(.25)) for j in seq]
            grids.append({'case':list(key),'height_steps_km':[float(np.diff(common.mat(Path(j['grid_mat']))['height_km'])[0]) for j in seq],'successive_changes':[metric(a,b) for a,b in zip(values[:-1],values[1:])],'numerical_uncertainty_budget_contribution':False})
        except Exception as e:grids.append({'case':list(key),'error':repr(e)})
    exp.write(out/'summary.json',{'role':'single_factor_diagnostics','jobs':rows,'fd':fd,'grid':grids,'included_in_reconstruction_budget':False,'source_sha256':common.sha(Path(__file__))})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--diagnostic-scans',action='store_true');a=p.parse_args();diagnostic_scans(a.run.resolve()) if a.diagnostic_scans else analyze(a.run.resolve())
