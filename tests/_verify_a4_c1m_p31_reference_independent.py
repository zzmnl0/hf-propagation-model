"""Independent cheap re-evaluation from frozen native states and MP records.

No certification helper or production ray solver is imported; no MP integration reruns.
"""
import argparse
import hashlib
import json
import shutil
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.interpolate import CubicHermiteSpline, PPoly
from scipy.io import loadmat


LIMITS={'rms_m':10.,'max_m':50.,'landing_m':50.,'group_rel':1e-5}
C=80.61638494312905


def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False),encoding='utf-8')
def exact(x):
    n,d=float(x).as_integer_ratio();return mp.mpf(n)/d
def metric(x,g,y,gy):
    error=np.linalg.norm(x-y,axis=1)*1000
    return dict(rms_m=float(np.sqrt(np.mean(error**2))),max_m=float(max(error)),
                landing_m=float(error[-1]),group_rel=float(abs(g[-1]-gy[-1])/abs(gy[-1])))
def gate(values,ratio=.2):return all(np.isfinite(values[k]) and values[k]<=ratio*LIMITS[k] for k in LIMITS)
def near(a,b,k):return bool(np.isclose(a,b,rtol=5e-8,atol=1e-16 if 'rel' in k else 5e-10))
def native_valid(t):
    return bool(len(t['s'])>=2 and t['s'][0]==0 and t['s'][-1]>=100 and np.all(np.diff(t['s'])>0)
                and all(np.isfinite(t[k]).all() for k in t))
def initial_valid(t,xyz,p):
    return bool(np.linalg.norm(t['xyz'][0]-xyz)*1000<=1e-5 and np.linalg.norm(t['p'][0]-p)<=1e-10)
def budget_valid(actual,budget):
    return gate(actual,1.) and all(np.isfinite(budget[k]) and 0<=actual[k]<=budget[k] for k in LIMITS)


def support_momenta(records,field,x0,p0,f):
    # Reconstruct momentum algebraically from stored MP positions and the frozen invariant.
    with mp.workdps(170):
        x0=mp.matrix([exact(v) for v in x0]);p0=mp.matrix([exact(v) for v in p0])
        r0=mp.norm(x0);er0=x0/r0;radial=mp.fdot(p0,er0);pt=p0-radial*er0;et0=pt/mp.norm(pt)
        J=r0*mp.norm(pt);a=exact(C)/(exact(f)*1000000)**2
        h0=exact(field['h'][0]);g=exact(field['coefficients'][2,0]);n0=exact(field['coefficients'][3,0])
        E=mp.fdot(p0,p0)+a*(n0+g*(r0-6371-h0));rows=[]
        for record in records:
            x=mp.matrix([mp.mpf(v) for v in record[:3]]);r=mp.norm(x);er=x/r
            sin=mp.fdot(er,et0);cos=mp.fdot(er,er0);et=-sin*er0+cos*et0
            n2=E-a*(n0+g*(r-6371-h0));p=mp.sqrt(n2-(J/r)**2)*er+J/r*et
            rows.append([float(v) for v in p])
        return np.asarray(rows)


def native_traces(parent,job,old):
    pp=parent/'experiments/pharlap_raw'/(job['job_id']+'.mat');raw=loadmat(pp,simplify_cells=True);st=raw['ray_state_vec']
    p=np.column_stack([st[k] for k in ('dir_x','dir_y','dir_z')]);ph={'s':np.asarray(st['geometrical_path'])/1000.,
        'xyz':np.column_stack([st[k] for k in ('pos_x','pos_y','pos_z')])/1000.,'p':p,'group':np.asarray(st['group_path'])/1000.,
        'tangent':p/np.linalg.norm(p,axis=1)[:,None]}
    ap=Path(old['levels']['ah'][-1]['raw_path']);ah=dict(np.load(ap));meta=read(ap.with_suffix('.json'))
    full={'ah':ah,'ph':ph};prefix={}
    for side,t in full.items():
        r=np.searchsorted(t['s'],100.)
        if r>=len(t['s']):raise ValueError('Missing actual 100 km right support')
        prefix[side]={k:t[k][:min(r+2,len(t['s']))] for k in ('s','xyz','p','group','tangent')}
    return full,prefix,raw,meta,pp,ap


def curvature_max(spline,formal_only=True):
    eligible=spline.x[:-1]<=100 if formal_only else np.ones(len(spline.x)-1,dtype=bool)
    width=(np.minimum(spline.x[1:],100)-spline.x[:-1] if formal_only else np.diff(spline.x))[eligible]
    left=2*spline.c[1,eligible]
    right=left+6*spline.c[0,eligible]*width.reshape((-1,)+(1,)*(left.ndim-1))
    if left.ndim==1:return float(max(np.max(abs(left)),np.max(abs(right))))
    return float(max(np.max(np.linalg.norm(left,axis=1)),np.max(np.linalg.norm(right,axis=1))))


def verify_case(parent,run,analysis,job):
    case=job['job_id'].rsplit('_t',1)[0];result=read(analysis/case/'result.json')
    old=read(parent/'analysis/p3_refined_v1'/(case+'.json'));previous=run/'reference_r3_v2'/case
    field=dict(np.load(previous/'frozen_field.npz'));truth=dict(np.load(previous/'reference_dps160_n32.npz'))
    full,traces,raw,meta,pp,ap=native_traces(parent,job,old)
    support=read(analysis/case/'support_dps160_n32.json');sx=np.array([[float(v) for v in row[:3]] for row in support['decimal']]);sg=np.array([float(row[3]) for row in support['decimal']])
    sp=support_momenta(support['decimal'],field,full['ph']['xyz'][0],full['ph']['p'][0],job['frequency_MHz'])
    for key,val in (('s',np.asarray(support['s'])),('xyz',sx),('p',sp),('group',sg)):truth[key]=np.concatenate([truth[key],val])
    ix=np.argsort(truth['s']);truth={k:v[ix] for k,v in truth.items()}
    q=np.linspace(0,100,1601);qi=np.searchsorted(truth['s'],q)
    if not np.array_equal(truth['s'][qi],q):raise ValueError('Reference samples do not cover formal query nodes')
    formal={k:truth[k][qi] for k in ('xyz','p','group')};rv=result['reference_budget']
    a=C/(job['frequency_MHz']*1e6)**2;g=abs(field['coefficients'][2,0]);n2min=np.min(np.sum(truth['p']**2,axis=1))
    M2=a*g/(2*n2min);G2=a*g/(2*n2min**1.5)
    sides={};bounds={};envelopes={};vectors={};checks={};continuous_group={};curvature_domains={}
    checks['affine_exact']=bool(np.all(field['coefficients'][:2]==0) and np.all(field['coefficients'][2]==field['coefficients'][2,0]) and
        all(Fraction(float(field['coefficients'][3,i]))==Fraction(float(field['coefficients'][3,0]))+Fraction(float(field['coefficients'][2,0]))*(Fraction(float(h))-Fraction(float(field['h'][0]))) for i,h in enumerate(field['h'][:-1])))
    radius=np.linalg.norm(truth['xyz'],axis=1)
    checks['continuous_curvature_assumptions']=bool(np.all(np.diff(radius)>=-1e-9) and n2min>0 and np.isfinite(M2) and np.isfinite(G2))
    checks['curvature_M2_record_matches']=bool(near(M2,result['reconstruction_proof']['position_M2'],'M2') and near(G2,result['reconstruction_proof']['group_M2'],'M2'))
    for side,t in traces.items():
        ids=np.searchsorted(truth['s'],t['s'])
        if not np.array_equal(truth['s'][ids],t['s']):raise ValueError('Missing true native support')
        reference={k:truth[k][ids] for k in ('xyz','p','group')}
        direction=t['p']/np.linalg.norm(t['p'],axis=1)[:,None]
        tangent=reference['p']/np.linalg.norm(reference['p'],axis=1)[:,None]
        checks[side+'_stored_tangent_is_physical']=bool(np.max(np.linalg.norm(t['tangent']-direction,axis=1))<=1e-14)
        checks[side+'_native_valid']=native_valid(t)
        xy=CubicHermiteSpline(t['s'],t['xyz'],t['tangent'])(q)
        group=CubicHermiteSpline(t['s'],t['group'],1/np.linalg.norm(t['p'],axis=1))(q)
        vectors[side]=(xy,group);sides[side]=metric(xy,group,formal['xyz'],formal['group'])
        ideal=CubicHermiteSpline(t['s'],reference['xyz'],tangent)
        ig=CubicHermiteSpline(t['s'],reference['group'],1/np.linalg.norm(reference['p'],axis=1))
        curvature_domains[side]={'formal_position_second_max':curvature_max(ideal),
                                'all_support_position_second_max':curvature_max(ideal,False),
                                'formal_group_second_max':curvature_max(ig),
                                'all_support_group_second_max':curvature_max(ig,False),
                                'support_s_end_km':float(t['s'][-1])}
        err=np.linalg.norm(ideal(q)-formal['xyz'],axis=1)
        de=np.linalg.norm(ideal.derivative()(q)-formal['p']/np.linalg.norm(formal['p'],axis=1)[:,None],axis=1)
        fill=.0625/2*max(de)+.0625**2/8*(curvature_max(ideal)+M2)
        ge=abs(ig(q)-formal['group']);gd=abs(ig.derivative()(q)-1/np.linalg.norm(formal['p'],axis=1))
        gf=.0625/2*max(gd)+.0625**2/8*(curvature_max(ig)+G2)
        bounds[side]={'rms_m':float((max(err)+fill)*1000+rv['rms_m']),'max_m':float((max(err)+fill)*1000+rv['max_m']),
                      'landing_m':float(err[-1]*1000+rv['landing_m']),'group_rel':float(ge[-1]/formal['group'][-1]+rv['group_rel'])}
        continuous_group[side]=float((max(ge)+gf)/formal['group'][-1]+rv['group_rel'])
        checks[side+'_continuous_group_with_reference']=continuous_group[side]<=2e-6
        dx=np.linalg.norm(t['xyz']-reference['xyz'],axis=1);dv=np.linalg.norm(direction-tangent,axis=1)
        dg=abs(t['group']-reference['group']);dgd=abs(1/np.linalg.norm(t['p'],axis=1)-1/np.linalg.norm(reference['p'],axis=1))
        eligible=np.flatnonzero(t['s'][:-1]<100);hs=np.diff(t['s'])[eligible]
        position=float(max(np.maximum(dx[eligible],dx[eligible+1])+hs/4*np.maximum(dv[eligible],dv[eligible+1]))*1000)
        group_error=float(max(np.maximum(dg[eligible],dg[eligible+1])+hs/4*np.maximum(dgd[eligible],dgd[eligible+1]))/formal['group'][-1])
        envelopes[side]={k:(group_error if k=='group_rel' else position)+rv[k] for k in LIMITS}
    actual=metric(*vectors['ah'],*vectors['ph'])
    sampled=[]
    for step in old['resampling']['steps_km']:
        query=np.linspace(0,100,int(np.ceil(100/step))+1);values={}
        for side,t in full.items():
            values[side]=(CubicHermiteSpline(t['s'],t['xyz'],t['tangent'])(query),
                          CubicHermiteSpline(t['s'],t['group'],1/np.linalg.norm(t['p'],axis=1))(query))
        sampled.append(metric(*values['ah'],*values['ph']))
    sampling={k:abs(sampled[-1][k]-sampled[-2][k]) for k in LIMITS}
    sampling_monotonic={k:sampling[k]<=abs(sampled[-2][k]-sampled[-3][k]) for k in LIMITS}
    checks['sampling_reproduced']=all(near(sampling[k],result['sampling'][k],k) for k in LIMITS)
    checks['sampling_under_limit']=gate(sampling)
    budget={k:3*(sum(envelopes[s][k]+bounds[s][k] for s in traces)+sampling[k]) for k in LIMITS}
    checks['sources_unchanged']=all(sha(item['path'])==item['sha256'] for item in result['sources'])
    checks['raw_identity']=bool(raw['job_sha256']==job['job_sha256'] and raw['grid_sha256']==job['grid_sha256'] and sha(ap)==meta['raw_sha256'])
    checks['actual_reproduced']=all(near(actual[k],result['actual'][k],k) for k in LIMITS)
    checks['three_way_reproduced']=all(near(sides[s][k],result['formal_errors'][s][k],k) for s in sides for k in LIMITS)
    checks['budget_reproduced']=all(near(budget[k],result['budget'][k],k) for k in LIMITS)
    checks['reconstruction_reproduced']=all(near(bounds[s][k],result['reconstruction'][s][k],k) for s in sides for k in LIMITS)
    checks['propagation_reproduced']=all(near(envelopes[s][k],result['side_envelopes'][s][k],k) for s in sides for k in LIMITS)
    checks['budgets_and_limits']=budget_valid(actual,budget) and all(gate(b) for b in list(bounds.values())+list(envelopes.values())+list(sides.values()))
    xyz0=6471*np.array([-np.sqrt(3)/4,.75,.5]);north=np.array([.25,-np.sqrt(3)/4,np.sqrt(3)/2])
    ne0=field['coefficients'][3,0]+field['coefficients'][2,0]*(100-field['h'][0]);el=np.deg2rad(job['elevation_deg'])
    p0=np.sqrt(1-a*ne0)*(np.sin(el)*xyz0/6471+np.cos(el)*north)
    checks['initial_position_and_wave_normal']=all(initial_valid(t,xyz0,p0) for t in traces.values())
    checks['initial_dispersion']=all(abs(t['p'][0]@t['p'][0]-(1-a*ne0))<=1e-10 for t in traces.values())
    grid=loadmat(job['grid_mat'],simplify_cells=True);path=raw['ray_path_data'];m=meta['metadata']
    local=full['ph']['s']<=100
    phheight=np.asarray(path['height']);phdensity=np.asarray(path['electron_density'])*1e6
    expected_ne=lambda height:field['coefficients'][3,0]+field['coefficients'][2,0]*(height-field['h'][0])
    checks['frequency_and_constant']=bool(m['frequency_MHz']==job['frequency_MHz'] and np.all(np.asarray(path['frequency'])==job['frequency_MHz']) and m['plasma_coefficient']==C and job['launch_k2']==C)
    checks['sphere']=bool(grid['earth_radius_m']==6371000 and m['earth_radius_km']==6371)
    checks['B0_collisionless']=bool(all(np.all(np.asarray(grid[k])==0) for k in ('Bx_T','By_T','Bz_T','collision_freq_Hz')) and m['gyrofrequency_MHz']==0 and m['collision_frequency_s^-1']==0 and all(np.all(np.asarray(path[k])==0) for k in ('geomag_x','geomag_y','geomag_z','collision_frequency')))
    checks['density_actual_local']=bool(np.max(abs(full['ah']['ne']-expected_ne(full['ah']['llh'][:,2])))<=.1 and np.max(abs(phdensity[local]-expected_ne(phheight[local])))<=.1)
    checks['termination_scope']=bool(meta['termination']=='fixed_geometric_path_completed' and abs(full['ah']['s'][-1]-100)<1e-8 and full['ph']['s'][-1]>=100 and old['scope']=='first_100_km')
    checks['physical_flags_all_true']=all(result['physical_checks'].values())
    checks['old_strict_status_unchanged']=result['strict_convergence_passed']==old['strict_passed']
    return {'case_id':case,'checks':checks,'passed':all(checks.values()),'actual':actual,'formal_errors':sides,
            'budget':budget,'reconstruction':bounds,'side_envelopes':envelopes,
            'continuous_group_reconstruction_bound_with_reference':continuous_group,
            'curvature_domain_diagnostic':curvature_domains,
            'sampling':sampling,'sampling_monotonic_diagnostic':sampling_monotonic,
            'mathematical_assumptions':'Exact affine radial medium; ascending radius and positive n squared. |x second derivative| <= alpha|g|/(2 n_min^2); |G second derivative| <= alpha|g|/(2 n_min^3). Cubic spline second derivative is linear, so interval endpoint norms bound it. MP quadrature remainder remains a numerical estimate.',
            'recorded_certification':result['reference_certified_passed'],'result_sha256':sha(analysis/case/'result.json')}


def negative_controls(parent,run,job):
    case=job['job_id'].rsplit('_t',1)[0];old=read(parent/'analysis/p3_refined_v1'/(case+'.json'))
    _,traces,_,_,_,_=native_traces(parent,job,old)
    field=dict(np.load(run/'reference_r3_v2'/case/'frozen_field.npz'))
    xyz=6471*np.array([-np.sqrt(3)/4,.75,.5]);north=np.array([.25,-np.sqrt(3)/4,np.sqrt(3)/2])
    ne=field['coefficients'][3,0]+field['coefficients'][2,0]*(100-field['h'][0]);a=C/(job['frequency_MHz']*1e6)**2
    elevation=np.deg2rad(job['elevation_deg']);p=np.sqrt(1-a*ne)*(np.sin(elevation)*xyz/6471+np.cos(elevation)*north)
    source=traces['ah'];wrong_position={k:v.copy() for k,v in source.items()};wrong_position['xyz'][0,0]+=.01
    wrong_p={k:v.copy() for k,v in source.items()};wrong_p['p'][0,0]+=1e-5
    early={k:v[source['s']<100-1e-6].copy() for k,v in source.items()}
    q=np.linspace(0,100,1601);values={}
    for side,t in traces.items():
        values[side]=(CubicHermiteSpline(t['s'],t['xyz'],t['tangent'])(q),
                      CubicHermiteSpline(t['s'],t['group'],1/np.linalg.norm(t['p'],axis=1))(q))
    actual=metric(*values['ah'],*values['ph']);adequate={k:2*actual[k] for k in LIMITS}
    insufficient={k:actual[k]/2 for k in LIMITS}
    if not any(actual[k]>0 for k in LIMITS):raise ValueError('Need a nonidentical real trace pair for the budget negative control')
    fixture=PPoly(np.array([[0.,1.,0.],[1.,2.,1e9],[0.,0.,0.],[0.,0.,0.]]),[0.,99.,101.,102.])
    checks={'baseline_native_and_initial_accepted':native_valid(source) and initial_valid(source,xyz,p),
            'wrong_initial_position_rejected':not initial_valid(wrong_position,xyz,p),
            'wrong_initial_momentum_rejected':not initial_valid(wrong_p,xyz,p),
            'early_termination_rejected':not native_valid(early),
            'adequate_cross_budget_accepted':budget_valid(actual,adequate),
            'insufficient_cross_budget_rejected':not budget_valid(actual,insufficient),
            'outside_support_curvature_excluded':curvature_max(fixture)==10. and curvature_max(fixture,False)==2e9,
            'shared_position_bias_rejected':not gate(dict(rms_m=3.,max_m=3.,landing_m=3.,group_rel=0.)),
            'shared_group_bias_rejected':not gate(dict(rms_m=0.,max_m=0.,landing_m=0.,group_rel=3e-6)),
            'nonfinite_rejected':not gate(dict(rms_m=float('nan'),max_m=0.,landing_m=0.,group_rel=0.))}
    return checks,{'source_case':case,'position_mutation_km':.01,'momentum_mutation':1e-5,
                   'early_final_s_km':float(early['s'][-1]),'actual_cross_metric':actual,
                   'adequate_budget':adequate,'insufficient_budget':insufficient}


def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--parent',type=Path,required=True)
    p.add_argument('--analysis',default='reference_r3_v4');p.add_argument('--output',default='reference_independent_v1');a=p.parse_args()
    out=a.run/a.output;out.mkdir(exist_ok=True);source=out/'source.py'
    if source.exists() and sha(source)!=sha(__file__):raise ValueError('New verifier version requires new output')
    if not source.exists():shutil.copy2(__file__,source)
    rows=[];verified_jobs=[]
    for job in read(a.parent/'experiments/manifest.json')['jobs']:
        if job['role']!='R3' or job['level']!=3:continue
        case=job['job_id'].rsplit('_t',1)[0]
        if not (a.run/a.analysis/case/'result.json').exists():continue
        result=verify_case(a.parent,a.run,a.run/a.analysis,job);rows.append(result);verified_jobs.append(job);save(out/(case+'.json'),result)
        print(case,result['passed'],[k for k,v in result['checks'].items() if not v],flush=True)
    if not verified_jobs:raise ValueError('No completed cases to verify')
    negative,negative_inputs=negative_controls(a.parent,a.run,verified_jobs[0])
    save(out/'verification.json',{'case_count':len(rows),'passed_count':sum(r['passed'] for r in rows),
         'all_16_reproduced':len(rows)==16 and all(r['passed'] for r in rows) and all(negative.values()),
         'negative_controls':negative,'negative_control_inputs':negative_inputs,'cases':rows,'source_sha256':sha(__file__),
         'notes':'Independent native/MP-record recomputation; no expensive MP integration rerun. Numerical error estimates are not interval proofs.'})
    if not all(r['passed'] for r in rows) or not all(negative.values()):raise SystemExit(1)


if __name__=='__main__':main()
