"""Independent re-evaluation from frozen native states and precision records.

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
from scipy.interpolate import CubicHermiteSpline, PPoly, PchipInterpolator
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


PHYSICAL_LIMITS={'rms_m':2.,'max_m':10.,'landing_m':10.,'group_rel':2e-6}


def physics_gate(values, fraction=1.):
    return all(np.isfinite(values[k]) and 0 <= values[k] <= fraction*v
               for k,v in PHYSICAL_LIMITS.items())


def physics_splines(trace):
    return (CubicHermiteSpline(trace['s'],trace['xyz'],trace['tangent']),
            CubicHermiteSpline(trace['s'],trace['group'],1/np.linalg.norm(trace['p'],axis=1)))


def physics_cross(a,b):
    end=min(float(a['s'][-1]),float(b['s'][-1]))
    q=np.unique(np.r_[np.arange(0.,end,.5),end])
    ax,_=physics_splines(a);bx,_=physics_splines(b)
    error=np.linalg.norm(ax(q)-bx(q),axis=1)*1000
    return dict(rms_m=float(np.sqrt(np.mean(error**2))),max_m=float(error.max()),
                landing_m=float(np.linalg.norm(a['xyz'][-1]-b['xyz'][-1])*1000),
                group_rel=float(abs(a['group'][-1]-b['group'][-1])/b['group'][-1]))


def physics_native_checks(trace,case,field,coefficient,residual_limit):
    s,xyz,p,group=(trace[k] for k in ('s','xyz','p','group'))
    pp=PPoly(field['coefficients'],field['h'],extrapolate=False)
    height=np.linalg.norm(xyz,axis=1)-6371.
    ne=pp(height)
    for edge,mask in ((0,height < field['h'][0]),(-1,height > field['h'][-1])):
        if np.any(mask) and float(pp(field['h'][edge])) == 0.:
            ne[mask]=0.
    residual=np.sum(p*p,axis=1)-(1-coefficient*ne/(case['frequency_MHz']*1e6)**2)
    lat,lon=np.deg2rad([case.get('tx_lat_deg',30.),case.get('tx_lon_deg',120.)])
    up=np.array([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)])
    north=np.array([-np.sin(lat)*np.cos(lon),-np.sin(lat)*np.sin(lon),np.cos(lat)])
    east=np.array([-np.sin(lon),np.cos(lon),0.])
    elevation,bearing=np.deg2rad([case['elevation_deg'],case.get('bearing_deg',0.)])
    origin=case['origin_height_km']
    n2=1-coefficient*float(pp(origin))/(case['frequency_MHz']*1e6)**2
    expected_xyz=(6371.+origin)*up
    expected_p=np.sqrt(n2)*(np.sin(elevation)*up+np.cos(elevation)*(np.cos(bearing)*north+np.sin(bearing)*east))
    angular=np.cross(xyz,p)
    norm=float(np.linalg.norm(angular[0]))
    roots=pp.derivative().roots(extrapolate=False)
    probes=np.r_[field['h'],roots[np.isfinite(roots)]]
    local=case['scope']=='first_100_km'
    termination=str(np.asarray(trace['termination']).item())
    radial=np.sum(xyz*p,axis=1)
    checks=dict(finite=bool(all(np.isfinite(trace[k]).all() for k in ('s','xyz','p','group','tangent'))),
                native_shapes=bool(xyz.shape==p.shape==trace['tangent'].shape==(len(s),3) and group.shape==s.shape),
                monotonic_arc_group=bool(s[0]==0 and group[0]==0 and np.all(np.diff(s)>0) and np.all(np.diff(group)>0)),
                initial_position=bool(np.linalg.norm(xyz[0]-expected_xyz)*1000<=1e-5),
                initial_momentum=bool(np.linalg.norm(p[0]-expected_p)<=1e-10),
                initial_dispersion=bool(abs(residual[0])<=1e-10),
                dispersion_residual=bool(np.isfinite(residual).all() and np.max(abs(residual))<=residual_limit),
                tangent=bool(np.max(abs(trace['tangent']-p/np.linalg.norm(p,axis=1)[:,None]))<=1e-12),
                coverage_termination=bool(abs(s[-1]-100.)<1e-8 if local else termination=='ground'),
                topology=bool(np.all(radial>0) if local else
                              (radial[0]>0 and radial[-1]<0 and abs(height[-1])<=1e-6
                               and np.count_nonzero((radial[:-1]>0)&(radial[1:]<=0))==1)),
                nonnegative_segments=bool(np.isfinite(pp(probes)).all() and np.min(pp(probes))>=0.),
                smooth_ground_entry=bool(local or (pp(0.)==0 and pp.derivative()(0.)==0)),
                radial_invariant=bool(norm>0 and np.max(np.linalg.norm(angular-angular[0],axis=1)/norm)<=2e-6),
                plane=bool(norm>0 and np.max(abs(xyz@(angular[0]/norm)))*1000<=.2))
    return checks


def physics_precision_error(record,reference,control_error,control_group):
    rows=record['precision_quadrature_rows']
    if [(r['dps'],r['order']) for r in rows] != [(80,32),(120,32),(160,32),(160,48),(160,64)]:
        raise ValueError('Missing registered precision or quadrature levels')
    with mp.workdps(170):
        changes=[]
        for a,b in zip(rows[:-1],rows[1:]):
            if len(a['cumulative'])!=len(b['cumulative']) or not a['cumulative']:
                raise ValueError('Precision records use inconsistent knot partitions')
            changes.append([float(max(abs(mp.mpf(x[k])-mp.mpf(y[k]))
                                  for x,y in zip(a['cumulative'],b['cumulative']))) for k in range(3)])
        final=np.array([float(mp.mpf(v)) for v in rows[-1]['integrals']])
    from _a4_c1m_p31_fullpath import LocalRadialPath
    scipy_value=(reference.integrals(reference.r0,reference.r_end)
                 if isinstance(reference,LocalRadialPath) else 2*reference.half)
    difference=abs(final-scipy_value)
    estimated=2*np.max(changes,axis=0)+difference
    position=float((estimated[0]+reference.r0*estimated[1])*1000)+control_error
    group=float(estimated[2]/reference.group_end)+control_group
    return {k:position if k!='group_rel' else group for k in PHYSICAL_LIMITS},changes,difference.tolist()


def physics_case_verify(run,case,manifest,contract):
    from _a4_c1m_p31_fullpath import LocalRadialPath,RadialFullPath
    folder=run/'cases'/case['case_id']
    paths=[folder/name for name in ('ah_l1.npz','ah_l2.npz','ah_l3.npz','ah_fd0.25.npz',
                                   'ah_fd0.0625.npz','reference.npz','cartesian_control.npz',
                                   'reference_precision.json','result.json')]
    paths.append(run/case['medium_path'])
    hashes=[dict(path=str(path.relative_to(run)),sha256=sha(path)) for path in paths]
    def load(path):
        with np.load(path,allow_pickle=False) as data:return {k:data[k].copy() for k in data.files}
    traces=[load(path) for path in paths[:3]]
    fd=[load(paths[3]),traces[-1],load(paths[4])]
    truth,control=load(paths[5]),load(paths[6])
    precision,reported=read(paths[7]),read(paths[8])
    field=load(paths[-1]);field={'h':field['h'],'coefficients':field['coefficients'],'Ne':field['Ne']}
    if (field['Ne'].shape!=field['h'].shape or
            not np.array_equal(PchipInterpolator(field['h'],field['Ne'],extrapolate=False).c,field['coefficients'])):
        raise ValueError('Actual frozen polynomial does not reproduce registered radial PCHIP nodes')
    final=traces[-1]
    constructor=LocalRadialPath if case['scope']=='first_100_km' else RadialFullPath
    reference=constructor(field,final['xyz'][0],final['p'][0],manifest['coefficient'],case['frequency_MHz'])
    checks={}
    registration=read(run/'reference_registration.json')
    registered=next(row for row in registration['cases'] if row['case_id']==case['case_id'])
    checks['reference_preregistered']=bool(registered['reference_registered'] and registration['AH_execution_started'] is False)
    checks['registered_initial_position']=bool(np.array_equal(final['xyz'][0],registered['initial_xyz']))
    checks['registered_initial_momentum']=bool(np.array_equal(final['p'][0],registered['initial_p']))
    checks['registered_topology']=registered['topology']==('local_upward' if case['scope']=='first_100_km' else 'single_turn_ground')
    checks['registered_natural_path']=bool(np.isclose(reference.s_end,registered['natural_s_end'],rtol=1e-12,atol=1e-9))
    checks['registered_natural_group']=bool(np.isclose(reference.group_end,registered['natural_group_end'],rtol=1e-12,atol=1e-9))
    checks['registered_natural_endpoint']=bool(np.max(abs(reference.endpoint-registered['natural_endpoint']))<=1e-8)
    hashes.append(dict(path='reference_registration.json',sha256=sha(run/'reference_registration.json')))
    for name,trace in zip(('l1','l2','l3','fd_coarse','fd_fine'),traces+fd[::2]):
        for key,value in physics_native_checks(trace,case,field,manifest['coefficient'],contract['dispersion_residual_limit']).items():
            checks[name+'_'+key]=value
    q=truth['s'];recomputed=reference.query(q)
    checks['reference_raw_positions']=bool(np.max(np.linalg.norm(truth['xyz']-recomputed['xyz'],axis=1))*1000<=1e-4)
    checks['reference_raw_group']=bool(np.max(abs(truth['group']-recomputed['group']))/reference.group_end<=1e-10)
    checks['reference_raw_momenta']=bool(np.max(abs(truth['p']-recomputed['p']))<=1e-10)
    checks['reference_raw_tangents']=bool(np.max(abs(truth['tangent']-recomputed['tangent']))<=1e-10)
    checks['natural_coverage']=bool(abs(final['s'][-1]-reference.s_end)*1000<=10.)
    end=min(float(final['s'][-1]),reference.s_end)
    native=final['s'][final['s']<=end];mids=(native[:-1]+native[1:])/2
    def at(values):
        index=np.searchsorted(q,values)
        if np.any(index>=len(q)) or not np.allclose(q[index],values,rtol=0,atol=1e-10):
            raise ValueError('Reference raw data lacks required native, midpoint or sample nodes')
        return index
    ix,mx=at(native),at(mids)
    x,g=physics_splines(final)
    refx=CubicHermiteSpline(native,recomputed['xyz'][ix],recomputed['tangent'][ix])
    refg=CubicHermiteSpline(native,recomputed['group'][ix],1/np.linalg.norm(recomputed['p'][ix],axis=1))
    native_error=float(np.max(np.linalg.norm(final['xyz'][:len(native)]-recomputed['xyz'][ix],axis=1))*1000)
    reconstruction_position=float(np.max(np.linalg.norm(refx(mids)-recomputed['xyz'][mx],axis=1))*1000)
    reconstruction_group=float(np.max(abs(refg(mids)-recomputed['group'][mx]))/reference.group_end)
    reconstruction={k:reconstruction_position for k in ('rms_m','max_m')}
    reconstruction.update(landing_m=0.,group_rel=reconstruction_group)
    samples=[]
    for spacing in contract['sampling_steps_km']:
        sample_q=np.unique(np.r_[np.arange(0.,end,spacing),end]);index=at(sample_q)
        error=np.linalg.norm(x(sample_q)-recomputed['xyz'][index],axis=1)*1000
        samples.append(dict(step_km=spacing,rms_m=float(np.sqrt(np.mean(error**2))),max_m=float(error.max())))
    sampling={k:max(abs(samples[i][k]-samples[i-1][k]) for i in (1,2)) for k in ('rms_m','max_m')}
    sampling.update(landing_m=0.,group_rel=0.)
    metrics={k:samples[-1][k] for k in ('rms_m','max_m')}
    metrics.update(landing_m=float(np.linalg.norm(final['xyz'][-1]-reference.endpoint)*1000),
                   group_rel=float(abs(final['group'][-1]-reference.group_end)/reference.group_end))
    propagation=dict(metrics,rms_m=native_error,max_m=native_error)
    cq=control['s'];ctruth=reference.query(cq);cy=control['state']
    control_error=float(np.max(np.linalg.norm(cy[:3].T-ctruth['xyz'],axis=1))*1000)
    control_group=float(np.max(abs(cy[6]-ctruth['group']))/reference.group_end)
    control_endpoint_error=float(np.linalg.norm(control['natural_endpoint']-reference.endpoint)*1000)
    control_total_group=float(abs(float(control['natural_group_end'])-reference.group_end)/reference.group_end)
    expected_termination='fixed_geometric_path_completed' if case['scope']=='first_100_km' else 'ground'
    checks['cartesian_natural_termination']=str(control['termination'].item())==expected_termination
    checks['cartesian_natural_coverage']=bool(abs(float(control['natural_s_end'])-reference.s_end)*1000<=10.)
    checks['cartesian_complete_support']=bool(abs(float(cq[-1])-min(float(control['natural_s_end']),reference.s_end))<=1e-8)
    checks['cartesian_raw_reference_positions']=bool(np.max(abs(control['reference_xyz']-ctruth['xyz']))<=1e-7)
    checks['cartesian_raw_reference_group']=bool(np.max(abs(control['reference_group']-ctruth['group']))<=1e-7)
    reference_error,precision_changes,precision_difference=physics_precision_error(
        precision,reference,max(control_error,control_endpoint_error),max(control_group,control_total_group))
    total={k:sum(part[k] for part in (propagation,reconstruction,sampling,reference_error)) for k in PHYSICAL_LIMITS}
    increments=[physics_cross(traces[i-1],traces[i]) for i in (1,2)]
    fd_increments=[physics_cross(fd[i-1],fd[i]) for i in (1,2)]
    monotonic=all(increments[1][k]<=increments[0][k] for k in PHYSICAL_LIMITS)
    strict=bool(monotonic and physics_gate(increments[-1]) and physics_gate(metrics)
                and all(metrics[k]<=3*(increments[-1][k]+reconstruction[k]+sampling[k]+reference_error[k]) for k in PHYSICAL_LIMITS))
    budget=bool(physics_gate(total) and all(physics_gate(part,.1) for part in (reconstruction,sampling,reference_error)))
    fd_passed=bool(physics_gate(fd_increments[-1],.2) and all(fd_increments[1][k]<=fd_increments[0][k] for k in PHYSICAL_LIMITS))
    computed=dict(metrics=metrics,propagation_envelope=propagation,reconstruction=reconstruction,
                  sampling=sampling,reference_error=reference_error,total_error_estimate=total,
                  integration_increments=increments,fd_increments=fd_increments,
                  reference_natural_s_end_km=reference.s_end,reference_natural_group_km=reference.group_end,
                  reference_endpoint=reference.endpoint.tolist())
    differences={}
    for key,value in computed.items():
        if key not in reported:
            differences[key]='missing';continue
        def same(a,b):
            if isinstance(a,dict):return set(a)<=set(b) and all(same(v,b[k]) for k,v in a.items())
            if isinstance(a,list):return len(a)==len(b) and all(same(x,y) for x,y in zip(a,b))
            return bool(np.isclose(a,b,rtol=1e-7,atol=1e-12))
        if not same(value,reported[key]):differences[key]=dict(recomputed=value,reported=reported[key])
    state=all(v for key,v in checks.items() if key.startswith(('l1_','l2_','l3_','fd_')))
    checks['result_numerical_reproduced']=not differences
    checks['fd_monotonic_and_bounded']=fd_passed
    checks['reference_reconstruction_sampling_budgets']=all(physics_gate(part,.1) for part in (reference_error,reconstruction,sampling))
    checks['total_budget']=physics_gate(total)
    checks['acceptance_channel']=strict or budget
    booleans=dict(medium_initial_state_passed=state,strict_convergence_passed=strict,
                  reference_certification_passed=budget,fd_passed=fd_passed,passed=bool(state and (strict or budget)))
    for key,value in booleans.items():
        if reported.get(key) is not value:differences[key]=dict(recomputed=value,reported=reported.get(key))
    checks['result_boolean_reproduced']=all(reported.get(k) is v for k,v in booleans.items())
    return dict(case_id=case['case_id'],passed=bool(all(checks.values())),checks=checks,
                recomputed=computed,channels=booleans,differences=differences,raw_sha256=hashes,
                precision_quadrature_changes=precision_changes,binary64_integral_difference=precision_difference,
                fd_platform_explained=False,reference_estimates_are_interval_proof=False,
                high_precision_integration_rerun=False)


def physics_negative_controls():
    zero={k:0. for k in PHYSICAL_LIMITS}
    checks={'zero_accepted':physics_gate(zero)}
    for name,key,value in (('natural_endpoint','landing_m',10.01),('natural_group','group_rel',2.01e-6)):
        checks[name+'_rejected']=not physics_gate(dict(zero,**{key:value}))
    s=np.array([0.,100.]);momentum=np.array([np.sin(np.deg2rad(22.)),0.,np.cos(np.deg2rad(22.))])
    baseline=dict(s=s,xyz=np.array([6471.,0.,0.])+s[:,None]*momentum,p=np.tile(momentum,(2,1)),
                  group=s.copy(),tangent=np.tile(momentum,(2,1)),termination='max_path')
    case=dict(scope='first_100_km',frequency_MHz=10.,tx_lat_deg=0.,tx_lon_deg=0.,
              origin_height_km=100.,elevation_deg=22.,bearing_deg=0.)
    field=dict(h=np.array([0.,200.]),coefficients=np.zeros((4,1)))
    checks['vacuum_native_control_accepted']=all(physics_native_checks(baseline,case,field,80.6404,1e-10).values())
    units=dict(baseline,xyz=baseline['xyz']*1000.)
    wrong=dict(baseline,p=baseline['p'].copy());wrong['p'][0,0]+=.01
    checks['wrong_units_rejected']=not all(physics_native_checks(units,case,field,80.6404,1e-10).values())
    checks['wrong_initial_state_rejected']=not physics_native_checks(wrong,case,field,80.6404,1e-10)['initial_momentum']
    checks['nonmonotonic_convergence_rejected']=not all(b<=a for a,b in zip([1.,.1],[.5,.2]))
    checks['nonfinite_rejected']=not physics_gate(dict(zero,max_m=float('nan')))
    checks['reference_budget_over_10_percent_rejected']=not physics_gate(dict(zero,max_m=1.01),.1)
    return checks


def physics_verify(run):
    run=run.resolve();out=run/'physics_independent_verification'
    if out.exists():raise ValueError('Refuse independent verification overwrite; preserve existing evidence')
    out.mkdir(parents=True);shutil.copy2(__file__,out/'source.py')
    manifest=read(run/'manifest.json');registry=read(run/'case_registry.json');contract=read(run/'budget_contract.json')
    identity=[]
    root=Path(__file__).resolve().parents[1]
    for kind,base,records in (('source',root,manifest['source_records']),('input',run,manifest['input_records'])):
        for record in records:
            path=(base/record['path']).resolve()
            actual=sha(path) if path.is_file() and path.is_relative_to(base) else None
            row=dict(kind=kind,path=record['path'],expected_sha256=record['sha256'],actual_sha256=actual,passed=actual==record['sha256'])
            if kind=='source':
                snapshot=run/'source_before'/record['path']
                row['snapshot_sha256']=sha(snapshot) if snapshot.is_file() else None
                row['passed']=row['passed'] and row['snapshot_sha256']==record['sha256']
            identity.append(row)
    for record in manifest.get('source_sha256',[]):
        path=Path(record['path']).resolve()
        if path.is_relative_to(root):continue
        snapshot=(run/record['snapshot_path']).resolve()
        actual=sha(path) if path.is_file() else None
        frozen=sha(snapshot) if snapshot.is_relative_to(run) and snapshot.is_file() else None
        identity.append(dict(kind='external_plan',path=str(path),expected_sha256=record['sha256'],
                             actual_sha256=actual,snapshot_sha256=frozen,
                             passed=actual==frozen==record['sha256']))
    for name,key in (('case_registry.json','case_registry_sha256'),('budget_contract.json','budget_contract_sha256'),('exposure_audit.json','exposure_audit_sha256')):
        actual=sha(run/name) if (run/name).is_file() else None
        identity.append(dict(kind='registry',path=name,expected_sha256=manifest.get(key),actual_sha256=actual,passed=actual==manifest.get(key)))
    import sys
    if str(root) not in sys.path:sys.path.insert(0,str(root))
    from config import K_FP
    configuration=(contract['physical_limits']==PHYSICAL_LIMITS and manifest['coefficient']==K_FP**2
                   and manifest.get('compatibility_mode') is None and manifest.get('interpolation')=='radial_pchip'
                   and manifest.get('earth_radius_km')==6371. and manifest.get('gyro_MHz')==0.
                   and manifest.get('collision_s')==0. and len(manifest['integration_profiles'])>=3
                   and type(manifest['projection']) is bool
                   and all(contract[key]==.1 for key in ('reference_fraction','reconstruction_fraction','sampling_fraction'))
                   and contract['fd_last_change_fraction']==.2)
    rows=[]
    for case in registry['cases']:
        try:row=physics_case_verify(run,case,manifest,contract)
        except Exception as exc:
            folder=run/'cases'/case['case_id'];artifacts=[p for p in folder.glob('*') if p.is_file()]
            row=dict(case_id=case['case_id'],passed=False,status='missing_or_invalid_evidence',
                     reason=type(exc).__name__+': '+str(exc),
                     raw_sha256=[dict(path=str(p.relative_to(run)),sha256=sha(p)) for p in sorted(artifacts)])
        rows.append(row);save(out/(case['case_id']+'.json'),row)
        print(case['case_id'],row['passed'],flush=True)
    negatives=physics_negative_controls()
    passed=bool(configuration and all(r['passed'] for r in identity) and len(rows)==26
                and len({r['case_id'] for r in rows})==26 and all(r['passed'] for r in rows) and all(negatives.values()))
    result=dict(schema='p31-physics-independent-verification-v1',case_count=len(rows),
                passed_count=sum(r['passed'] for r in rows),passed=passed,independent_verification_passed=passed,
                configuration_passed=bool(configuration),identity_records=identity,cases=rows,
                negative_controls=negatives,source_sha256=sha(__file__),
                high_precision_integration_rerun=False,
                notes='Raw AH/reference/Cartesian records, frozen polynomial, budgets and precision-record arithmetic independently recomputed; production solver and execution summary booleans are not used as acceptance evidence.')
    save(out/'verification.json',result)
    if not passed:raise SystemExit(1)


def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path);p.add_argument('--parent',type=Path)
    p.add_argument('--physics-run',type=Path)
    p.add_argument('--analysis',default='reference_r3_v4');p.add_argument('--output',default='reference_independent_v1');a=p.parse_args()
    if a.physics_run:
        physics_verify(a.physics_run);return
    if a.run is None or a.parent is None:p.error('--run and --parent are required for historical verification')
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
