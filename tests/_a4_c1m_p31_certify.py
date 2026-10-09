"""P3.1 R3 independent-reference certification; never edits P3 artifacts."""
from __future__ import annotations
import argparse,json,os,shutil,time
from pathlib import Path
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import numpy as np
import mpmath as mp
import _a4_c1m_p3_compare as comp
from _a4_c1m_p31_reference import AffineReference,exact,self_check

LIMITS=comp.LIMITS
def write(path,data):comp.common.write_json(path,data)
def metric(xyz,group,reference):
    distance=np.linalg.norm(xyz-reference['xyz'],axis=1)*1000
    return {'rms_m':float(np.sqrt(np.mean(distance**2))),'max_m':float(max(distance)),'landing_m':float(distance[-1]),'group_rel':float(abs(group[-1]-reference['group'][-1])/abs(reference['group'][-1]))}
def gate(errors):return all(errors[k]<=.2*v for k,v in LIMITS.items())
def decimal_difference(a,b,dps=170):
    with mp.workdps(dps):
        pos=mp.mpf(0);group=mp.mpf(0)
        for x,y in zip(a,b):
            pos=max(pos,mp.sqrt(sum((mp.mpf(x[i])-mp.mpf(y[i]))**2 for i in range(3)))*1000)
            group=max(group,abs(mp.mpf(x[3])-mp.mpf(y[3])))
        return {'position_max_m':mp.nstr(pos,80),'group_max_km':mp.nstr(group,80)}
def certify_case(parent,out,job):
    case=job['job_id'].rsplit('_t',1)[0]; target=out/case;target.mkdir()
    oldpath=parent/'analysis/p3_refined_v1'/(case+'.json');old=comp.exp.read(oldpath)
    d=parent/'experiments';ph=comp.read_trace(d,job,'ph')
    last=old['levels']['ah'][-1];rawpath=Path(last['raw_path'])
    if comp.common.sha(rawpath)!=last['raw_sha256']:raise ValueError('AH final hash mismatch')
    with np.load(rawpath) as f:ah={k:np.array(f[k]) for k in ('s','xyz','p','group','tangent')}
    traces={'ah':ah,'ph':ph};q=comp.query(.0625)
    native={side:np.flatnonzero(t['s']<=100.) for side,t in traces.items()}
    queries=np.unique(np.r_[q,*[t['s'][native[side]] for side,t in traces.items()]])
    ode=comp.p2.load_ode(comp.exp.COEFFICIENT);field=ode.medium(job,False)
    np.savez_compressed(target/'frozen_field.npz',h=field['h'],coefficients=field['coefficients'])
    source_records=[{'path':str(p.resolve()),'sha256':comp.common.sha(p)} for p in (oldpath,rawpath,Path(ph['raw_path']),Path(job['grid_mat']))]
    write(target/'input_identity.json',{'job':job,'sources':source_records,'query_count':len(queries),'scope':'first_100_km','initial_reference':'PH frozen native state, not fitted'})
    results={};timings={};proof=None;reuse_records=[]
    for dps,order in ((80,32),(120,32),(160,16),(160,24),(160,32)):
        start=time.time();reference=AffineReference(field,ph['xyz'][0],ph['p'][0],comp.exp.COEFFICIENT,job['frequency_MHz'],dps,order)
        prior=out.parent/'reference_r3_v1'/case
        old_npz=prior/f'reference_dps{dps}_n{order}.npz';old_decimal=prior/f'reference_dps{dps}_n{order}_decimal.json'
        source_ok=(prior.parent/'source/_a4_c1m_p31_reference.py').exists() and comp.common.sha(prior.parent/'source/_a4_c1m_p31_reference.py')==comp.common.sha(Path(__file__).with_name('_a4_c1m_p31_reference.py'))
        identity_ok=(prior/'input_identity.json').exists() and comp.exp.read(prior/'input_identity.json')==comp.exp.read(target/'input_identity.json')
        if old_npz.exists() and old_decimal.exists() and source_ok and identity_ok:
            with np.load(old_npz) as raw:
                if not np.array_equal(raw['s'],queries):raise ValueError('Reuse query mismatch')
                result={k:np.array(raw[k]) for k in ('xyz','p','group')}
            result['decimal']=comp.exp.read(old_decimal)
            reuse_records += [{'path':str(x.resolve()),'sha256':comp.common.sha(x)} for x in (old_npz,old_decimal)]
        else:result=reference.query(queries)
        results[(dps,order)]=result;proof=reference.proof
        np.savez_compressed(target/f'reference_dps{dps}_n{order}.npz',s=queries,xyz=result['xyz'],p=result['p'],group=result['group'])
        write(target/f'reference_dps{dps}_n{order}_decimal.json',result['decimal']);timings[f'{dps}_{order}']=time.time()-start
    final=results[(160,32)];precision=[decimal_difference(results[(a,32)]['decimal'],results[(b,32)]['decimal']) for a,b in ((80,120),(120,160))]
    quadrature=[decimal_difference(results[(160,a)]['decimal'],results[(160,b)]['decimal']) for a,b in ((16,24),(24,32))]
    # Conversion is measured under sufficient precision, never folded into sampling.
    with mp.workdps(170):converted=[[mp.nstr(exact(x),170) for x in (*pos,g)] for pos,g in zip(final['xyz'],final['group'])]
    rounding=decimal_difference(final['decimal'],converted)
    ref_est=max(float(precision[-1]['position_max_m']),float(quadrature[-1]['position_max_m']))+float(rounding['position_max_m'])
    group_est=max(float(precision[-1]['group_max_km']),float(quadrature[-1]['group_max_km']))+float(rounding['group_max_km'])
    iq=np.searchsorted(queries,q);formal={k:final[k][iq] for k in ('xyz','p','group')}
    errors={};native_errors={};reconstructed={}
    for side,t in traces.items():
        ids=native[side];ix=np.searchsorted(queries,t['s'][ids]);r={k:final[k][ix] for k in ('xyz','p','group')}
        native_errors[side]=metric(t['xyz'][ids],t['group'][ids],r)
        native_errors[side]['group_rel']=float(np.max(np.abs(t['group'][ids]-r['group']))/abs(formal['group'][-1]))
        v=comp.at(t,q);reconstructed[side]=v;errors[side]=metric(v['xyz'],np.array([v['g100']]),formal)
    actual=comp.metric(reconstructed['ah'],reconstructed['ph'])
    reconstruction={side:old['reconstruction']['sides'][side]['bounds'] for side in traces}
    sampling=old['resampling']['final_bound']
    refbudget={'rms_m':ref_est,'max_m':ref_est,'landing_m':ref_est,'group_rel':group_est/float(formal['group'][-1])}
    envelopes={side:{k:(native_errors[side]['max_m'] if k in ('rms_m','max_m','landing_m') else native_errors[side][k])+refbudget[k] for k in LIMITS} for side in traces}
    budget={k:3*(sum(envelopes[s][k]+reconstruction[s][k] for s in traces)+sampling[k]) for k in LIMITS}
    reference_stability_limits={'position_max_m':1e-12,'group_max_km':1e-15}
    stable=all(float(precision[-1][k])<=float(precision[0][k]) and all(float(delta[k])<=reference_stability_limits[k] for delta in quadrature) for k in precision[0])
    checks={'initial_state':bool(old['checks']['initial_state']),'scope':bool(old['checks']['scope'] and old['topology']['coverage_complete']),'reference_stable':stable,'reference_error':gate(refbudget),'side_errors':all(gate(v) for v in errors.values()),'native_errors':all(gate(v) for v in native_errors.values()),'reconstruction':bool(old['reconstruction']['validated']) and all(gate(v) for v in reconstruction.values()),'sampling':gate(sampling),'cross_thresholds':all(actual[k]<=v for k,v in LIMITS.items()),'budget':all(actual[k]<=budget[k] for k in LIMITS)}
    result={'case_id':case,'strict_convergence_passed':bool(old['strict_passed']),'reference_certified_passed':bool(all(checks.values())),'checks':checks,'actual':actual,'native_errors':native_errors,'formal_errors':errors,'side_envelopes':envelopes,'reconstruction':reconstruction,'sampling':sampling,'reference_budget':refbudget,'budget':budget,'precision_changes':precision,'quadrature_changes':quadrature,'float_conversion':rounding,'proof':proof,'sources':source_records,'timing_s':timings,'interval_certified':False,'uncertainty_kind':'independent numerical estimate, not rigorous interval proof','reconstruction_source':str(oldpath.resolve())}
    result['reused_quadrature']=reuse_records
    result['reference_stability_limits']=reference_stability_limits
    comp.p2.write_csv(target/'formal_points.csv',[{'s_km':float(s),'ah_reference_m':float(np.linalg.norm(reconstructed['ah']['xyz'][i]-formal['xyz'][i])*1000),'ph_reference_m':float(np.linalg.norm(reconstructed['ph']['xyz'][i]-formal['xyz'][i])*1000),'cross_m':float(np.linalg.norm(reconstructed['ah']['xyz'][i]-reconstructed['ph']['xyz'][i])*1000)} for i,s in enumerate(q)])
    write(target/'result.json',result);return result
def main():
    p=argparse.ArgumentParser();p.add_argument('--parent',type=Path,required=True);p.add_argument('--run',type=Path,required=True);p.add_argument('--case');args=p.parse_args()
    out=args.run/'reference_r3_v2';out.mkdir(exist_ok=True)
    source=out/'source';source.mkdir(exist_ok=True)
    for name in ('_a4_c1m_p31_reference.py','_a4_c1m_p31_certify.py','_a4_c1m_p3_compare.py','_a4_c1m_p2_analysis_core_v2.py','_a4_c1m_p3_experiments.py','_run_a4_c1m_attribution.py'):
        original=Path(__file__).with_name(name);dest=source/name
        if dest.exists() and comp.common.sha(dest)!=comp.common.sha(original):raise ValueError('Source changed; new run required')
        if not dest.exists():shutil.copy2(original,dest)
    write(out/'source_manifest.json',[{'path':str(x.resolve()),'sha256':comp.common.sha(x)} for x in source.iterdir()])
    tests=self_check();tests['shared_bias_rejected']=not gate({'rms_m':3.,'max_m':3.,'landing_m':3.,'group_rel':0.});assert all(tests.values());write(out/'self_check.json',tests)
    jobs=[j for j in comp.exp.read(args.parent/'experiments/manifest.json')['jobs'] if j['role']=='R3' and j['level']==3]
    rows=[]
    for j in jobs:
        case=j['job_id'].rsplit('_t',1)[0]
        if args.case and args.case!=case:continue
        path=out/case/'result.json'
        if path.exists():rows.append(comp.exp.read(path));continue
        print('START '+case,flush=True);r=certify_case(args.parent,out,j);rows.append(r);print('DONE '+case+' '+str(r['reference_certified_passed']),flush=True)
    write(out/'summary.json',{'case_count':len(rows),'reference_certified_count':sum(x['reference_certified_passed'] for x in rows),'all_16_certified':len(rows)==16 and all(x['reference_certified_passed'] for x in rows),'cases':rows,'self_checks':tests})
if __name__=='__main__':main()
