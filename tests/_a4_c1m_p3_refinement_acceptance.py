"""Re-evaluate allowed refinements without replacing original failures."""
from pathlib import Path
import argparse,shutil
import numpy as np
import _a4_c1m_p3_compare as c

def evaluate(run):
    d=run/'experiments';extra=d/'refinements_v1';out=run/'analysis/p3_refined_v1'
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True);shutil.copy2(Path(__file__),out/'analysis_source.py')
    old=c.exp.read(run/'analysis/p3_v2/summary.json');manifest=c.exp.read(d/'manifest.json');ejobs=c.exp.read(extra/'manifest.json')['jobs'];rows=[]
    for original in old['cases']:
        if original['role']!='R3':rows.append(dict(original,source_analysis=str((run/'analysis/p3_v2').resolve()),formal_query_step_km=.25));continue
        row=dict(original);jobs=[j for j in manifest['jobs'] if j['job_id'] in row['jobs']];jobs.sort(key=lambda j:j['level']);added=[j for j in ejobs if j['job_id'].startswith(row['case_id']+'_t')];added.sort(key=lambda j:j['level'])
        ahjobs=jobs[-1:]+added if added else jobs
        ah=[c.read_trace(extra if j['level']>3 else d,j,'ah') for j in ahjobs];ph=[c.read_trace(d,j,'ph') for j in jobs]
        steps=(1.,.5,.25) if row['checks']['sampling'] and not added else (.25,.125,.0625)
        q=c.query(steps[-1]);aval=[c.at(t,q) for t in ah];pval=[c.at(t,q) for t in ph];actual=c.metric(aval[-1],pval[-1]);integration={'ah':c.controls(aval),'ph':c.controls(pval)}
        sampled=[c.metric(c.at(ah[-1],c.query(s)),c.at(ph[-1],c.query(s))) for s in steps]
        delta12={k:abs(sampled[1][k]-sampled[0][k]) for k in c.LIMITS};sampling={k:abs(sampled[2][k]-sampled[1][k]) for k in c.LIMITS};monotonic={k:sampling[k]<=delta12[k] for k in c.LIMITS}
        ode=c.p2.load_ode(c.exp.COEFFICIENT);field,proof_reduction=c.equivalent_reference_field(ode.medium(jobs[-1],False));max_s=max(float(c.p2.prefix(t)['s'][-1]) for t in ah+ph)
        ref=c.p2.RadialReference(ph[-1],field,ode,c.exp.COEFFICIENT,jobs[-1]['frequency_MHz'],5,max_s);proof=c.p2.field_bounds(field,ref.alpha,ref.E,max_s)
        reconstruction={side:c.p2.reconstruction_bound(t,'derivative_hermite',ref,proof,side) for side,t in [('ah',ah[-1]),('ph',ph[-1])]};bounds={k:sum(r['bounds'][k] for r in reconstruction.values()) for k in c.LIMITS}
        references=[];reference_warnings=[]
        for level in (3,4):
            lower=c.p2.RadialReference(ph[-1],field,ode,c.exp.COEFFICIENT,jobs[-1]['frequency_MHz'],level,max_s);references.append(lower.query(q));reference_warnings.extend(lower.warnings)
            c.exp.write(out/(row['case_id']+f'_quadrature_l{level}.json'),{'estimates':list(lower.uncertainty_records.values()),'warnings':lower.warnings})
        references.append(ref.query(q));reference_controls=c.controls(references)
        reference_valid=reference_controls['passed'] and not reference_warnings and not ref.warnings
        budget={k:3*(integration['ah']['delta23'][k]+integration['ph']['delta23'][k]+sampling[k]+bounds[k]) for k in c.LIMITS}
        initial_pos=max(float(np.linalg.norm(a['xyz'][0]-b['xyz'][0])*1000) for a,b in zip(ah,ph));initial_p=max(float(np.linalg.norm(a['p'][0]-b['p'][0])) for a,b in zip(ah,ph))
        checks={'thresholds':all(actual[k]<=c.LIMITS[k] for k in c.LIMITS),'budget':all(actual[k]<=budget[k] for k in c.LIMITS),'integration':all(v['passed'] for v in integration.values()),'sampling':all(monotonic[k] and sampling[k]<=.2*c.LIMITS[k] for k in c.LIMITS),'reconstruction':reference_valid and all(bounds[k]<=.2*c.LIMITS[k] for k in c.LIMITS),'initial_state':initial_pos<=1e-3 and initial_p<=1e-10,'scope':all(abs(t['s'][-1]-100)<1e-8 and t['termination']=='fixed_geometric_path_completed' for t in ah) and all(t['s'][-1]>=100 for t in ph)}
        levels={side:[{'level':j['level'],'job_id':j['job_id'],'raw_path':t['raw_path'],'raw_sha256':t['raw_sha256'],'termination':t['termination'],'native_s_end_km':float(t['s'][-1]),'native_group_end_km':float(t['group'][-1])} for j,t in zip(js,ts)] for side,js,ts in [('ah',ahjobs,ah),('ph',jobs,ph)]}
        row.update(actual=actual,integration=integration,levels=levels,initial_position_m=initial_pos,initial_p_difference=initial_p,resampling={'steps_km':steps,'metrics':sampled,'delta12':delta12,'monotonic':monotonic,'final_bound':sampling},reconstruction={'validated':reference_valid,'sides':reconstruction,'bound':bounds},checks=checks,budget=budget,strict_passed=all(checks.values()),original_strict_passed=original['strict_passed'],original_record=str((run/'analysis/p3_v2'/(row['case_id']+'.json')).resolve()),AH_levels_used=[j['level'] for j in ahjobs],PH_levels_used=[1,2,3],additional_sampling_levels_used=0 if steps==c.STEPS else 2)
        row.update(formal_query_step_km=steps[-1],source_analysis=str(out.resolve()),independent_reference=dict(original['independent_reference'],convergence=reference_controls,AH=c.metric(aval[-1],references[-1]),PH=c.metric(pval[-1],references[-1])))
        c.exp.write(out/(row['case_id']+'_quadrature.json'),{'estimates':list(ref.uncertainty_records.values()),'warnings':ref.warnings,'field_reduction':proof_reduction})
        c.exp.write(out/(row['case_id']+'.json'),row);rows.append(row);print(row['case_id'],row['strict_passed'],flush=True)
    summary=dict(old,cases=rows,formal_passed_count=sum(r['strict_passed'] for r in rows),p3_scope_passed=all(r['strict_passed'] for r in rows),recommend_p4=False,analysis_sha256=c.common.sha(Path(__file__)),parent_summary=str((run/'analysis/p3_v2/summary.json').resolve()))
    c.exp.write(out/'summary.json',summary)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();evaluate(a.run.resolve())
