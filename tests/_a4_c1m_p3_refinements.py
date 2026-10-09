"""At most two pre-existing legal AH levels; retain original P3 judgments."""
from pathlib import Path
import argparse,json
from concurrent.futures import ProcessPoolExecutor,as_completed
import _a4_c1m_p3_compare as c

def refine(run):
    d=run/'experiments';out=d/'refinements_v1'
    if out.exists():raise FileExistsError(out)
    out.mkdir();manifest=c.exp.read(d/'manifest.json');summary=c.exp.read(run/'analysis/p3_v2/summary.json');sources=c.exp.read(d/'ah_source_sha256.json')
    if any(c.common.sha(Path(p))!=h for p,h in sources.items()):raise ValueError('Production execution sources changed')
    jobs=[];selected=[]
    for row in summary['cases']:
        if row['role']!='R3' or row.get('integration',{}).get('ah',{}).get('passed',True):continue
        selected.append(row['case_id']);base=next(j for j in manifest['jobs'] if j['job_id']==row['jobs'][-1])
        for level in (4,5):
            j=dict(base,job_id=row['case_id']+f'_t{level}',level=level,role='R3_extra',parent_job_sha256=base['job_sha256'])
            j.pop('job_sha256');j['job_sha256']=c.exp.digest(j);jobs.append(j)
    c.exp.write(out/'manifest.json',{'jobs':jobs,'reason':'Original AH self-refinement nonmonotonic, no threshold waiver','source_sha256':c.common.sha(Path(__file__)),'selected':selected,'PH_extra_levels':0,'maximum_additional_AH_levels':2})
    with ProcessPoolExecutor(max_workers=2) as pool:
        fs=[pool.submit(c.exp.ah_one,j,str(out/'ah'),sources) for j in jobs]
        for f in as_completed(fs):r=f.result();print(r['job_id'],r['status'],flush=True)
    rows=[]
    for row in summary['cases']:
        if row['role']!='R3':continue
        base=next(j for j in manifest['jobs'] if j['job_id']==row['jobs'][-1]);ah=c.read_trace(d,base,'ah');ph=c.read_trace(d,base,'ph')
        result={'case_id':row['case_id'],'original_strict_passed':row['strict_passed'],'PH_original_controls':row['integration']['ph']}
        extra=[j for j in jobs if j['job_id'].startswith(row['case_id']+'_t')]
        if extra:
            values=[c.at(ah,c.query(.25))]+[c.at(c.read_trace(out,j,'ah'),c.query(.25)) for j in extra]
            result['AH_final_three_levels']=[3,4,5];result['AH_final_controls']=c.controls(values)
        steps=(1.,.5,.25,.125,.0625);samples=[c.metric(c.at(ah,c.query(h)),c.at(ph,c.query(h))) for h in steps]
        result['sampling_steps_km']=steps;result['sampling_metrics']=samples;result['sampling_last_deltas']=[{k:abs(a[k]-b[k]) for k in c.LIMITS} for a,b in zip(samples[:-1],samples[1:])]
        result['original_record_replaced']=False;result['strict_passed']=False;rows.append(result)
    c.exp.write(out/'results.json',{'rows':rows,'raw_records_replaced':False,'uncertainty_waiver':False})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();refine(a.run.resolve())
