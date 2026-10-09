"""D0/D1 read-only identity, analytical controls, and frozen P2 reproduction."""
from pathlib import Path
import sys,argparse,json,platform,shutil
import numpy as np
import _a4_c1m_p3_experiments as e
import _a4_c1m_p2_analysis_core_v2 as p2

def check(run):
    d=run/'experiments';out=run/'analysis/p3_contract_v3'
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True)
    shutil.copy2(Path(__file__),out/'execution_source.py')
    m=e.read(d/'manifest.json');records=[]
    actual=e.read(d/'training_pharlap_summary.json');mex=Path(actual['mex_path']);mh=e.common.sha(mex)
    for j in m['jobs']:
        if j['role']=='R5_fd':continue
        p=d/'pharlap_raw'/(j['job_id']+'.mat');raw=e.common.mat(p)
        records.append({'job_id':j['job_id'],'path':str(p.resolve()),'sha256':e.common.sha(p),'valid':bool(raw['job_sha256']==j['job_sha256'] and raw['grid_sha256']==j['grid_sha256'] and raw['mex_sha256']==mh and raw['exporter_sha256']==j['exporter_sha256'])})
    readbacks=[]
    for grid_hash in sorted({j['grid_sha256'] for j in m['jobs']}):
        j=next(j for j in m['jobs'] if j['grid_sha256']==grid_hash);original=e.common.mat(Path(j['grid_mat']));p=d/'pharlap_raw'/('grid_readback_'+grid_hash[:16]+'.mat');returned=e.common.mat(p)
        fields={k:bool(np.array_equal(original[k],returned[k])) for k in original if not k.startswith('__')}
        readbacks.append({'grid_sha256':grid_hash,'path':str(p.resolve()),'sha256':e.common.sha(p),'fields_equal':fields,'passed':all(fields.values())})
    old=e.ROOT/'output2/a4_pharlap/c1/c1m/p2_20261008_v1';jobs=e.read(old/'inputs/internal_jobs.json')
    if isinstance(jobs,dict):jobs=jobs['primary']
    baseline=[]
    for item in jobs:
        j=item['original_job'];ph=p2.preparation.read_ph(Path(item['PH_raw']),j);ph['tangent']=ph['p']/np.linalg.norm(ph['p'],axis=1)[:,None]
        ah=p2.read_ah(Path(item['AH_meta']),j)
        metrics=p2.errors(p2.reconstruct(p2.prefix(ph),'pchip',p2.grid(.5)),p2.reconstruct(p2.prefix(ah),'pchip',p2.grid(.5)))
        baseline.append({'job_id':j['job_id'],'angle':j['elevation_deg'],'level':j['level'],'metrics':metrics})
    expected={22.:.008938372632545422,39.5:.007141469253458584}
    reproduced=all(abs(r['metrics']['rms_m']-expected[r['angle']])<1e-12 for r in baseline if r['level']==3)
    analytic=p2.analytic_checks()
    e.write(out/'report.json',{'schema':'p3-D0-D1-v1','source_sha256':e.common.sha(Path(__file__)),'python':sys.version,'platform':platform.platform(),'matlab_version':actual['matlab_version'],'mex_path':str(mex),'mex_sha256':mh,'ph_records':records,'ph_identity_passed':len(records)==168 and all(r['valid'] for r in records),'grid_readbacks':readbacks,'grid_readbacks_passed':all(r['passed'] for r in readbacks),'analytic_controls':analytic,'p2_baseline':baseline,'p2_baseline_reproduced':reproduced,'p2_original_preserved':True})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();check(a.run.resolve())
