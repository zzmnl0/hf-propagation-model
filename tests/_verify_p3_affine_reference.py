"""Exact affine-knot reduction cross-check, preserving both numerical estimates."""
from pathlib import Path
import argparse,shutil
import numpy as np
import _a4_c1m_p3_compare as c

def verify(run):
    d=run/'experiments';out=run/'analysis/p3_affine_reference_equivalence_v2'
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True);shutil.copy2(Path(__file__),out/'verification_source.py')
    j=next(j for j in c.exp.read(d/'manifest.json')['jobs'] if j['role']=='R3' and j['level']==3)
    ph=c.read_trace(d,j,'ph');ode=c.p2.load_ode(c.exp.COEFFICIENT);field=ode.medium(j,False);reduced,proof=c.equivalent_reference_field(field);q=c.query(.25);values=[];records=[]
    for label,f in [('original',field),('reduced',reduced)]:
        ref=c.p2.RadialReference(ph,f,ode,c.exp.COEFFICIENT,j['frequency_MHz'],5,105.);v=ref.query(q);values.append(v)
        np.savez_compressed(out/(label+'.npz'),s=q,xyz=v['xyz'],p=v['p'],group=v['group'])
        records.append({'label':label,'position_uncertainty_m':ref.max_pos_uncertainty,'group_uncertainty_km':ref.max_group_uncertainty,'warnings':ref.warnings})
        c.exp.write(out/(label+'_quad.json'),list(ref.uncertainty_records.values()))
    delta=c.metric(*values);limits={'rms_m':1e-7,'max_m':1e-7,'landing_m':1e-7,'group_rel':1e-12}
    c.exp.write(out/'report.json',{'job_id':j['job_id'],'proof':proof,'comparison':delta,'limits':limits,'uncertainty':records,'passed':all(delta[k]<=v for k,v in limits.items()),'source_sha256':{str(p.resolve()):c.common.sha(p) for p in [Path(__file__),Path(c.__file__),Path(c.p2.__file__),c.p2.preparation.ODE_SOURCE]}})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();verify(a.run.resolve())
