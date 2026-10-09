"""Read-only raw-MAT recheck; does not import the boundary analyzer."""
import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
from scipy.io import loadmat


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def interpolate(h, y, q):
    # The inputs have equal height spacing. Explicit Hermite basis is independent of SciPy's spline.
    step=np.diff(h)
    if not np.all(step==step[0]):raise ValueError('Expected equal frozen grid spacing')
    delta=float(step[0]);m=np.empty_like(y)
    m[1:-1]=(y[2:]-y[:-2])/(2*delta)
    m[0]=(-3*y[0]+4*y[1]-y[2])/(2*delta)
    m[-1]=(3*y[-1]-4*y[-2]+y[-3])/(2*delta)
    i=np.clip(np.searchsorted(h,q,side='right')-1,0,len(h)-2)
    u=(q-h[i])/delta
    return (2*u**3-3*u**2+1)*y[i]+(u**3-2*u**2+u)*delta*m[i]+(-2*u**3+3*u**2)*y[i+1]+(u**3-u**2)*delta*m[i+1]


def verify(run):
    d=run/'boundary';manifest=json.loads((d/'manifest.json').read_text())
    conclusion=json.loads((d/'diagnostics/conclusion.json').read_text())
    out=run/'boundary_independent_check'
    if out.exists():raise FileExistsError(out)
    out.mkdir();shutil.copy2(__file__,out/'source.py')
    errors=[];rows=[];mismatch=[];sources=[]
    for name,digest in manifest['source_sha256'].items():
        sources.append({'file':str(d/'source'/name),'sha256':sha(d/'source'/name),'matched':sha(d/'source'/name)==digest})
    for j in manifest['jobs']:
        path=d/'pharlap_raw'/(j['job_id']+'.mat');raw=loadmat(path,simplify_cells=True)
        grid=loadmat(j['grid_mat'],simplify_cells=True)
        v=raw['ray_state_vec']; rp=raw['ray_path_data']
        s=np.asarray(v['geometrical_path'])/1000
        xyz=np.column_stack([v[k] for k in ('pos_x','pos_y','pos_z')])/1000
        p=np.column_stack([v[k] for k in ('dir_x','dir_y','dir_z')])
        local=(s>=0)&(s<=20)
        h=np.linalg.norm(xyz,axis=1)-6371
        ne=interpolate(grid['height_km'],grid['ne_profile_m3'],h)
        output=np.asarray(rp['electron_density'])*1e6
        n2=np.asarray(rp['refractive_index'])**2
        p2=np.einsum('ij,ij->i',p,p)
        difference=abs(output[local]-np.maximum(ne[local],0))
        jid=dict(j);expected=jid.pop('job_sha256')
        digest=hashlib.sha256(json.dumps(jid,sort_keys=True,separators=(',',':')).encode()).hexdigest()
        checks={'job_identity':digest==expected==raw['job_sha256'],
                'grid_identity':sha(j['grid_mat'])==j['grid_sha256']==raw['grid_sha256'],
                'exporter_identity':raw['exporter_sha256']==j['exporter_sha256']==sha(d/'source/a4_c1m_uniform_export.m'),
                'full_twenty_km_coverage':bool(s[0]==0 and s[-1]>=20 and np.all(np.diff(s)>0)),
                'initial_p_exact':bool(np.array_equal(p[0],j['initial_p'])),
                'initial_position':bool(np.linalg.norm(xyz[0]-j['initial_xyz_km'])<1e-9),
                'density_output_clipped_hermite':bool(np.max(difference)<=max(1e-3,1e-12*np.max(abs(ne[local])))),
                'uniform_input':bool(np.array_equal(grid['Ne_m3'],np.broadcast_to(grid['ne_profile_m3'],grid['Ne_m3'].shape))),
                'zero_B_collision':bool(all(np.all(grid[k]==0) for k in ('Bx_T','By_T','Bz_T','collision_freq_Hz')))}
        for key,ok in checks.items():
            if not ok:errors.append({'job_id':j['job_id'],'check':key})
        if j['level']==3:
            for index in np.flatnonzero(local&(ne<0)&(abs(p2-n2)>1e-13)):
                mismatch.append({'job_id':j['job_id'],'index':int(index),'role':j['role'],
                                 'raw_ne_m3':float(ne[index]),'p2_minus_one':float(p2[index]-1),
                                 'output_n2_minus_one':float(n2[index]-1),
                                 'next_p2_minus_one':float(p2[index+1]-1),
                                 'next_is_reset':bool(abs(p2[index+1]-1)<=1e-13)})
        rows.append({'job_id':j['job_id'],'checks':checks,'raw_sha256':sha(path),
                     'max_density_output_error_m3':float(np.max(difference)),
                     'local_native_nodes':int(sum(local)),'raw_s_end_km':float(s[-1]),
                     'mex_sha256':raw['mex_sha256']})
    mex=Path('D:/pharlap/pharlap_4.7.4/mex/raytrace_3d_sp.mexw64')
    mex_ok=all(row['mex_sha256']==sha(mex) for row in rows)
    reported=json.loads((d/'diagnostics/state_reset_evidence.json').read_text())
    identities=lambda x:sorted((row['job_id'],row['index']) for row in x)
    check_summary={'all_sources_unchanged':all(x['matched'] for x in sources),'mex_identity':mex_ok,
                   'all_raw_case_checks':not errors,'same_thirty_native_mismatches':len(mismatch)==30 and identities(mismatch)==identities(reported),
                   'twelve_holdout_mismatches':sum(r['role']=='holdout' for r in mismatch)==12,
                   'mismatches_have_vacuum_output':all(abs(r['output_n2_minus_one'])<=1e-13 for r in mismatch),
                   'retains_unresolved_mechanism':conclusion['boundary_convention_identified'] is False and conclusion['production_repair_permitted'] is False}
    report={'schema':'p31-boundary-independent-check-v1','checks':check_summary,'passed':all(check_summary.values()),
            'job_count':len(rows),'state_output_mismatch_count':len(mismatch),
            'holdout_mismatch_count':sum(r['role']=='holdout' for r in mismatch),
            'next_sample_reset_count':sum(r['next_is_reset'] for r in mismatch),
            'sources':sources,'cases':rows,'errors':errors,'mismatches':mismatch,
            'propagation_rule_identified':False,'source_sha256':sha(Path(__file__)),
            'limits':'Finite native observations do not identify closed-MEX reset location or continuous internal gradient.'}
    (out/'verification.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    with (out/'mismatch_points.csv').open('w',newline='',encoding='utf-8') as stream:
        w=csv.DictWriter(stream,fieldnames=list(mismatch[0]));w.writeheader();w.writerows(mismatch)
    print(json.dumps({'passed':report['passed'],'checks':check_summary,'errors':errors}))
    if not report['passed']:raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--run',required=True,type=Path)
    verify(parser.parse_args().run)
