"""Nonlinear C2 piecewise-polynomial control for the full reflected reference."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
from _a4_c1m_p31_fullpath import RadialFullPath, cartesian_control, compare_natural, LIMITS


def main(run):
    out=run/'reference_fullpath_piecewise'
    if out.exists():raise FileExistsError(out)
    out.mkdir();shutil.copy2(__file__,out/'source.py')
    reference=Path(__file__).with_name('_a4_c1m_p31_fullpath.py')
    shutil.copy2(reference,out/'reference_source.py')
    C,f,R=80.61638494312905,10.,6371.
    h=np.arange(0.,462.,2.);v=h[:-1];active=v>=80;u=np.maximum(v-80,0)
    a,b=4e-6,1e-8
    coeff=(f*1e6)**2/C*np.vstack([b*active,a+3*b*u,2*a*v+3*b*u*u,a*v*v+b*u**3])
    cases=[]
    for angle in (22.,39.5):
        p=np.array([np.sin(np.deg2rad(angle)),np.cos(np.deg2rad(angle)),0.])
        ref=RadialFullPath({'h':h,'coefficients':coeff},[R,0,0],p,C,f)
        ode=cartesian_control(ref)
        q=np.unique(np.r_[np.linspace(0,min(ref.s_end,ode.t[-1]),257),ref.half[0]])
        exact=ref.query(q);y=ode.sol(q)
        raw=dict(s=q,xyz=y[:3].T,termination='ground',natural_s_end=float(ode.t[-1]),
                 natural_group_end=float(ode.y[6,-1]),natural_endpoint=ode.y[:3,-1])
        metric=compare_natural(exact,raw)
        np.savez(out/f'e{angle:g}.npz',height=h,coefficients=coeff,s=q,xyz=exact['xyz'],
                 group=exact['group'],cartesian_state=y,cartesian_native_s=ode.t,cartesian_native_state=ode.y)
        cases.append({'elevation_deg':angle,'apex_km':ref.rt-R,'s_end_km':ref.s_end,'group_end_km':ref.group_end,
                      'cartesian_agreement':metric,'passed':all(metric[k]<.2*LIMITS[k] for k in LIMITS),
                      'polynomial':'Ne=f_Hz^2/C * (4e-6*h^2 + 1e-8*max(h-80,0)^3)',
                      'quad_error_estimate':ref.half_quad_errors.tolist()})
    report={'schema':'p31-fullpath-piecewise-control-v1','cases':cases,'passed':all(c['passed'] for c in cases),
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'reference_sha256':hashlib.sha256(reference.read_bytes()).hexdigest(),
            'formal_smooth_fullpath_verified':False}
    (out/'verification.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report))
    if not report['passed']:raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--run',required=True,type=Path)
    main(parser.parse_args().run)
