"""Independent spherical B=0 full-path quadrature, with natural endpoints.

The medium is the supplied piecewise polynomial, not a new fitted profile.
No production ray solver is imported. All uncertainty reports are estimates.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import CubicHermiteSpline, PPoly
from scipy.optimize import brentq


LIMITS = {'rms_m': 10., 'max_m': 50., 'landing_m': 50., 'group_rel': 1e-5}


class RadialFullPath:
    """Single ascending turn and return to the initial radius.

    field: h (km), coefficients (SciPy PPoly descending powers, m^-3).
    p0 is ck/omega. Clipping is an explicit hypothesis, never the default.
    """

    def __init__(self, field, x0, p0, coefficient, frequency_MHz,
                 *, radius=6371., density_policy='strict', eps=2e-10):
        self.h = np.asarray(field['h'], float)
        self.c = np.asarray(field['coefficients'], float)
        if (self.h.ndim != 1 or self.c.shape != (4, len(self.h)-1)
                or not np.isfinite(self.h).all() or not np.all(np.diff(self.h) > 0)
                or not np.isfinite(self.c).all()):
            raise ValueError('Invalid frozen polynomial')
        if density_policy not in ('strict', 'clip_zero'):
            raise ValueError('Unknown density policy')
        self.policy = density_policy
        self.poly = PPoly(self.c, self.h, extrapolate=False)
        self.radius = float(radius)
        self.alpha = float(coefficient)/(float(frequency_MHz)*1e6)**2
        self.eps = eps
        self.x0 = np.asarray(x0, float)
        self.p0 = np.asarray(p0, float)
        self.r0 = float(np.linalg.norm(self.x0))
        self.er = self.x0/self.r0
        radial = float(self.p0 @ self.er)
        pt = self.p0-radial*self.er
        if radial <= 0 or np.linalg.norm(pt) == 0:
            raise ValueError('Reference requires upward nonvertical launch')
        self.et = pt/np.linalg.norm(pt)
        self.J = self.r0*np.linalg.norm(pt)
        self.E = 1.
        if abs(float(self.p0 @ self.p0)-self.n2(self.r0)) > 1e-10:
            raise ValueError('Initial dispersion mismatch')
        # Preserve the actual initial invariant, including binary64 roundoff.
        self.E = float(self.p0 @ self.p0)+(1.-self.n2(self.r0))
        extrema = self.poly.derivative().roots(extrapolate=False)
        check = np.r_[self.h, extrema[np.isfinite(extrema)]]
        if self.policy == 'strict' and np.min(self.poly(check)) < 0:
            raise ValueError('Negative density in supplied polynomial')
        zeros = self.poly.roots(extrapolate=False)
        self.breaks = np.unique(np.r_[self.h, zeros[np.isfinite(zeros)]])
        self.rt = self._turn()
        self.umax = math.sqrt(self.rt-self.r0)
        cuts = self.radius+self.breaks
        self.ucuts = np.sort(np.r_[0., np.sqrt(self.rt-cuts[(cuts > self.r0) & (cuts < self.rt)]), self.umax])
        self.quad_errors = []
        segments = [self._integral(a, b) for a,b in zip(self.ucuts[:-1], self.ucuts[1:])]
        self.half_quad_errors = np.sum(self.quad_errors,axis=0)
        self.cumulative = np.vstack([np.zeros(3), np.cumsum(segments, axis=0)])
        self.half = self.cumulative[-1]
        self.s_end = float(2*self.half[0])
        self.group_end = float(2*self.half[2])
        self.endpoint = self.r0*(math.cos(2*self.half[1])*self.er+math.sin(2*self.half[1])*self.et)

    def n2(self, r):
        v = float(self.poly(r-self.radius))
        if self.policy == 'clip_zero':
            v = max(0., v)
        return self.E-self.alpha*v

    def dn2(self, r):
        h = r-self.radius
        if self.policy == 'clip_zero' and float(self.poly(h)) <= 0:
            return 0.
        return -self.alpha*float(self.poly.derivative()(h))

    def _turn(self):
        q = lambda r: self.n2(r)*r*r-self.J*self.J
        samples = np.unique(np.r_[self.r0, self.radius+self.breaks[self.breaks > self.r0-self.radius],
                                  np.linspace(self.r0, self.radius+self.h[-1], 2001)])
        values = np.array([q(r) for r in samples])
        ix = np.flatnonzero((values[:-1] > 0) & (values[1:] <= 0))
        if not len(ix):
            raise ValueError('No ascending turn: penetration is not a full reflected ray')
        i = int(ix[0])
        guess = brentq(q, samples[i], samples[i+1], xtol=1e-11)
        self.turn_segment = min(len(self.h)-2, int(np.searchsorted(self.h, guess-self.radius)-1))
        with mp.workdps(80):
            k = self.turn_segment
            c = [mp.mpf(float(v)) for v in self.c[:, k]]
            base = mp.mpf(self.radius)+mp.mpf(float(self.h[k]))
            alpha, J, E = mp.mpf(self.alpha), mp.mpf(self.J), mp.mpf(self.E)
            def exact_q(r):
                t = r-base
                ne = ((c[0]*t+c[1])*t+c[2])*t+c[3]
                return (E-alpha*ne)*r*r-J*J
            root = mp.findroot(exact_q, (mp.mpf(samples[i]), mp.mpf(samples[i+1])), solver='secant')
            if exact_q(root+mp.mpf('1e-6')) >= 0:
                raise ValueError('Degenerate or wrong turning root')
            self.turn_taylor = np.array([float(mp.diff(exact_q, root, k)/mp.factorial(k)) for k in range(1,6)])
            self.turn_root_residual = str(abs(exact_q(root)))
            return float(root)

    def _components(self, u):
        r = self.rt-u*u
        if r >= self.radius+self.h[self.turn_segment]:
            delta = -u*u
            q = delta*np.polynomial.polynomial.polyval(delta, self.turn_taylor)
        else:
            q = self.n2(r)*r*r-self.J*self.J
        if u == 0:
            factor = 2/math.sqrt(-self.turn_taylor[0])
        else:
            if q <= 0:
                raise ValueError('Nonpositive radicand before selected turn')
            factor = 2*u/math.sqrt(q)
        n = math.sqrt(self.n2(r))
        return factor*np.array([n*r, self.J/r, r])

    def _integral(self, a, b):
        vals, errs = [], []
        for k in range(3):
            v,e = quad(lambda u: self._components(u)[k], a,b,
                       epsabs=self.eps, epsrel=self.eps, limit=150)
            vals.append(v); errs.append(e)
        self.quad_errors.append(errs)
        return np.array(vals)

    def _from_turn(self, u):
        i = min(len(self.ucuts)-2, max(0, int(np.searchsorted(self.ucuts,u)-1)))
        return self.cumulative[i]+self._integral(self.ucuts[i],u)

    def query(self, s):
        s = np.atleast_1d(np.asarray(s,float))
        if np.any(s < 0) or np.any(s > self.s_end+1e-8):
            raise ValueError('Query outside natural trajectory')
        xyz, group, tangent, p = [], [], [], []
        for distance in s:
            down = distance > self.half[0]
            target = abs(distance-self.half[0])
            if target < 1e-11:
                u = 0.
            elif abs(target-self.half[0]) < 1e-10:
                u = self.umax
            else:
                i = min(len(self.ucuts)-2, max(0, int(np.searchsorted(self.cumulative[:,0],target)-1)))
                u = brentq(lambda v: self._from_turn(v)[0]-target,
                           self.ucuts[i],self.ucuts[i+1],xtol=2e-12)
            integrals = self._from_turn(u)
            sign = 1 if down else -1
            theta = self.half[1]+sign*integrals[1]
            r = self.rt-u*u
            er = math.cos(theta)*self.er+math.sin(theta)*self.et
            et = -math.sin(theta)*self.er+math.cos(theta)*self.et
            n = math.sqrt(self.n2(r))
            tangential = self.J/(n*r)
            radial = (-1 if down else 1)*math.sqrt(max(0,1-tangential*tangential))
            direction = radial*er+tangential*et
            xyz.append(r*er); tangent.append(direction); p.append(n*direction)
            group.append(self.half[2]+sign*integrals[2])
        return {'s':s,'xyz':np.asarray(xyz),'group':np.asarray(group),
                'tangent':np.asarray(tangent),'p':np.asarray(p),
                'termination':'ground','natural_s_end':self.s_end,
                'natural_group_end':self.group_end,'natural_endpoint':self.endpoint.copy()}


def compare_natural(a, b):
    """Never replace independent natural endpoints by common-prefix endpoints."""
    if a['termination'] != 'ground' or b['termination'] != 'ground':
        raise ValueError('Early termination or non-ground topology')
    if not np.array_equal(a['s'],b['s']):
        raise ValueError('Different absolute arc-length sample grids')
    if (len(a['s']) < 2 or abs(float(a['s'][0])) > 1e-10
            or float(a['s'][-1]) < min(a['natural_s_end'],b['natural_s_end'])-1e-8):
        raise ValueError('Truncated common prefix cannot certify full trajectories')
    d = np.linalg.norm(a['xyz']-b['xyz'],axis=1)*1000
    return {'rms_m':float(np.sqrt(np.mean(d*d))), 'max_m':float(d.max()),
            'landing_m':float(np.linalg.norm(a['natural_endpoint']-b['natural_endpoint'])*1000),
            'group_rel':float(abs(a['natural_group_end']-b['natural_group_end'])/abs(b['natural_group_end']))}


def cartesian_control(reference, rtol=2e-12, max_step=.5):
    """Independent analytic-gradient ODE; no AH production helper or finite difference."""
    def rhs(s,y):
        r = np.linalg.norm(y[:3]); h = r-reference.radius
        # Event bracketing stages can cross the ground; the returned path ends at it.
        ne = float(reference.poly(h,extrapolate=True))
        gradient = float(reference.poly.derivative()(h,extrapolate=True))
        if reference.policy == 'clip_zero' and ne <= 0:
            ne,gradient = 0.,0.
        n = math.sqrt(reference.E-reference.alpha*ne)
        return np.r_[y[3:6]/n, -reference.alpha*gradient*y[:3]/(2*n*r), 1/n]
    def ground(s,y):
        return np.linalg.norm(y[:3])-reference.r0 if s > 1e-7 else 1e-7
    ground.terminal=True; ground.direction=-1
    result = solve_ivp(rhs,(0,reference.s_end*1.05),np.r_[reference.x0,reference.p0,0.],
                       method='DOP853',rtol=rtol,atol=rtol*.01,max_step=max_step,
                       events=ground,dense_output=True)
    if not len(result.t_events[0]):
        raise ValueError('Cartesian control did not land')
    return result


def analytic_control_report(output=None):
    C,f,R = 80.61638494312905,10.,6371.
    h = np.arange(0.,462.,2.)
    a = .002
    slope = a*(f*1e6)**2/C
    c = np.vstack([np.zeros(len(h)-1),np.zeros(len(h)-1),np.full(len(h)-1,slope),slope*h[:-1]])
    field={'h':h,'coefficients':c}
    x0=np.array([R,0.,0.]);p0=np.array([.5,math.sqrt(.75),0.])
    ref=RadialFullPath(field,x0,p0,C,f)
    other=RadialFullPath(field,x0,p0,C,f,eps=2e-12)
    ode=cartesian_control(ref)
    common_end=min(ref.s_end,float(ode.t[-1]))
    # Pointwise reference queries also examine off-node positions, including the turn.
    q=np.unique(np.r_[np.linspace(0,common_end,401),ref.half[0]])
    exact=ref.query(q); y=ode.sol(q)
    control={'s':q,'xyz':y[:3].T,'group':y[6], 'termination':'ground',
             'natural_s_end':float(ode.t[-1]),'natural_group_end':float(ode.y[6,-1]),
             'natural_endpoint':ode.y[:3,-1]}
    agreement=compare_natural(exact,control)
    refinement=compare_natural(exact,other.query(q))
    native=ref.query(ode.t[ode.t <= ref.s_end])
    nn=len(native['s'])
    native_error=float(np.max(np.linalg.norm(native['xyz']-ode.y[:3,:nn].T,axis=1))*1000)
    spline=CubicHermiteSpline(ode.t,ode.y[:3].T,
                             ode.y[3:6].T/np.sqrt([ref.n2(np.linalg.norm(v)) for v in ode.y[:3].T])[:,None])
    # Compare reconstruction of exact reference nodes separately from native propagation.
    ref_spline=CubicHermiteSpline(native['s'],native['xyz'],native['tangent'])
    ref_group_spline=CubicHermiteSpline(native['s'],native['group'],1/np.linalg.norm(native['p'],axis=1))
    mids=(native['s'][:-1]+native['s'][1:])/2
    midpoint_exact=ref.query(mids)
    point_error=float(np.max(np.linalg.norm(ref_spline(mids)-midpoint_exact['xyz'],axis=1))*1000)
    point_group_error=float(np.max(abs(ref_group_spline(mids)-midpoint_exact['group']))/ref.group_end)
    if output is not None:
        np.savez(output/'analytic_control_raw.npz', height=h,coefficients=c,
                 ode_s=ode.t,ode_state=ode.y,reference_s=q,reference_xyz=exact['xyz'],
                 reference_group=exact['group'],native_reference_s=native['s'],
                 native_reference_xyz=native['xyz'],native_reference_group=native['group'],
                 midpoint_s=mids,midpoint_reference_xyz=midpoint_exact['xyz'],
                 midpoint_reference_group=midpoint_exact['group'])
        import csv
        with (output/'analytic_control_points.csv').open('w',newline='',encoding='utf-8') as stream:
            writer=csv.writer(stream);writer.writerow(['s_km','reference_height_km','position_error_m','group_error_km'])
            writer.writerows(zip(q,np.linalg.norm(exact['xyz'],axis=1)-R,
                                 np.linalg.norm(exact['xyz']-control['xyz'],axis=1)*1000,
                                 exact['group']-control['group']))
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(2,1,figsize=(8,6),sharex=True)
        axes[0].plot(q,np.linalg.norm(exact['xyz'],axis=1)-R);axes[0].set_ylabel('Height (km)')
        axes[1].plot(q,np.linalg.norm(exact['xyz']-control['xyz'],axis=1)*1000)
        axes[1].set_ylabel('Cartesian control error (m)');axes[1].set_xlabel('Absolute geometric path (km)')
        fig.tight_layout();fig.savefig(output/'analytic_control.png',dpi=140);plt.close(fig)
    steps=[]
    for ds in (1.,.5,.25):
        grid=np.linspace(0,common_end,int(math.ceil(common_end/ds))+1)
        truth=ref.query(grid)
        reconstructed=dict(truth,xyz=spline(grid),natural_endpoint=control['natural_endpoint'],
                           natural_group_end=control['natural_group_end'])
        steps.append({'step_km':ds,'metric':compare_natural(reconstructed,truth)})
    endpoint_bad=dict(exact,natural_endpoint=exact['natural_endpoint']+np.array([.1,0,0]))
    group_bad=dict(exact,natural_group_end=exact['natural_group_end']*1.001)
    rejected=False
    try: compare_natural(dict(exact,termination='early'),exact)
    except ValueError: rejected=True
    truncated_rejected=False
    try: compare_natural(dict(exact,s=exact['s'][:-1],xyz=exact['xyz'][:-1]),
                         dict(exact,s=exact['s'][:-1],xyz=exact['xyz'][:-1]))
    except ValueError: truncated_rejected=True
    t=np.array([0.,7.,20.]);xyz=np.c_[t,np.zeros((3,2))];v=np.tile([1.,0.,0.],(3,1))
    straight=CubicHermiteSpline(t,xyz,v)(np.linspace(0,20,31))
    tests={'straight_reconstruction':bool(np.max(abs(straight[:,0]-np.linspace(0,20,31)))<1e-12),
           'cartesian_curved_control':bool(all(agreement[k]<.2*LIMITS[k] for k in LIMITS)),
           'quadrature_refinement':bool(all(refinement[k]<.01*.2*LIMITS[k] for k in LIMITS)),
           'natural_endpoint_mismatch_detected':bool(compare_natural(endpoint_bad,exact)['landing_m']>50),
           'natural_group_mismatch_detected':bool(compare_natural(group_bad,exact)['group_rel']>1e-5),
           'early_termination_rejected':rejected,
           'truncated_common_prefix_rejected':truncated_rejected,
           'turn_point_finite':bool(np.isfinite(ref.query([ref.half[0]])['xyz']).all()),
           'point_reconstruction_small':bool(point_error<2.),
           'point_group_reconstruction_small':bool(point_group_error<2e-6)}
    return {'schema':'p31-fullpath-independent-controls-v1','checks':tests,'passed':all(tests.values()),
            'fullpath_reference_implementation_verified':all(tests.values()),
            'formal_smooth_fullpath_verified':False,'formal_blocker':'Boundary propagation convention not yet identified',
            'scope':'Analytic affine density, B=0, spherical single ascending turn',
            'reference':{'C':C,'frequency_MHz':f,'radius_km':R,'density_policy':ref.policy,
                         'initial_invariant_E':ref.E,'s_end_km':ref.s_end,'group_end_km':ref.group_end,'apex_km':ref.rt-R,
                         'root_residual_80dps':ref.turn_root_residual,'quadrature_interval_certified':False,
                         'half_quad_error_estimate_ds_theta_group':ref.half_quad_errors.tolist()},
            'cartesian_agreement':agreement,'quadrature_refinement':refinement,
            'native_propagation_max_m':native_error,'exact_node_midpoint_reconstruction_max_m':point_error,
            'exact_node_midpoint_group_reconstruction_rel':point_group_error,
            'sampling':steps,'sampling_deltas':[{k:abs(steps[i]['metric'][k]-steps[i-1]['metric'][k]) for k in LIMITS} for i in (1,2)]}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output-name',default='reference_fullpath')
    args=parser.parse_args();out=args.run/args.output_name
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True)
    shutil.copy2(__file__,out/'source.py')
    report=analytic_control_report(out)
    report['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (out/'controls.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps({'passed':report['passed'],'output':str(out)},ensure_ascii=True),flush=True)
    if not report['passed']:raise SystemExit(1)


if __name__=='__main__':main()
