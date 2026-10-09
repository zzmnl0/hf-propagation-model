"""P2 fixed-100-km numerical attribution, executed from a frozen source copy."""
from __future__ import annotations
import os
for _key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_key] = "1"
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import multiprocessing
from pathlib import Path
import sys
import types
import warnings
import numpy as np
from scipy.integrate import quad, IntegrationWarning
from scipy.interpolate import PchipInterpolator, CubicHermiteSpline
from scipy.optimize import brentq
import _a4_c1m_p2_analysis as preparation
import _run_a4_c1m_attribution as common

METHODS = ("pchip", "derivative_hermite")
QUERY_STEPS = (1., .5, .25)
STRIDES = (4, 2, 1)
VALIDATION_STEP = .0625
QUAD_PROFILES = ((1e-5, 1e-5, 1e-3), (1e-8, 1e-8, 1e-6), (1e-11, 1e-11, 1e-9))
KEYS = tuple(preparation.LIMITS)


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def grid(step: float) -> np.ndarray:
    return np.linspace(0., 100., round(100./step)+1)


def errors(a: dict, b: dict) -> dict:
    m = common.norm_metrics(a["xyz"]-b["xyz"])
    return {"rms_m": m["rms_m"], "max_m": m["max_m"],
        "landing_m": float(np.linalg.norm(a["xyz"][-1]-b["xyz"][-1])*1000),
        "group_rel": abs(a["g100"]-b["g100"])/max(abs(a["g100"]),1e-15)}


def controls(values: list[dict]) -> dict:
    d12, d23 = errors(values[0],values[1]), errors(values[1],values[2])
    return {"delta12": d12, "delta23": d23,
        "monotonic": {k: d23[k] <= d12[k] for k in KEYS},
        "last_change_small": {k: d23[k] <= .2*preparation.LIMITS[k] for k in KEYS},
        "passed": all(d23[k] <= d12[k] and d23[k] <= .2*preparation.LIMITS[k] for k in KEYS),
        "roundoff_waiver_applied": False}


def load_ode(C: float):
    # Executing the immutable P1 text avoids creating bytecode in the sealed P1 tree.
    module = types.ModuleType("_p2_frozen_P1_ODE")
    module.__file__ = str(preparation.ODE_SOURCE)
    sys.modules[module.__name__] = module
    exec(compile(preparation.ODE_SOURCE.read_text(encoding="utf-8"),module.__file__,"exec"),module.__dict__)
    module.C = C
    return module


def prefix(trace: dict, stride: int=1) -> dict:
    s = trace["s"]
    right = int(np.searchsorted(s,100.,side="left"))
    if right >= len(s):
        raise ValueError("100-km interval is incomplete")
    indices = np.unique(np.r_[np.arange(0,right+1,stride),0,max(0,right-1),right]).astype(int)
    if np.any(np.diff(s[indices]) <= 0):
        raise ValueError("Invalid native geometric coordinate")
    result = dict(trace)
    for key in ("s","xyz","p","tangent","group"):
        if key in trace and trace[key] is not None:
            result[key] = np.asarray(trace[key])[indices]
    result["native_indices"] = indices
    return result


def splines(trace: dict, method: str):
    if method == "pchip":
        position = PchipInterpolator(trace["s"],trace["xyz"],axis=0)
        group = PchipInterpolator(trace["s"],trace["group"]) if trace.get("group") is not None else None
    else:
        position = CubicHermiteSpline(trace["s"],trace["xyz"],trace["tangent"],axis=0)
        group = CubicHermiteSpline(trace["s"],trace["group"],1/np.linalg.norm(trace["p"],axis=1)) if trace.get("group") is not None else None
    return position,group


def reconstruct(trace: dict, method: str, sample: np.ndarray) -> dict:
    position,group = splines(trace,method)
    return {"xyz": position(sample), "g100": float(group(100.)) if group is not None else trace["native_group100"]}


def read_ah(meta_path: Path, job: dict) -> dict:
    meta = preparation._json(meta_path)
    if meta["status"] != "completed":
        return {"status": meta["status"],"meta": meta}
    raw_path = meta_path.with_suffix(".npz")
    if meta["job_sha256"] != job["job_sha256"] or meta["grid_sha256"] != job["grid_sha256"] or common.sha(raw_path) != meta["raw_npz_sha256"]:
        raise ValueError("AH frozen identity mismatch")
    with np.load(raw_path) as raw:
        result = {"s": np.array(raw["geom_path_km"]),"xyz":np.array(raw["xyz_km"]),
            "p":np.array(raw["wavevector_km_inv"])/(2*np.pi*job["frequency_MHz"]*1e6/299792.458),
            "tangent":np.array(raw["group_tangent_ecef"]),"group":None,
            "native_group100":float(meta["group_path_km_native"]),"status":"completed","meta":meta}
    if abs(result["s"][-1]-100.) > 1e-10 or np.any(np.diff(result["s"]) <= 0):
        raise ValueError("AH endpoint is not the complete native 100-km endpoint")
    return result


class RadialReference:
    """Spherical Snell invariant quadratures, with no fitted quantities."""
    def __init__(self,ph: dict,field: dict,ode,C: float,frequency: float,level: int,max_s: float):
        self.r0 = np.array(ph["xyz"][0]); self.p0 = np.array(ph["p"][0]); self.radius0=float(np.linalg.norm(self.r0))
        self.field,self.ode=field,ode; self.alpha=C/(frequency*1e6)**2
        self.E=float(self.p0@self.p0+self.alpha*ode.density_gradient(self.radius0-6371.,field)[0])
        self.L=float(np.linalg.norm(np.cross(self.r0,self.p0)))
        self.er=self.r0/self.radius0
        tangential=self.p0-(self.p0@self.er)*self.er
        self.et=tangential/np.linalg.norm(tangential)
        self.epsabs,self.epsrel,self.xtol=QUAD_PROFILES[level-1]
        self.cache={}; self.warnings=[]; self.max_pos_uncertainty=0.; self.max_group_uncertainty=0.
        fraction=float((self.p0@self.er)/np.linalg.norm(self.p0))
        self.upper=self.radius0+max_s*fraction+5.
        for _ in range(40):
            if self._integral(self.upper,0)[0] >= max_s:break
            self.upper+=.5
        else:raise ValueError("Independent radial arc cannot be bracketed")
    def local(self,r: float):
        Ne=self.ode.density_gradient(r-6371.,self.field)[0]
        n2=self.E-self.alpha*Ne; q=n2-(self.L/r)**2
        if n2<=0 or q<=0:raise ValueError("Independent radial interval reaches a turning/cutoff")
        return np.sqrt(n2),np.sqrt(q)
    def _integral(self,r: float,kind: int):
        if r==self.radius0:return 0.,0.
        points=(6371.+self.field["h"])
        points=points[(points>self.radius0)&(points<r)]
        def f(radius):
            n,q=self.local(radius)
            return (n/q,self.L/(radius**2*q),1/q)[kind]
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always",IntegrationWarning)
            value,error=quad(f,self.radius0,r,points=points,epsabs=self.epsabs,epsrel=self.epsrel,limit=max(80,len(points)+20))
        self.warnings.extend(str(w.message) for w in caught)
        return value,error
    def query(self,sample: np.ndarray) -> dict:
        rows=[]
        for sample_s in np.asarray(sample,float):
            key=float(sample_s)
            if key not in self.cache:
                if key==0.:
                    self.cache[key]=(self.r0.copy(),self.p0.copy(),0.,self.p0/np.linalg.norm(self.p0))
                else:
                    radius=brentq(lambda r:self._integral(r,0)[0]-key,self.radius0,self.upper,xtol=self.xtol,rtol=4*np.finfo(float).eps)
                    phi,phi_error=self._integral(radius,1); G,G_error=self._integral(radius,2)
                    _,s_error=self._integral(radius,0)
                    n,q=self.local(radius)
                    er=np.cos(phi)*self.er+np.sin(phi)*self.et
                    et=-np.sin(phi)*self.er+np.cos(phi)*self.et
                    p=q*er+self.L/radius*et
                    dr_error=self.xtol+s_error
                    position_error=dr_error*(1+self.L/(radius*q))+radius*phi_error
                    group_error=G_error+dr_error/q
                    self.max_pos_uncertainty=max(self.max_pos_uncertainty,position_error*1000.)
                    self.max_group_uncertainty=max(self.max_group_uncertainty,group_error)
                    self.cache[key]=(radius*er,p,G,p/n)
            rows.append(self.cache[key])
        xyz=np.array([r[0] for r in rows]); p=np.array([r[1] for r in rows]); G=np.array([r[2] for r in rows])
        return {"xyz":xyz,"p":p,"group":G,"tangent":np.array([r[3] for r in rows]),"g100":float(self.cache[100.][2]) if 100. in self.cache else float(G[-1])}


def field_bounds(field: dict,alpha: float,E: float,max_s: float) -> dict:
    lower,upper=100.-1e-8,100.+max_s+1e-8
    ne_max,grad_max=0.,0.
    for i in range(len(field["h"])-1):
        lo=max(lower,field["h"][i]); hi=min(upper,field["h"][i+1])
        if hi<lo:continue
        a,b,c,d=field["coefficients"][:,i]; u,v=lo-field["h"][i],hi-field["h"][i]
        candidates=[u,v]
        roots=np.roots([3*a,2*b,c]) if a!=0 else ([-c/(2*b)] if b!=0 else [])
        candidates += [float(x.real) for x in roots if np.isreal(x) and u<=x.real<=v]
        ne_max=max(ne_max,*[float(((a*x+b)*x+c)*x+d) for x in candidates])
        slope_points=[u,v]+([-b/(3*a)] if a!=0 and u<=-b/(3*a)<=v else [])
        grad_max=max(grad_max,*[abs(float((3*a*x+2*b)*x+c)) for x in slope_points])
    n2_min=E-alpha*ne_max
    if n2_min<=0:raise ValueError("No positive radial n2 bound")
    return {"alpha":alpha,"height_min_km":lower,"height_max_km":upper,"Ne_max_m3":ne_max,
        "Neprime_max_m3_per_km":grad_max,"n2_min":n2_min,
        "position_M2_per_km":alpha*grad_max/(2*n2_min),"group_M2_per_km":alpha*grad_max/(2*n2_min**1.5),
        "proof":"t'=-alpha*(gradNe-(gradNe dot t)t)/(2*p^2); G''=alpha*(gradNe dot t)/(2*|p|^3). Piecewise Hermite Ne is C1; only bounded second derivatives are used, across knots as well.",
        "h4_remainder_used":False}


def polynomial_second_max(spline) -> float:
    maxima=[]
    for i in range(len(spline.x)-1):
        if spline.x[i]>100.:continue
        width=min(100.,spline.x[i+1])-spline.x[i]
        for x in (0.,width):
            value=6*spline.c[0,i]*x+2*spline.c[1,i]
            maxima.append(float(np.linalg.norm(np.atleast_1d(value))))
    return max(maxima,default=0.)


def reconstruction_bound(raw: dict,method: str,ref: RadialReference,proof: dict,side: str) -> dict:
    trace=prefix(raw); exact=ref.query(trace["s"])
    model={"s":trace["s"],"xyz":exact["xyz"],"p":exact["p"],"tangent":exact["tangent"],"group":exact["group"]}
    pos,G=splines(model,method)
    q=grid(VALIDATION_STEP); truth=ref.query(q)
    actual=np.linalg.norm(pos(q)-truth["xyz"],axis=1)*1000.
    derivative=np.linalg.norm(pos.derivative()(q)-truth["tangent"],axis=1)
    second=polynomial_second_max(pos)
    fill=VALIDATION_STEP/2*float(np.max(derivative))+VALIDATION_STEP**2/8*(second+proof["position_M2_per_km"])
    uncertainty=8*ref.max_pos_uncertainty
    bound=float(np.max(actual))+fill*1000.+uncertainty
    endpoint=float(actual[-1])+uncertainty
    group_error=np.abs(G(q)-truth["group"])
    group_derivative=np.abs(G.derivative()(q)-1/np.linalg.norm(truth["p"],axis=1))
    group_second=polynomial_second_max(G)
    group_fill=VALIDATION_STEP/2*float(np.max(group_derivative))+VALIDATION_STEP**2/8*(group_second+proof["group_M2_per_km"])
    group_uncertainty=8*ref.max_group_uncertainty
    group100=float(group_error[-1])+group_uncertainty
    denominator=abs(float(reconstruct(trace,method,np.array([100.]))["g100"]))
    raw_query=reconstruct(trace,method,q)
    node_propagation=np.linalg.norm(raw_query["xyz"]-pos(q),axis=1)*1000.
    full_proof=dict(proof,C=ref.alpha*1e14,R_km=6371.,node_maxgap_km=float(np.max(np.diff(trace["s"]))),
        validation_step_km=VALIDATION_STEP,poly_second_derivative_max=second,fill_bound_m=fill*1000.,
        quad_uncertainty_m=uncertainty,actual_reference_interpolation_rms_m=float(np.sqrt(np.mean(actual**2))),
        actual_reference_interpolation_max_m=float(np.max(actual)),actual_reference_interpolation_endpoint_m=float(actual[-1]),
        group_poly_second_derivative_max=group_second,group_fill_bound_km=group_fill,group_quad_uncertainty_km=group_uncertainty,
        actual_reference_group_interpolation_max_km=float(np.max(group_error)),
        actual_reference_group_interpolation_endpoint_km=float(group_error[-1]),
        propagated_node_error_query_rms_m=float(np.sqrt(np.mean(node_propagation**2))),
        propagated_node_error_query_max_m=float(np.max(node_propagation)),node_error_in_reconstruction_budget=False,
        reference_samples_are_exact_model_queries=True,quadrature_errors_are_numerical_estimates=True,
        continuous_fill_proof="Nearest uniform validation point is <=step/2 away. Taylor with bounded e'' gives max(sample |e|)+step/2*max(sample |e'|)+step^2/8*(max |I''|+M2).")
    return {"bounds":{"rms_m":bound,"max_m":bound,"landing_m":endpoint,"group_rel":group100/denominator if side=="ph" else 0.},
        "proof":full_proof,"bound_established":True,
        "support_method":"Physical radial reference sampled at original native knots, exact derivatives, C2 interpolation error plus continuous fill bound; native propagation errors excluded",
        "AH_group_is_native_total":side=="ah"}


def analytic_checks() -> dict:
    results=[];q=grid(.0625)
    for case in ("straight","circle"):
        def exact(s):
            if case=="straight":return np.c_[s,np.zeros_like(s),np.zeros_like(s)],np.tile([1.,0.,0.],(len(s),1))
            angle=s/200.;return np.c_[200*np.sin(angle),200*(1-np.cos(angle)),np.zeros_like(s)],np.c_[np.cos(angle),np.sin(angle),np.zeros_like(s)]
        target,_=exact(q)
        for method in METHODS:
            values=[]
            for step in (2.,1.,.5):
                nodes=grid(step);xyz,tangent=exact(nodes)
                trace={"s":nodes,"xyz":xyz,"tangent":tangent,"p":tangent,"group":nodes}
                reconstructed=reconstruct(trace,method,q)
                error=np.linalg.norm(reconstructed["xyz"]-target,axis=1)*1000.
                values.append({"step_km":step,"rms_m":float(np.sqrt(np.mean(error**2))),"max_m":float(np.max(error)),"group_max_km":float(np.max(abs(splines(trace,method)[1](q)-q)))})
            results.append({"case":case,"method":method,"levels":values,
                "error_monotonic":all(values[2][k]<=values[1][k]<=values[0][k] for k in ("rms_m","max_m","group_max_km")),
                "implementation_passed": all(v["max_m"]<1. for v in values),"analytic_control_only":True})
    return {"schema":"a4-c1m-p2-analytic-controls-v1","cases":results,"no_raw_trajectory_bound_inferred_from_controls":True,
        "passed":all(r["implementation_passed"] for r in results)}


def candidate_worker(job,out,C,hashes):
    return common._a2_ah_one(job,out,C,"hermite",.125,hashes)


def invariant_rows(trace: dict,ode,field,C: float,freq: float):
    density=np.array([ode.density_gradient(h,field)[0] for h in np.linalg.norm(trace["xyz"],axis=1)-6371.])
    E=np.sum(trace["p"]**2,axis=1)+C/(freq*1e6)**2*density
    L=np.cross(trace["xyz"],trace["p"]);L0=L[0]
    return {"energy_initial":float(E[0]),"energy_max_drift":float(np.max(abs(E-E[0]))),
        "angular_momentum_vector_max_relative_drift":float(np.max(np.linalg.norm(L-L0,axis=1))/np.linalg.norm(L0)),
        "angular_momentum_initial":L0.tolist(),"initial_position_km":trace["xyz"][0].tolist(),"initial_p":trace["p"][0].tolist()}


def analyze(run: Path,output: Path,manifest: dict) -> None:
    C=manifest["constants"]["C"]
    ode=load_ode(C);jobs=preparation._json(run/"inputs/internal_jobs.json")["primary"]
    preparation._write(output/"analytic_reconstruction.json",analytic_checks())
    sources={"C":C,"P1_ODE_sha256":common.sha(preparation.ODE_SOURCE),"producer_sha256":common._a2_ah_source_hashes(),
        "candidate_workers":2,"BLAS_threads":1,"DOP_profiles":ode.PROFILES,"quadrature_profiles":QUAD_PROFILES,
        "C_override_is_in_separate_reference_module":True,"production_K_FP_modified_in_analysis_process":False}
    preparation._write(output/"numerical_sources.json",sources)
    candidate=run/"ah/fd0.125";candidate.mkdir(parents=True)
    hashes=common._a2_ah_source_hashes(); worker_rows=[]
    with ProcessPoolExecutor(max_workers=2,mp_context=multiprocessing.get_context("spawn")) as pool:
        futures={pool.submit(common._a2_ah_one,item["original_job"],str(candidate),C,"hermite",.125,hashes):item for item in jobs}
        for future in as_completed(futures):
            worker_rows.append(future.result())
            preparation._write(output/"AH_candidate_checkpoint.json",{"completed":len(worker_rows),"total":6,"rows":worker_rows})
            print("P2 AH candidate "+str(len(worker_rows))+"/6",flush=True)
    preparation._write(output/"AH_candidate_summary.json",{"rows":worker_rows,"PH_recomputed":False})
    reference_dir=output/"reference";reference_dir.mkdir()
    native_rows=[];native_summary=[];radial_controls=[];dop_controls=[];gates=[];causes=[];reference_status=[]
    for angle in (22.,39.5):
        items=sorted([j for j in jobs if j["original_job"]["elevation_deg"]==angle],key=lambda j:j["original_job"]["level"])
        ph=[preparation.read_ph(Path(j["PH_raw"]),j["original_job"]) for j in items]
        for t in ph:t["tangent"]=t["p"]/np.linalg.norm(t["p"],axis=1)[:,None]
        baseline=[read_ah(Path(j["AH_meta"]),j["original_job"]) for j in items]
        proposed=[read_ah(candidate/(j["original_job"]["job_id"]+".json"),j["original_job"]) for j in items]
        field=ode.medium(items[2]["original_job"],False)
        max_s=max(prefix(t)["s"][-1] for t in ph)
        reflevels=[RadialReference(ph[2],field,ode,C,10.,level,max_s) for level in (1,2,3)]
        q=grid(.25);references=[r.query(q) for r in reflevels]
        rcontrol=controls(references);rcontrol.update(angle_deg=angle,method="spherical_invariant_quad_root",profiles=QUAD_PROFILES)
        radial_controls.append(rcontrol)
        for level,r in enumerate(references,1):
            np.savez_compressed(reference_dir/f"e{angle:g}_level{level}.npz",s_km=q,xyz_km=r["xyz"],p=r["p"],group_km=r["group"],tangent=r["tangent"])
        fine=reflevels[2]
        dopjob=dict(items[2]["original_job"],scope="local",restart_kind="P2_internal_common")
        dops=[ode.integrate(dopjob,ph[2],field,False,level) for level in (1,2,3)]
        dopvalues=[]
        for solved in dops:
            xyz,p,G=ode.query(solved,q);dopvalues.append({"xyz":xyz,"p":p,"group":G,"g100":float(G[-1])})
        dcontrol=controls(dopvalues);dcontrol.update(angle_deg=angle,method="frozen_P1_DOP853",profiles=ode.PROFILES,not_in_budget=True)
        dop_controls.append(dcontrol)
        cross=errors(references[2],dopvalues[2])
        cross_pass=cross["max_m"]<=.01 and cross["group_rel"]<=1e-8
        ref_invariants=invariant_rows({"xyz":references[2]["xyz"],"p":references[2]["p"]},ode,field,C,10.)
        established=rcontrol["passed"] and cross_pass and ref_invariants["energy_max_drift"]<=1e-10 and ref_invariants["angular_momentum_vector_max_relative_drift"]<=1e-10 and not fine.warnings
        reference_status.append({"angle_deg":angle,"established":established,"quadrature_control":rcontrol,"DOP_control":dcontrol,
            "cross_check":cross,"cross_check_limits":{"max_m":.01,"group_rel":1e-8},"invariants":ref_invariants,"quad_warnings":fine.warnings})
        proof=field_bounds(field,C/1e14,fine.E,max_s)
        trace_sets=[("ph",0.,ph),("ah",.5,baseline),("ah",.125,proposed)]
        for side,fd,traces in trace_sets:
            for item,trace in zip(items,traces):
                if trace.get("status","completed")!="completed":continue
                native=prefix(trace);sample=native["s"]
                truth=fine.query(sample)
                position_error=np.linalg.norm(native["xyz"]-truth["xyz"],axis=1)*1000.
                p_error=np.linalg.norm(native["p"]-truth["p"],axis=1)
                direction_error=np.linalg.norm(native["tangent"]-truth["tangent"],axis=1)
                group_error=native["group"]-truth["group"] if side=="ph" else None
                for i,s in enumerate(sample):
                    native_rows.append({"job_id":item["original_job"]["job_id"],"side":side,"fd_km":fd,"angle_deg":angle,
                        "level":item["original_job"]["level"],"native_index":int(native["native_indices"][i]),"s_km":float(s),
                        "position_error_m":float(position_error[i]),"p_vector_error":float(p_error[i]),"direction_error":float(direction_error[i]),
                        "native_group_error_km":float(group_error[i]) if group_error is not None else None,"within_100km":bool(s<=100.)})
                use=sample<=100.
                row={"job_id":item["original_job"]["job_id"],"side":side,"fd_km":fd,"angle_deg":angle,"level":item["original_job"]["level"],
                    "native_count_in_100":int(np.sum(use)),"native_position_rms_m":float(np.sqrt(np.mean(position_error[use]**2))),
                    "native_position_max_m":float(np.max(position_error[use])),"native_p_error_max":float(np.max(p_error[use])),
                    "native_direction_error_max":float(np.max(direction_error[use])),"invariants":invariant_rows(native,ode,field,C,10.)}
                if side=="ph":
                    intervals=np.diff(native["s"]);slope=np.diff(native["group"])/intervals
                    group_interval_ref=np.diff(truth["group"])/intervals
                    row.update(native_group_error_max_km=float(np.max(abs(group_error[use]))),
                        ng_minus_inverse_p_max=float(np.max(abs(trace["ng_out"][native["native_indices"]]-1/np.linalg.norm(native["p"],axis=1)))),
                        native_group_interval_minus_reference_max=float(np.max(abs(slope-group_interval_ref))),
                        interval_slope_minus_endpoint_average_inverse_p_max=float(np.max(abs(slope-.5*(1/np.linalg.norm(native["p"][:-1],axis=1)+1/np.linalg.norm(native["p"][1:],axis=1))))))
                else:row.update(native_group100_error_km=float(trace["native_group100"]-fine.query(np.array([100.]))["g100"]))
                native_summary.append(row)
        for fd,ah,role in ((.5,baseline,"original_fd0.5"),(.125,proposed,"candidate_fd0.125")):
            for method in METHODS:
                group_id=f"identified|e{angle:g}|fd{fd:g}|{method}"
                levels=[{"level":j["original_job"]["level"],"PH_raw":j["PH_raw"],
                    "AH_meta":j["AH_meta"] if fd==.5 else str(candidate/(j["original_job"]["job_id"]+".json")),
                    "AH_raw":j["AH_raw"] if fd==.5 else str(candidate/(j["original_job"]["job_id"]+".npz")),
                    "status":t.get("status","completed")} for j,t in zip(items,ah)]
                if any(t.get("status")!="completed" for t in ah):
                    gates.append({"group_id":group_id,"angle_deg":angle,"fd_km":fd,"role":role,"method":method,"arc_km":100.,
                        "levels":levels,"eligible":False,"strict_passed":False,"status":"not_evaluable","reasons":["AH candidate execution/domain rejection retained"]})
                    continue
                phv=[reconstruct(prefix(t),method,q) for t in ph];ahv=[reconstruct(prefix(t),method,q) for t in ah]
                actual=errors(phv[2],ahv[2]);integration={"ph":controls(phv),"ah":controls(ahv)}
                sampling={"ph":controls([reconstruct(prefix(ph[2],s),method,q) for s in STRIDES]),
                          "ah":controls([reconstruct(prefix(ah[2],s),method,q) for s in STRIDES])}
                rb={"ph":reconstruction_bound(ph[2],method,fine,proof,"ph"),"ah":reconstruction_bound(ah[2],method,fine,proof,"ah")}
                for r in rb.values():r["bound_established"]=bool(established)
                eligible=established and all(c["passed"] for c in [*integration.values(),*sampling.values()])
                resampling={k:rb["ph"]["bounds"][k]+rb["ah"]["bounds"][k] for k in KEYS}
                d12={side:integration[side]["delta12"] for side in ("ph","ah")};d23={side:integration[side]["delta23"] for side in ("ph","ah")}
                decision=common.numeric_decision_v2(actual,d12,d23,resampling,eligible=eligible)
                reasons=[]
                if not established:reasons.append("independent_reference_not_established")
                for kind,cs in (("integration",integration),("sampling",sampling)):
                    for side,c in cs.items():
                        if not c["passed"]:reasons.append(side+"_"+kind+"_strict_control_failed")
                gate={"group_id":group_id,"angle_deg":angle,"fd_km":fd,"role":role,"method":method,"arc_km":100.,"levels":levels,
                    "status":"diagnostic_complete","metrics":actual,"integration":integration,"sampling":sampling,
                    "query_sampling":[{"step_km":s,"metrics":errors(reconstruct(prefix(ph[2]),method,grid(s)),reconstruct(prefix(ah[2]),method,grid(s)))} for s in QUERY_STEPS],
                    "reconstruction":rb,"resampling":resampling,"budget":{k:d23["ph"][k]+d23["ah"][k]+resampling[k] for k in KEYS},
                    "eligible":bool(eligible),"strict_passed":bool(decision["passed"]),"decision":decision,"reasons":reasons,
                    "method_difference_in_budget":False,"FD_difference_in_budget":False,"AH_group_kind":"native_total_at_exact_100km",
                    "PH_group_kind":"native_G_nodes_reconstructed_at_100_per_method","node_error_in_reconstruction_budget":False}
                gates.append(gate)
                causes.append({"group_id":group_id,**actual,"PH_vs_reference":errors(phv[2],fine.query(q)),"AH_vs_reference":errors(ahv[2],fine.query(q)),
                    "PH_interpolation_only":{k:rb["ph"]["proof"][k] for k in ["actual_reference_interpolation_rms_m","actual_reference_interpolation_max_m","propagated_node_error_query_rms_m"]},
                    "AH_interpolation_only":{k:rb["ah"]["proof"][k] for k in ["actual_reference_interpolation_rms_m","actual_reference_interpolation_max_m","propagated_node_error_query_rms_m"]},
                    "strict_passed":gate["strict_passed"],"FD_bias_is_separate":True})
        all_s=np.array(sorted(fine.cache))
        all_ref=fine.query(all_s)
        np.savez_compressed(reference_dir/f"native_queries_e{angle:g}.npz",s_km=all_s,xyz_km=all_ref["xyz"],p=all_ref["p"],group_km=all_ref["group"],tangent=all_ref["tangent"])
        preparation._write(output/"checkpoint.json",{"angles_complete":angle,"strict_groups":gates,"reference_status":reference_status})
        print("P2 independent angle complete "+str(angle),flush=True)
    write_csv(output/"native_reference.csv",native_rows)
    preparation._write(output/"native_semantics.json",{"cases":native_summary,"C":C,"PH_native_s_not_resampled":True,
        "AH_pointwise_group_available":False,"ng_endpoint_average_is_only_interval_diagnostic":True})
    preparation._write(output/"radial_ODE_controls.json",{"quadrature":radial_controls,"DOP853":dop_controls,"reference_status":reference_status})
    preparation._write(output/"strict_gates.json",{"schema":"a4-c1m-p2-strict-gates-v1","limits":preparation.LIMITS,
        "sample_steps_km":QUERY_STEPS,"node_sampling_strides":STRIDES,"groups":gates,
        "original_baseline_never_replaced_by_candidate":True,"roundoff_waiver_applied":False,
        "method_selection":"Both methods are reported with their own controls and bounds; neither replaces a failed group",
        "diagnostic_computation_complete":True,"physics_passed":False,"global_passed":False})
    preparation._write(output/"cause_matrix.json",{"groups":causes,"node_propagation_and_interpolation_separated":True,"C":C})
    plot_causes(output,causes)
    completion=run/"exit/completion_index_v3.json"
    preparation._write(output/"checks.json",{"diagnostic_complete":True,"source_verified":True,"prepared_identity_verified":True,
        "exit_completed_reference":{"path":str(completion),"sha256":common.sha(completion)} if completion.exists() else None,
        "C":C,"reference_established_by_angle":{str(r["angle_deg"]):r["established"] for r in reference_status},
        "strict_pass_count":sum(g.get("strict_passed",False) for g in gates),"strict_group_count":len(gates),
        "physics_passed":False,"global_passed":False})
    preparation._write(output/"sha256_manifest.json",{"records":[{"path":str(p),"sha256":common.sha(p)} for p in sorted(output.rglob("*")) if p.is_file()]})
    print("P2 independent analysis complete; strict="+str(sum(g.get("strict_passed",False) for g in gates))+"/"+str(len(gates)),flush=True)


def plot_causes(output: Path,rows: list[dict]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.2),layout="constrained")
    for ax,side in zip(axes,("PH_vs_reference","AH_vs_reference")):
        for method,marker in zip(METHODS,("o","s")):
            for fd,style in ((.5,"-"),(.125,"--")):
                selected=[r for r in rows if f"fd{fd:g}|{method}" in r["group_id"]]
                if not selected:continue
                angles=[float(r["group_id"].split("|")[1][1:]) for r in selected]
                ax.plot(angles,[r[side]["rms_m"] for r in selected],marker=marker,linestyle=style,label=f"{method}, FD={fd:g}")
        ax.set(xlabel="Launch elevation (deg)",ylabel="Fixed 100-km RMS relative to radial reference (m)",title=side.replace("_"," "))
        ax.grid(alpha=.25);ax.legend(fontsize=7)
    fig.savefig(output/"independent_residual_sources.png",dpi=180)
    fig.savefig(output/"independent_residual_sources.pdf")
    plt.close(fig)
