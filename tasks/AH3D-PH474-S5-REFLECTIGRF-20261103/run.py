"""New-task bounded engineering/reference/AH diagnostics; no legacy entry points."""
from pathlib import Path
import csv,hashlib,importlib.util,json,os,subprocess,sys,time,warnings
from datetime import datetime,timezone
import numpy as np
from scipy.io import savemat,loadmat
from scipy.interpolate import PchipInterpolator,interp1d
import reference as ref
HERE=Path(__file__).resolve().parent
TASK_ROOT=HERE
ROOT=HERE.parents[1]
TASK=HERE.name
OLD=ROOT/'tasks/AH3D-PH474-S5-REFLECTOX-20261102/refinement'
PREVIOUS=ROOT/'tasks/AH3D-PH474-S5-TRANSNE-20261102'
MEX=Path('D:/pharlap/pharlap_4.7.4/mex/raytrace_3d_sp.mexw64')
MEX_SHA='5D052459354ED105B967274EE9859E4EF0B6AAE2DA14E7B5353049105FC3AF2C'
LAUNCH=dict(lat_deg=30.,lon_deg=120.,height_km=0.,elevation_deg=25.,azimuth_deg=0.)
PROFILES=dict(L1=dict(rtol=1e-9,atol_r_km=1e-7,atol_p=1e-11,atol_t=1e-12,ds_min=.005,ds_max=1.),L2=dict(rtol=1e-10,atol_r_km=1e-8,atol_p=1e-12,atol_t=1e-12,ds_min=.002,ds_max=.5))
PH_TOL=dict(L1=[1e-9,.001,1.],L2=[1e-10,.001,1.])
FD_BASE_KM=.03125
sys.path.insert(0,str(ROOT))
from models.ne_provider_geocentric import GeocentricGridNeProvider3D
class TransverseProvider(GeocentricGridNeProvider3D):
    def get_Ne_geocentric(self,lat_deg,lon_deg,h_km,t=0.):
        return super().get_Ne_geocentric(lat_deg,lon_deg,h_km,t)*(1.+.006*(lon_deg-120.))
    def get_metadata(self):
        result=super().get_metadata();result['transverse_factor']='1+0.006*(lon_deg-120)';return result
def make_provider():
    ne=np.broadcast_to(ref.NE,(len(ref.LAT),len(ref.LON),len(ref.H))).copy()
    provider=TransverseProvider(ne,ref.LAT,ref.LON,ref.H,ref.R,interpolation='radial_pchip')
    assert np.array_equal(provider._radial_spline.c,ref.COEFFICIENTS)
    return provider
def remaining():
    clock=load(HERE/'execution-clock.json')
    elapsed=(datetime.now(timezone.utc)-datetime.fromisoformat(clock['execution_started_at'])).total_seconds()
    excluded=sum((datetime.fromisoformat(x['end'])-datetime.fromisoformat(x['start'])).total_seconds() for x in load(HERE/clock['independent_review_intervals_file']))
    return clock['execution_budget_seconds']-elapsed+excluded
def identity():
    baseline=load(HERE/'identity-before.json')
    assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()==baseline['head']
    assert all(sha(ROOT/x['path'])==x['sha256'] for x in baseline['protected'])
    assert sha(PREVIOUS/'science-manifest.json')=='B8405BB178C848200649BC5A528A1ECB9C0265AFCF3766FEDC4635FB25705711'
    assert all(sha(ROOT/x['path'])==x['sha256'] for x in load(PREVIOUS/'science-manifest.json')['files'])
    gate=load(HERE/'plan-acceptance.json');assert gate['decision']=='PLAN_ACCEPT'
    assert gate['plan_sha256']==sha(HERE.parent/(TASK+'-plan.md'))=='C75B99D3E3EF7A072EDFCCE02DD1F83F9C9A3FE03341A215AE51C2F85DB48918'
    assert gate['review_sha256']==sha(ROOT/gate['review_file'])
def actual_gate(endpoint,mode):
    gate=load(HERE/'engineering-accept.json');assert gate['decision']=='ACCEPT'
    assert gate['engineering_manifest_sha256']==sha(HERE/'engineering-manifest.json')
    assert endpoint+'_'+mode in gate['accepted_dependencies']
    assert all(sha(ROOT/x['path'])==x['sha256'] for x in load(HERE/'engineering-manifest.json')['files'])


LAUNCH = dict(tx_lat_deg=30., tx_lon_deg=120., tx_height_km=0., elevation_deg=25., bearing_deg=0.)
PROFILES = {
    'L1': dict(rtol=1e-9, atol_r_km=1e-7, atol_p=1e-11, ds_min_km=.005, ds_max_km=1.),
    'L2': dict(rtol=1e-10, atol_r_km=1e-8, atol_p=1e-12, ds_min_km=.002, ds_max_km=.5)}
PH_TOL = {'L1': [1e-9, .001, 1.], 'L2': [1e-10, .001, 1.]}
FD_BASE_KM = .03125
STATE_FIELDS = ['pos_x', 'pos_y', 'pos_z', 'dir_x', 'dir_y', 'dir_z', 'group_path',
                'geometrical_path', 'phase_path', 'absorption', 'indep_var', 'ODE_step_size']

def model_for(mode, level='L1', provider=None):
    sys.path.insert(0, str(ROOT))
    import config
    from models.ray_tracer_ah_3d_ecef import RefractiveIndexAH3DECEF
    from igrf_field import production_provider
    if provider is None:
        provider=make_provider()
    profile = TASK+'_'+level
    settings = {k: v for k, v in PROFILES[level].items() if k != 'atol_p'}
    settings['atol_k'] = ref.K0*PROFILES[level]['atol_p']
    config.UPDATE2['AH3D_INTEGRATION'][profile] = settings
    b_provider = production_provider()
    launch_b = b_provider.get_B_ecef_T(ref.XYZ0)
    model = RefractiveIndexAH3DECEF(provider, earth_radius_km=ref.R, freq_MHz=12., wave_mode=mode,
        geomagnetic_provider=b_provider, gyro_MHz=ref.CONSTANTS['AH'][1]*np.linalg.norm(launch_b),
        spatial_fd_base_km=FD_BASE_KM, integration_profile=profile, plasma_coefficient=ref.CONSTANTS['AH'][0],
        compatibility_mode=None)
    assert model.k0 == ref.K0 and model.compatibility_mode is None
    return model, provider, profile

def utc():
    return datetime.now(timezone.utc).isoformat()

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest().upper()

def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, np.ndarray):
        return clean(value.tolist())
    if isinstance(value, np.generic):
        return clean(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value

def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(clean(value), indent=2, allow_nan=False), encoding='utf-8')

def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))

def event(kind, **kwargs):
    with (HERE/'execution-events.jsonl').open('a', encoding='utf-8') as stream:
        stream.write(json.dumps(clean(dict(at=utc(), kind=kind, **kwargs)))+'\n')

def error(actual, expected):
    return float(np.linalg.norm(np.asarray(actual)-expected)/max(np.linalg.norm(expected), 1e-30))

def configurations_pass(folders):
    return bool(folders) and all((p/'audit.json').exists() and load(p/'audit.json')['status'] == 'PASS' for p in folders)

def engineering():
    import engineering_checks
    engineering_checks.run()

def reference_one(endpoint, mode, level):
    identity()
    assert load(HERE/'engineering.json')['status'] == 'PASS'
    directory = HERE/'references'/f'{endpoint}_{mode}_{level}'
    directory.mkdir(parents=True, exist_ok=False)
    event('REFERENCE_START', endpoint=endpoint, mode=mode, level=level)
    started = time.monotonic()
    try:
        with warnings.catch_warnings(record=True) as recorded_warnings:
            warnings.simplefilter('always')
            sol, grid, curve, metadata = ref.integrate(endpoint, mode, level)
        metadata['scipy_warnings'] = [str(item.message) for item in recorded_warnings]
        metadata['actual_rtol'] = max(metadata['rtol'], 100*np.finfo(float).eps)
        metadata['actual_atol'] = np.full(7, metadata['atol'])
        # Save the complete returned raw before applying any numerical acceptance gate.
        np.savez_compressed(directory/'raw.npz', s=sol.t, state=sol.y.T, grid_s=grid, grid_state=curve)
        save(directory/'metadata.json', metadata)
        shells = []
        for xyz, p in zip(sol.y[:3].T, sol.y[3:6].T):
            ne = ref.density(xyz, endpoint)[0]
            b = ref.field(xyz, endpoint)[0]
            shells.append(abs(float(p@p)-ref.schur(ne, p, mode, endpoint, b=b)['n2']))
        q = np.array([ref.llh(row[:3]) for row in sol.y.T])
        radial = np.einsum('ij,ij->i', sol.y[:3].T, sol.y[3:6].T)/np.linalg.norm(sol.y[:3].T, axis=1)
        checks = dict(solver=metadata['success'], finite=bool(np.all(np.isfinite(sol.y))),
            shell=max(shells) <= 1e-10, group_monotone=bool(np.all(np.diff(sol.y[6]) > 0.)),
            nonnegative_ne=bool(all(ref.density(x, endpoint)[0] >= 0. for x in sol.y[:3].T)),
            support=bool(np.all(q[:, 2] >= -1e-8) and np.all(q[:, 2] <= 460.+1e-8)),
            natural_event=metadata['event'] in ('GROUND', 'TOP_OUT'),
            ground_turn=bool(metadata['event'] != 'GROUND' or (np.any(radial > 0.) and radial[-1] < 0. and np.max(q[:, 2]) > 1.)))
        audit = dict(status='PASS' if all(checks.values()) else 'FAIL', checks=checks,
            saved_nodes=len(sol.t), max_shell=max(shells), final_s_km=sol.t[-1], final_group_km=sol.y[6, -1],
            max_height_km=float(max(q[:, 2])), final_radial_p=float(radial[-1]), event=metadata['event'])
        save(directory/'audit.json', audit)
        event('REFERENCE_END', endpoint=endpoint, mode=mode, level=level, status=audit['status'], elapsed_s=time.monotonic()-started)
        print(json.dumps(clean(audit)))
    except Exception as exc:
        save(directory/'failure.json', dict(error=type(exc).__name__, message=str(exc), elapsed_s=time.monotonic()-started,
            raw_exists=(directory/'raw.npz').exists(), automatic_retry=False))
        event('REFERENCE_FAILED', endpoint=endpoint, mode=mode, level=level, error=type(exc).__name__)
        raise

def child(command, timeout):
    identity()
    duration = min(float(timeout), remaining()-5.)
    assert duration > 0.
    name = '_'.join(command)
    path = HERE/(name+'.log')
    assert not path.exists(), 'No automatic replay'
    started = time.monotonic()
    env = dict(os.environ)
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NRT_NUM_THREADS'):
        env[key] = '1'
    with path.open('xb') as output:
        try:
            completed = subprocess.run([sys.executable, '-u', str(HERE/'run.py'), *command], cwd=ROOT, env=env,
                stdout=output, stderr=subprocess.STDOUT, timeout=duration,
                creationflags=subprocess.CREATE_NO_WINDOW if os.name == 'nt' else 0)
            code = completed.returncode
        except subprocess.TimeoutExpired:
            code = 'TIMEOUT'
    receipt = dict(command=command, elapsed_s=time.monotonic()-started, timeout_seconds=duration,
        exit_code=code, log_sha256=sha(path), automatic_retry=False)
    save(HERE/(name+'.receipt.json'), receipt)
    return receipt

def ah_one(mode, level):
    identity()
    actual_gate('AH', mode)
    sys.path.insert(0, str(ROOT))
    from models import ray_tracer_ah_3d_ecef as ah
    directory = HERE/'AH'/f'{mode}_{level}'
    directory.mkdir(parents=True, exist_ok=False)
    model, provider, profile = model_for(mode, level)
    queries = []
    magnetic_queries = []
    current_s = [0.]
    start = time.monotonic()
    original_get = provider.get_Ne_geocentric
    original_geodetic = model.geomagnetic_provider._geodetic
    original_rhs = ah.ah3d_ecef_rhs

    def observed_get(lat, lon, height, t=0.):
        if time.monotonic()-start > 890.:
            raise TimeoutError('AH internal 890 s budget')
        queries.append([current_s[0], lat, lon, height])
        return original_get(lat, lon, height, t)

    def observed_rhs(s, state, *args, **kwargs):
        current_s[0] = s
        return original_rhs(s, state, *args, **kwargs)

    def observed_geodetic(xyz):
        q = original_geodetic(xyz)
        magnetic_queries.append([current_s[0], *q])
        return q

    provider.get_Ne_geocentric = observed_get
    model.geomagnetic_provider._geodetic = observed_geodetic
    ah.ah3d_ecef_rhs = observed_rhs
    save(directory/'request.json', dict(mode=mode, level=level, profile=profile,
        settings=PROFILES[level], actual_atol_k=ref.K0*PROFILES[level]['atol_p'], constants=ref.CONSTANTS['AH'],
        launch=LAUNCH, B_ecef_T=ref.field(ref.XYZ0, 'AH')[0], spatial_fd_base_km=FD_BASE_KM, initial=ah.launch_state_ah3d_ecef(model, **LAUNCH)))
    event('AH_START', mode=mode, level=level)
    try:
        ray = ah.integrate_fixed_path_ah3d_ecef(model, **LAUNCH, path_length_km=5000., profile=profile,
            h_stop_km=460., stop_on_physical_events=True, apply_dispersion_projection=True)
        arrays = {k: v for k, v in ray.items() if isinstance(v, np.ndarray)}
        np.savez_compressed(directory/'raw.npz', **arrays)
        save(directory/'metadata.json', {k: v for k, v in ray.items() if not isinstance(v, np.ndarray)})
        state = np.column_stack([ray['trajectory_ecef_km'], np.asarray(ray['wavevector_ecef_km^-1'])/ref.K0,
                                 ray['trajectory_group_path_km']])
        np.savez_compressed(directory/'normalized.npz', s=ray['trajectory_geom_path_km'], state=state)
        save(directory/'elapsed.json', dict(elapsed_s=time.monotonic()-start, termination=ray['termination_status']))
        event('AH_END', mode=mode, level=level, termination=ray['termination_status'])
        print('AH_COMPLETE', mode, level, ray['termination_status'])
    except Exception as exc:
        save(directory/'failure.json', dict(error=type(exc).__name__, message=str(exc), elapsed_s=time.monotonic()-start))
        event('AH_FAILED', mode=mode, level=level, error=type(exc).__name__)
        raise
    finally:
        provider.get_Ne_geocentric = original_get
        model.geomagnetic_provider._geodetic = original_geodetic
        ah.ah3d_ecef_rhs = original_rhs
        np.savez_compressed(directory/'provider_queries.npz', s_lat_lon_h=np.array(queries))
        np.savez_compressed(directory/'magnetic_queries.npz', s_wgs_lat_lon_h=np.array(magnetic_queries))

if __name__=='__main__':
    action,*args=sys.argv[1:]
    if action=='grid-build':
        identity()
        import igrf_field
        print(json.dumps(igrf_field.build_grids()), flush=True)
    elif action=='reference-one':reference_one(*args)
    elif action=='ah-one':ah_one(*args)
    else:raise ValueError(action)
