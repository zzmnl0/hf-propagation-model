"""Read-only stage dispatch for the frozen P3.1 run; never silently restart experiments."""
from __future__ import annotations
import argparse, hashlib, inspect, json, os, re, shutil, subprocess, sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PARENT=ROOT/'output2/a4_pharlap/c1/c1m/p3_20261009_v1'
PLAN=ROOT/'word/a4_c1m_p31_implementation_plan_20261009.txt'
class Blocked(RuntimeError): pass
def read(path):
    if not Path(path).is_file():raise Blocked(f'Required artifact missing: {path}')
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def require_hash(path,expected):
    if not Path(path).is_file() or sha(path)!=expected:raise ValueError(f'Identity mismatch: {path}')
def contained(run,value):
    path=(run/value).resolve()
    if not path.is_relative_to(run.resolve()):raise ValueError('Artifact escapes run')
    return path
def write(path,data):Path(path).write_text(json.dumps(data,indent=2,allow_nan=False),encoding='utf-8')
def execute(run,script,*arguments):
    source=ROOT/'tests'/script
    frozen=run/'source_before/tests'/script
    if not frozen.is_file():raise Blocked('Execution source was not frozen; prepare a new run: '+script)
    require_hash(source,sha(frozen))
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    result=subprocess.run([sys.executable,'-B','-X','utf8',str(source),*map(str,arguments)],cwd=ROOT,env=env)
    if result.returncode:raise Blocked(f'{script} exited {result.returncode}; preserve this failure and resume its independent subexperiment.')
def freeze_new(run):
    if run.exists():raise Blocked('Existing directory has no baseline freeze; select a new run rather than overwriting partial preparation.')
    delivery=read(PARENT/'delivery_sha256.json')
    for r in delivery['records']:require_hash(PARENT/r['path'],r['sha256'])
    for r in delivery['external_tracked_sources_and_report']:require_hash(r['path'],r['sha256'])
    if not PLAN.is_file():raise Blocked('Approved plan text missing')
    template=read(ROOT/'output2/a4_pharlap/c1/c1m/p31_20261009_v1/baseline_freeze.json')
    paths={r['path'] for r in template['records']}
    paths.update(str(p.relative_to(ROOT)) for p in (ROOT/'tests').glob('*p31*.py'))
    run.mkdir(parents=True);records=[]
    for name in sorted(paths):
        original=ROOT/name;destination=run/'source_before'/name;destination.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(original,destination);records.append(dict(path=name,sha256=sha(original)))
    shutil.copy2(PLAN,run/'approved_plan.txt')
    write(run/'plan_identity.json',dict(path=str(PLAN),sha256=sha(PLAN),source='Approved P3.1 plan; frozen by stage runner'))
    write(run/'baseline_freeze.json',dict(parent_run=str(PARENT),parent_manifest_sha256=sha(PARENT/'delivery_sha256.json'),parent_verified_files=len(delivery['records'])+len(delivery['external_tracked_sources_and_report']),parent_errors=[],head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),git_status=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True),python=sys.executable,records=records))
def prepare(run):
    if not (run/'baseline_freeze.json').exists():freeze_new(run)
    frozen=read(run/'baseline_freeze.json');parent=Path(frozen['parent_run']);manifest=parent/'delivery_sha256.json'
    require_hash(manifest,frozen['parent_manifest_sha256']);delivery=read(manifest)
    for record in delivery['records']:require_hash(parent/record['path'],record['sha256'])
    for record in delivery['external_tracked_sources_and_report']:require_hash(record['path'],record['sha256'])
    for record in frozen['records']:require_hash(run/'source_before'/record['path'],record['sha256'])
    plan=read(run/'plan_identity.json');require_hash(plan['path'],plan['sha256'])
    return dict(status='verified',baseline_sources=len(frozen['records']),parent_files=len(delivery['records']),production_modified=False)
def boundary(run):
    d=run/'boundary'
    if not (d/'manifest.json').exists():
        prepare(run);execute(run,'_a4_c1m_p31_boundary.py','prepare','--run',run)
    manifest=read(d/'manifest.json')
    for name,h in manifest['source_sha256'].items():require_hash(d/'source'/name,h)
    for name,h in manifest['input_sha256'].items():require_hash(d/f'{name}.mat',h)
    for j in manifest['jobs']:require_hash(j['grid_mat'],j['grid_sha256'])
    for role in ('training','holdout'):
        missing=any(not (d/'pharlap_raw'/f"{j['job_id']}.mat").is_file() for j in manifest['jobs'] if j['role']==role)
        if missing:execute(run,'_a4_c1m_p31_boundary.py','ph','--role',role,'--run',run)
        if not (d/f'{role}_analysis/metrics.json').exists():execute(run,'_a4_c1m_p31_boundary.py','analyze','--role',role,'--run',run)
        for record in read(d/f'{role}_analysis/metrics.json'):require_hash(d/'pharlap_raw'/f"{record['job_id']}.mat",record['raw_sha256'])
    if not (d/'diagnostics/conclusion.json').exists():execute(run,'_a4_c1m_p31_boundary.py','diagnose','--run',run)
    if not (d/'endpoint_supplement_v1/conclusion.json').exists():execute(run,'_a4_c1m_p31_boundary_endpoint.py','--run',run)
    if not (run/'boundary_independent_check/verification.json').exists():execute(run,'_verify_a4_c1m_p31_boundary_independent.py','--run',run)
    for item in ('diagnostics','endpoint_supplement_v1'):
        checks=read(d/item/'sha256.json');base=d if item=='diagnostics' else d/item
        for name,h in checks.items():require_hash(base/name,h)
    conclusion=read(d/'diagnostics/conclusion.json');independent=read(run/'boundary_independent_check/verification.json')
    if not independent['passed']:raise ValueError('Independent boundary verification failed')
    return dict(status='reused_verified',jobs=len(manifest['jobs']),boundary_convention_identified=conclusion['boundary_convention_identified'],formal_20km_endpoint_acceptance=read(d/'endpoint_supplement_v1/conclusion.json')['formal_20km_endpoint_acceptance'])
def reference(run):
    prepare(run)
    if not (run/'reference_certification_index.json').exists():
        frozen=read(run/'baseline_freeze.json');parent=Path(frozen['parent_run'])
        for folder,script in (('reference_r3_v2','_a4_c1m_p31_certify.py'),('reference_r3_v4','_a4_c1m_p31_certify_finalize.py')):
            summary=run/folder/'summary.json'
            if not summary.exists() or read(summary).get('case_count')!=16:execute(run,script,'--parent',parent,'--run',run)
        raise Blocked('R3 certification/finalization executed. Freeze independent verification and reference_certification_index.json before acceptance. Fullpath controls require _a4_c1m_p31_fullpath.py --run RUN and its reviewed version selection; this runner does not invent that index.')
    index=read(run/'reference_certification_index.json');directory=contained(run,index['analysis_directory']);verification=contained(run,index['independent_verification'])
    require_hash(directory/'summary.json',index['analysis_summary_sha256']);require_hash(verification,index['independent_verification_sha256'])
    summary=read(directory/'summary.json');independent=read(verification)
    if not independent.get('all_16_reproduced',False) or independent.get('case_count')!=16:raise ValueError('Independent reference verification failed or did not cover 16 cases')
    for case in summary['cases']:
        result=read(directory/case['case_id']/'result.json')
        for record in result['sources']:require_hash(record['path'],record['sha256'])
        if result['reference_certified_passed']!=case['reference_certified_passed']:raise ValueError('Reference case summary mismatch')
    full=read(run/'reference_fullpath_index.json');controls=read(contained(run,full['selected']));pieces=read(contained(run,full['nonlinear_control']))
    require_hash(contained(run,full['selected']).parent/'source.py',controls['source_sha256'])
    return dict(status='reused_verified',case_count=summary['case_count'],reference_certified_count=summary['reference_certified_count'],all_16_certified=summary['all_16_certified'],fullpath_analytic_controls=bool(controls['passed'] and pieces['passed']),formal_smooth_fullpath_verified=controls['formal_smooth_fullpath_verified'])
def repair(run):
    result=boundary(run)
    if not result['boundary_convention_identified']:raise Blocked('Boundary propagation convention unresolved: no production patch, formal smooth-background regression or new holdout is authorized by evidence.')
    raise Blocked('A confirmed mechanism requires the concrete reviewed production patch; this dispatcher does not implement a repair automatically.')
def verify(run):
    # Verification must never dispatch a missing scientific experiment.
    read(run/'baseline_freeze.json');read(run/'reference_certification_index.json')
    manifest=read(run/'boundary/manifest.json')
    for j in manifest['jobs']:
        if not (run/'boundary/pharlap_raw'/f"{j['job_id']}.mat").is_file():raise Blocked('Run boundary stage to complete missing raw data before verification')
    for name in ('training_analysis/metrics.json','holdout_analysis/metrics.json','diagnostics/conclusion.json','endpoint_supplement_v1/conclusion.json'):read(run/'boundary'/name)
    read(run/'boundary_independent_check/verification.json')
    result=dict(baseline=prepare(run),boundary=boundary(run),reference=reference(run))
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    check=subprocess.run([sys.executable,'-B','-X','utf8',str(ROOT/'tests/_verify_a4_c1m_p31_reference.py')],cwd=ROOT,env=env,capture_output=True,text=True)
    if check.returncode:raise ValueError('Precision verifier failed: '+check.stdout+check.stderr)
    result.update(status='checks_complete',precision_verifier_stdout=check.stdout,production_modified=False,recommend_p4=False)
    return result

def physics_exposure_audit(run):
    """Distinguish a recorded experiment from a proposed angle in prose."""
    candidates=((11.,27.,33.),(11.5,26.,32.))
    matches=[];checked=0
    def scan(value,path):
        if isinstance(value,dict):
            frequency=value.get('frequency_MHz',value.get('freq_MHz'))
            angle=value.get('elevation_deg',value.get('elev_deg'))
            for index,(f,a,b) in enumerate(candidates):
                if frequency==f and angle in (a,b):
                    proposal=path.name.endswith('exposure_audit.json') or path.name=='case_registry.json'
                    conditional=path.name=='holdout_conditional_topology.json'
                    classification='proposal_only' if proposal else ('conditional_reference_prediction_only' if conditional else 'historical_condition')
                    matches.append(dict(candidate=index,path=str(path),sha256=sha(path),frequency_MHz=f,elevation_deg=angle,
                                        classification=classification,counts_as_exposed=not(proposal or conditional),
                                        evidence='No AH execution or parameter selection in proposal/reference preregistration' if proposal or conditional else 'structured historical condition; conservatively counted as exposed'))
            for child in value.values():scan(child,path)
        elif isinstance(value,list):
            for child in value:scan(child,path)
    for path in sorted((ROOT/'output2/a4_pharlap').rglob('*.json')):
        if run in path.parents:continue
        checked+=1
        try:scan(read(path),path)
        except (ValueError,UnicodeError):continue
    text_records=[]
    texts=list((ROOT/'tests').glob('*.py'))+list((ROOT/'word').glob('*.txt'))+list((ROOT/'word').glob('*.md'))
    texts += [ROOT.parent/'plan_next_steps.md']
    for path in texts:
        content=path.read_text(encoding='utf-8-sig',errors='replace')
        for index,(f,a,b) in enumerate(candidates):
            for line_no,line in enumerate(content.splitlines(),1):
                # Script literals defining the audit itself are not experiments.
                if re.search(r'(?<![\d.])'+re.escape(str(f).rstrip('0').rstrip('.'))+r'(?![\d.])',line) and all(re.search(r'(?<!\d)'+str(int(x))+r'(?!\d)',line) for x in (a,b)):
                    text_records.append(dict(candidate=index,path=str(path),line=line_no,text=line[:500],
                                             classification='proposal_or_script_text_requires_context'))
    exposed=[any(row['candidate']==i and row['counts_as_exposed'] for row in matches) for i in range(2)]
    selected=next((i for i,value in enumerate(exposed) if not value),None)
    topology_path=ROOT/'output2/a4_pharlap/c1/c1m/p31_20261009_v1/holdout_conditional_topology.json'
    topology_context=None
    if topology_path.is_file():
        topology_context=dict(path=str(topology_path),sha256=sha(topology_path),record=read(topology_path),
                              old_assumption='strict supplied polynomial rejected negative density; clip_zero predictions were conditional and not valid formal registrations',
                              new_medium='nonnegative radial_pchip, strict density policy, current config.K_FP**2; old topology values will not be inherited')
    return dict(json_files_checked=checked,candidates=[dict(frequency_MHz=f,elevations_deg=[a,b],exposed=exposed[i]) for i,(f,a,b) in enumerate(candidates)],
                actual_or_conservative_structured_matches=matches,text_matches=text_records,selected_candidate=selected,
                historical_reference_context=topology_context,
                method='Structured historical conditions count as exposure; requests, prepared case registries and conditional reference preregistrations do not. Prose/script matches are saved for context audit.',
                source_scope='all JSON under output2/a4_pharlap; current test scripts, word documents and overall plan')

def physics_prepare(run):
    """Freeze new physical inputs without starting AH, PH or reference traces."""
    sys.path.insert(0,str(ROOT))
    import numpy as np
    import scipy
    import mpmath
    from scipy.io import loadmat
    from scipy.interpolate import PchipInterpolator
    import config
    import _run_a4_c1m_attribution as common
    from models.ray_tracer_ah_3d_ecef import (
        RefractiveIndexAH3DECEF, integrate_fixed_path_ah3d_ecef, trace_ray_ah3d_ecef)
    plan=ROOT/'word/ah3d_step1_p31_physical_acceptance.md'
    if run.parent!=ROOT/'output2/a4_pharlap/c1/c1m' or not run.name.startswith('p31_physics_'):
        raise ValueError('Physical preparation requires output2/a4_pharlap/c1/c1m/p31_physics_<run_id>')
    if run.exists():raise Blocked('Physical preparation never overwrites an existing run; use a new run identity')
    if not plan.is_file():raise Blocked('Physical plan missing')
    exposure=physics_exposure_audit(run)
    run.mkdir(parents=True)
    write(run/'exposure_audit.json',exposure)
    if exposure['selected_candidate'] is None:
        write(run/'preparation_failure.json',dict(reason='Both preregistered holdout sets exposed',recommend_p4_physics=False))
        raise Blocked('Both preregistered holdout sets exposed; preserve audit and register a new set')
    parent_manifest=PARENT/'experiments/manifest.json'
    original=read(parent_manifest)
    jobs=[j for j in original['jobs'] if j['role'] in ('R1','R2','R3','R4') and j['level']==1]
    if len(jobs)!=24:raise ValueError('Historical formal registry must contain 24 distinct physical cases before new holdouts')
    counts={role:sum(j['role']==role for j in jobs) for role in ('R1','R2','R3','R4')}
    if counts!={'R1':2,'R2':4,'R3':16,'R4':2}:raise ValueError('Historical role counts differ from approved physical matrix')
    if any(j['origin_height_km']!=100. or j['scope']!='first_100_km' for j in jobs if j['role'] in ('R1','R3')):
        raise ValueError('Historical local origin/scope conflicts with physical plan')
    profiles=[dict(value) for value in common.AH_PROFILES]
    if len(profiles)!=5:raise ValueError('Expected three historical integration levels and two legal candidate refinements')
    originals={};media={}
    def save_medium(name,height,density,lat,lon,source,step,kind,gradient=None):
        pp=PchipInterpolator(height,density,extrapolate=False)
        path=run/'inputs'/name/f'medium_d{step:g}'.replace('.','p')
        path=path.with_suffix('.npz');path.parent.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(path,h=height,height_km=height,coefficients=pp.c,Ne1d=density,Ne=density,
                            latitude_deg=lat,longitude_deg=lon,lat=lat,lon=lon)
        row=dict(path=str(path.relative_to(run)),sha256=sha(path),source=source,step_km=step,
                 height_domain_km=[float(height[0]),float(height[-1])],interpolation='radial_pchip',
                 continuous_background_kind=kind,boundary_semantics='reject_nonzero_density_or_horizontal_crossing; zero_density_edge_to_vacuum',
                 minimum_node_density_m3=float(np.min(density)))
        if gradient is not None:
            row['affine_gradient_m3_per_km']=gradient
            row['exact_affine_coefficients']=bool(np.all(pp.c[:2]==0) and np.all(pp.c[2]==gradient))
            row['coefficient_curvature_max']=float(np.max(np.abs(pp.c[:2])))
        return row
    for j in jobs:
        name=j['case_id']
        if name in media:continue
        source=Path(j['grid_mat']);require_hash(source,j['grid_sha256'])
        grid=loadmat(source,simplify_cells=True)
        h=np.asarray(grid['height_km'],float);ne=np.asarray(grid['Ne_m3'],float)
        one=ne[0,0]
        if not np.array_equal(ne,np.broadcast_to(one,ne.shape)):raise ValueError('Formal input is not radial')
        if any(np.any(grid[key]!=0) for key in ('Bx_T','By_T','Bz_T','collision_freq_Hz')):raise ValueError('Formal input is not B0 collisionless')
        if float(grid['earth_radius_m'])!=6371000.:raise ValueError('Formal input sphere mismatch')
        destination=run/'inputs'/name/'original.mat';destination.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,destination)
        originals[name]=dict(path=str(destination.relative_to(run)),sha256=sha(destination),original_path=str(source),original_sha256=j['grid_sha256'])
        kind='affine' if j['role']=='R3' else ('internal' if j['role']=='R1' else 'smooth')
        gradient=float((one[-1]-one[0])/(h[-1]-h[0])) if kind=='affine' else None
        expected=common._a2_profile(kind,h,gradient=gradient or 0.)
        if not np.array_equal(one,expected):raise ValueError('Continuous generating background does not reproduce frozen original nodes: '+name)
        spacing=float(h[1]-h[0]);steps=(2.,1.,.5,.25) if kind=='affine' else (2.,1.,.5)
        if kind!='affine' and (h[0]!=0 or h[-1]!=460 or spacing!=2):raise ValueError('Complete background baseline/domain changed')
        variants=[]
        for step in steps:
            heights=h.copy() if step==spacing else np.linspace(h[0],h[-1],int(round((h[-1]-h[0])/step))+1)
            values=one.copy() if step==spacing else common._a2_profile(kind,heights,gradient=gradient or 0.)
            variants.append(save_medium(name,heights,values,grid['latitude_deg'],grid['longitude_deg'],originals[name],step,kind,gradient))
        baseline=next(row for row in variants if row['step_km']==spacing)
        media[name]=dict(baseline=baseline,grid_variants=variants,continuous_background=dict(kind=kind,gradient_m3_per_km=gradient,
                       source='tests/_run_a4_c1m_attribution.py:_a2_profile',chapman_peak_m3=607635874822.0494,
                       chapman_peak_height_km=350.,chapman_scale_height_km=100.,smooth_entrance_km=[20.,80.],nodes_exactly_reproduced=True))
    cases=[]
    for j in jobs:
        case=dict(case_id=j['job_id'].rsplit('_t',1)[0],role=j['role'],frequency_MHz=j['frequency_MHz'],
                  elevation_deg=j['elevation_deg'],origin_height_km=j['origin_height_km'],scope=j['scope'],
                  tx_lat_deg=30.,tx_lon_deg=120.,bearing_deg=0.,medium_path=media[j['case_id']]['baseline']['path'],
                  medium_id=j['case_id'],source_job=j,source_manifest=str(parent_manifest),source_manifest_sha256=sha(parent_manifest))
        cases.append(case)
    chosen=exposure['candidates'][exposure['selected_candidate']]
    for angle in chosen['elevations_deg']:
        cases.append(dict(case_id=f"holdout_f{chosen['frequency_MHz']:g}_e{angle:g}".replace('.','p'),role='holdout',
                          frequency_MHz=chosen['frequency_MHz'],elevation_deg=angle,origin_height_km=0.,scope='natural_full_path',
                          tx_lat_deg=30.,tx_lon_deg=120.,bearing_deg=0.,medium_id='smooth_s5',medium_path=media['smooth_s5']['baseline']['path'],
                          source='preregistered physical plan holdout; exposure_audit.json'))
    for case in cases:
        case['frequency_Hz']=case['frequency_MHz']*1e6
        case['max_path_km']=5000.
        case['h_stop_km']=460.
        case['grid_variants']=media[case['medium_id']]['grid_variants']
        case['reference_topology']=dict(expected='local_upward' if case['scope']=='first_100_km' else 'reflected_single_hop',
                                        registration_status='awaiting_independent_reference',source='approved case scope; not AH results')
    if len({case['case_id'] for case in cases})!=26:raise ValueError('Physical registry is not 26 unique conditions')
    source_paths={plan,ROOT.parent/'plan_next_steps.md'}
    source_paths.update(ROOT.glob('*.py'))
    for folder in ('models','tests','analysis','viz'):
        source_paths.update(p for p in (ROOT/folder).rglob('*.py') if 'legacy' not in p.parts)
    source_records=[]
    for path in sorted(source_paths):
        name=path.relative_to(ROOT) if path.is_relative_to(ROOT) else Path('overall_plan_next_steps.md')
        destination=run/'source_before'/name;destination.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,destination)
        source_records.append(dict(path=str(path),snapshot_path=str(destination.relative_to(run)),sha256=sha(path)))
    for name,command in (('git_head.txt',['git','rev-parse','HEAD']),('git_status.txt',['git','status','--short']),
                         ('git_diff.patch',['git','diff','--binary','HEAD']),('git_diff_stat.txt',['git','diff','--stat','HEAD'])):
        result=subprocess.run(command,cwd=ROOT,capture_output=True)
        if result.returncode:raise ValueError('Git identity command failed: '+repr(command))
        (run/name).write_bytes(result.stdout)
    limits=dict(rms_m=2.,max_m=10.,landing_m=10.,group_rel=2e-6)
    write(run/'budget_contract.json',dict(physical_limits=limits,external_limits=dict(rms_m=10.,max_m=50.,landing_m=50.,group_rel=1e-5),
          reference_fraction=.1,reconstruction_fraction=.1,sampling_fraction=.1,fd_last_change_fraction=.2,
          total_rule='propagation + reference + reconstruction + sampling <= physical limit',strict_multiplier=3.,
          sampling_steps_km=[1.,.5,.25],maximum_extra_sampling_levels=2,maximum_extra_integration_levels=2,
          dispersion_residual_limit=1e-10,dispersion_residual_source='production launch_state and projection tolerance',
          source_plan_sha256=sha(plan),frozen_before_AH=True))
    write(run/'case_registry.json',dict(schema='p31-physics-registry-v1',formal_case_count=26,cases=cases,media=media,original_inputs=originals))
    signature=inspect.signature(RefractiveIndexAH3DECEF)
    fd_defaults={key:signature.parameters[key].default for key in ('spatial_fd_scales','spatial_fd_stability_rtol','spatial_fd_stability_atol_n2_per_km','near_cutoff_n2')}
    manifest=dict(schema='p31-physics-manifest-v1',status='prepared_not_executed',formal_case_count=26,plan=str(plan),plan_sha256=sha(plan),
          run_purpose='preparation_selfcheck' if '_prepare_selfcheck_' in run.name else 'formal_physical_matrix_preparation',formal_matrix_execution_started=False,
          parent_manifest=str(parent_manifest),parent_manifest_sha256=sha(parent_manifest),plasma_coefficient=config.K_FP**2,coefficient=config.K_FP**2,
          plasma_coefficient_hex=float(config.K_FP**2).hex(),compatibility_mode=None,interpolation='radial_pchip',
          earth_radius_km=6371.,gyro_MHz=0.,collision_s=0.,integration_profiles=profiles[:3],candidate_integration_profiles=profiles[3:],
          projection_enabled=inspect.signature(trace_ray_ah3d_ecef).parameters['apply_dispersion_projection'].default,
          projection=inspect.signature(trace_ray_ah3d_ecef).parameters['apply_dispersion_projection'].default,
          local_projection_enabled=inspect.signature(integrate_fixed_path_ah3d_ecef).parameters['apply_dispersion_projection'].default,
          spatial_fd_base_km=.125,spatial_fd_steps_km=[.25,.125,.0625],spatial_fd_defaults=fd_defaults,
          analytic_generator=dict(path='tests/_run_a4_c1m_attribution.py',sha256=sha(ROOT/'tests/_run_a4_c1m_attribution.py'),function='_a2_profile',
                                  snapshot_path='source_before/tests/_run_a4_c1m_attribution.py',
                                  chapman_formula='N_peak*exp(0.5*(1-u-exp(-u))), u=(h-350)/100',chapman_peak_m3=607635874822.0494,
                                  smooth_formula='Chapman*x^3*(10+x*(-15+6*x)), x=clip((h-20)/60,0,1)',
                                  affine_formula='1e11+gradient*(h-100)',parameter_source='case_registry.media continuous_background'),
          max_path_km=5000.,h_stop_km=460.,local_path_length_km=100.,local_stop_on_physical_events=False,
          event_parameters=dict(departure_tolerance_km=1e-6,bisection_tolerance_km=1e-9,group_atol=1e-12),
          environment=dict(python=sys.version,executable=sys.executable,numpy=np.__version__,scipy=scipy.__version__,mpmath=mpmath.__version__),
          source_sha256=source_records,
          source_records=[dict(path=str(Path(row['path']).relative_to(ROOT)),sha256=row['sha256']) for row in source_records if Path(row['path']).is_relative_to(ROOT)],
          input_records=[dict(path=str(path.relative_to(run)),sha256=sha(path)) for path in sorted((run/'inputs').rglob('*')) if path.is_file()],
          git_head=(run/'git_head.txt').read_text().strip(),
          git_artifacts={name:sha(run/name) for name in ('git_head.txt','git_status.txt','git_diff.patch','git_diff_stat.txt')},
          case_registry_sha256=sha(run/'case_registry.json'),budget_contract_sha256=sha(run/'budget_contract.json'),exposure_audit_sha256=sha(run/'exposure_audit.json'),
          reference_topology_registered=False,production_modified=False,recommend_p4_physics=False)
    write(run/'manifest.json',manifest)
    return dict(status='prepared_not_executed',run=str(run),formal_case_count=26,source_files=len(source_records),
                reference_topology_registered=False,production_modified=False,recommend_p4_physics=False)
def main():
    p=argparse.ArgumentParser();p.add_argument('stage',choices=('prepare','boundary','reference','repair-regression','verify','physics-prepare'));p.add_argument('--run',required=True,type=Path);args=p.parse_args()
    try:result={'prepare':prepare,'boundary':boundary,'reference':reference,'repair-regression':repair,'verify':verify,'physics-prepare':physics_prepare}[args.stage](args.run.resolve())
    except Blocked as error:print(json.dumps(dict(stage=args.stage,status='blocked',reason=str(error),recommend_p4=False)));return 2
    except (ValueError,KeyError,OSError) as error:print(json.dumps(dict(stage=args.stage,status='rejected',reason=str(error),recommend_p4=False)));return 1
    print(json.dumps(dict(stage=args.stage,**result),indent=2));return 0
if __name__=='__main__':raise SystemExit(main())
