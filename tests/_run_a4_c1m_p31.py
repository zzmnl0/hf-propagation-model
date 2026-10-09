"""Read-only stage dispatch for the frozen P3.1 run; never silently restart experiments."""
from __future__ import annotations
import argparse, hashlib, json, os, shutil, subprocess, sys
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
def main():
    p=argparse.ArgumentParser();p.add_argument('stage',choices=('prepare','boundary','reference','repair-regression','verify'));p.add_argument('--run',required=True,type=Path);args=p.parse_args()
    try:result={'prepare':prepare,'boundary':boundary,'reference':reference,'repair-regression':repair,'verify':verify}[args.stage](args.run.resolve())
    except Blocked as error:print(json.dumps(dict(stage=args.stage,status='blocked',reason=str(error),recommend_p4=False)));return 2
    except (ValueError,KeyError,OSError) as error:print(json.dumps(dict(stage=args.stage,status='rejected',reason=str(error),recommend_p4=False)));return 1
    print(json.dumps(dict(stage=args.stage,**result),indent=2));return 0
if __name__=='__main__':raise SystemExit(main())
