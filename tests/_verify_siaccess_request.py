"""One bounded SIACCESS engineering import: 20 real factory constructions, no physics."""

import argparse
import dataclasses
import sys
from pathlib import Path


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True)
    args=parser.parse_args();folder=Path(args.output).resolve();sys.path.insert(0,str(folder))
    from common import ROOT, HERE, TASK, load, save, sha, guard, now
    guard('E');assert folder==HERE
    save(HERE/'factory_production_import.json',{'task':TASK,'utc':now(),'production_import_process':1})
    sys.path.insert(0,str(ROOT))
    import numpy as np
    import config as cfg
    import models.unified_hybrid_pipeline as u
    import models.ray_tracer_ah_3d_ecef as ah
    from models.ne_provider_geocentric import GeocentricGridNeProvider3D
    from threadpoolctl import threadpool_info
    assert all(sha(p)==v for p,v in load(HERE/'source_registry.json')['hashes'].items())
    denied_calls=[];constructed=[];leaves=[];controls=[]
    def require(ok,message):
        if not ok:raise ValueError(message)
    def denied(*args,**kwargs):
        denied_calls.append('roots_launch_RHS_integral');raise RuntimeError('engineering physics forbidden')
    for name in ['roots','root','root_id','n2_raw']:
        if hasattr(ah.RefractiveIndexAH3DECEF,name):setattr(ah.RefractiveIndexAH3DECEF,name,denied)
    for name in ['launch_state_ah3d_ecef','ah3d_ecef_rhs','integrate_fixed_path_ah3d_ecef','trace_ray_ah3d_ecef']:
        setattr(ah,name,denied)
    real_class=u.RefractiveIndexAH3DECEF
    def counted(*args,**kwargs):
        require(len(constructed)<20,'constructor cap')
        item={'attempt':len(constructed)+1,'wave_mode':kwargs.get('wave_mode'),'profile':kwargs.get('integration_profile'),
            'requested_coefficient':str(kwargs.get('plasma_coefficient')),'compatibility':kwargs.get('compatibility_mode')}
        constructed.append(item)
        try:model=real_class(*args,**kwargs)
        except BaseException as exc:
            item['exception']=repr(exc);raise
        item.update({'actual_coefficient':model.plasma_coefficient,'actual_compatibility':model.compatibility_mode,
            'actual_zero_B':bool(np.all(model.geomagnetic_provider.b_ecef_T==0.))})
        return model
    u.RefractiveIndexAH3DECEF=counted
    cfg.GEOMAG['fH_MHz']=0.
    provider=GeocentricGridNeProvider3D(np.full((2,2,3),1e11),[20.,40.],[110.,130.],[80.,130.,180.],6371.,interpolation='radial_hermite')
    base={'scenario':TASK+'_engineering_only','provider':provider,'coordinate_model':ah.COORDINATE_MODEL,
        'earth_radius_km':6371.,'tx_geocentric':(30.,120.,100.),'rx_geocentric':(30.,120.,0.),'bearing_deg':90.,'freq_MHz':10.,
        'x_km':np.array([0.,1.]),'z_km':np.array([100.,110.]),'make_pd':False}
    SI=80.61638587963631
    choices=[('omitted',False,{},cfg.AH3D_PLASMA_COEFFICIENT_SI),('SI',False,{'plasma_coefficient':SI},SI),
        ('None_compat',True,{'plasma_coefficient':None},80.61638494312905),('SI_compat',True,{'plasma_coefficient':SI},SI),
        ('zero',False,{'plasma_coefficient':0.},None),('negative',False,{'plasma_coefficient':-1.},None),
        ('NaN',False,{'plasma_coefficient':float('nan')},None),('Inf',False,{'plasma_coefficient':float('inf')},None)]
    require(dataclasses.fields(u.UnifiedPropagationRequest)[-1].name=='plasma_coefficient','field not appended')
    for label,compat,kwargs,expected in choices:
        cfg.UPDATE2['PHARLAP_RADIAL_B0']['enable']=compat
        request=u.UnifiedPropagationRequest(**base,**kwargs)
        require(request.plasma_coefficient is None if not kwargs else True,'omission default')
        try:model=u._ah_model(request,'O','research')
        except ValueError:
            require(expected is None,'valid factory rejected')
        else:
            require(expected is not None and model.plasma_coefficient==expected,'factory coefficient/illegal acceptance')
            require(model.compatibility_mode==('pharlap_radial_b0' if compat else None),'compat identity')
        controls.append({'id':label,'decision':'PASS'})
    cfg.UPDATE2['PHARLAP_RADIAL_B0']['enable']=False
    cfg.UPDATE2['RAY_FAN_INCLUDE_X']=True
    require(not cfg.ABSORPTION.get('enable',False) and not cfg.UPDATE2['PE_COUPLING'].get('enable',False),'archived config active')
    # Only expensive propagation/power leaves are replaced; public dispatch, seeds and factory remain real.
    lat,lon=np.deg2rad([30.,120.]);landing=6371.*np.array([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)])
    def ray(model,elevation=35.,bearing=90.):
        return {'termination_status':'ground','landing_lat_lon_height':(30.,120.,0.),
            'trajectory_ecef_km':np.array([landing+np.array([0.,0.,100.]),landing]),'wave_mode':model.wave_mode,
            'elevation_tx_deg':elevation,'bearing_tx_deg':bearing,'apex_height_km':100.,'tau_ms':1.,'group_path_km':300.}
    def trace_leaf(model,**kwargs):
        leaves.append({'leaf':'trace','wave_mode':model.wave_mode,'profile':kwargs['profile'],'coefficient':model.plasma_coefficient})
        return ray(model,kwargs['elevation_deg'],kwargs['bearing_deg'])
    def p2p_leaf(model,**kwargs):
        require(len(kwargs['seeds'])==1,'real fan seed not reachable')
        leaves.append({'leaf':'p2p','wave_mode':model.wave_mode,'profile':kwargs['profile'],'coefficient':model.plasma_coefficient,'seeds':kwargs['seeds']})
        return [ray(model)]
    u.trace_ray_ah3d_ecef=trace_leaf;u.find_seeded_rays_ah3d_ecef=p2p_leaf
    u._add_power=lambda *args,**kwargs:None
    u._add_scatter_diagnostics=lambda *args,**kwargs:None
    u._apply_p3_energy_closure=lambda *args,**kwargs:None
    # Avoid actual provider sampling in a route-only fixture; no Ne derivative/root is exercised.
    provider.sample_link_plane=lambda *args,**kwargs:np.full((2,2),1e11)
    for label,kwargs,expected in [('default',{},cfg.AH3D_PLASMA_COEFFICIENT_SI),('SI',{'plasma_coefficient':SI},SI)]:
        start=len(constructed);leaf_start=len(leaves)
        request=u.UnifiedPropagationRequest(**base,**kwargs,
            rt_config={**cfg.RT_3D,'beta_min':35.,'beta_max':35.,'n_fan':1},
            p2p_config={**cfg.P2P_3D,'solve_profile':'research','final_retrace_profile':'reference','candidate_nearest':1})
        clone=u.UnifiedPropagationRequest(**{**request.__dict__,'apply_scatter_to_pd':False})
        require(clone.plasma_coefficient==request.plasma_coefficient,'request.__dict__ copy lost coefficient')
        result=u.run_unified_pipeline(request)
        require(len(constructed)==start+6 and len(result.rays3d)==2 and len(result.modes)==2,'route coverage count')
        observed=[(r['wave_mode'],r['profile']) for r in constructed[start:]]
        require(sorted(observed)==sorted([(m,p) for m in ['O','X'] for p in ['smoke','research','reference']]),'fan/solve/retrace mode/profile coverage')
        require(all(r['actual_coefficient']==expected and r['actual_zero_B'] for r in constructed[start:]),'route coefficient/field drift')
        require(len(leaves[leaf_start:])==6 and sum(r['leaf']=='p2p' for r in leaves[leaf_start:])==2,'expensive leaves/roots not reached')
        controls.append({'id':'public_route_'+label,'decision':'PASS','constructors':6,'modes':['O','X'],'profiles':observed,
            'leaf_scope':'synthetic propagation/power only; factory/dispatch/seed selection real'})
    require(len(constructed)==20 and not denied_calls,'engineering count/forbidden calls')
    require(all(p['num_threads']==1 for p in threadpool_info()),'thread pool count')
    modules=[]
    for name,module in sorted(sys.modules.items()):
        path=getattr(module,'__file__',None)
        if path and Path(path).is_file():modules.append({'module':name,'path':str(Path(path).resolve()),'sha256':sha(path)})
    save(HERE/'factory_result.json',{'task':TASK,'decision':'PASS','constructors':constructed,'controls':controls,
        'constructor_attempts':20,'production_import_processes':1,'forbidden_calls':denied_calls,
        'roots_calls':0,'launch_calls':0,'RHS_calls':0,'integrals':0,'leaf_calls':leaves,'module_files':modules,'threadpools':threadpool_info(),
        'claim':'parameter parsing and public routing only; synthetic leaves are not trajectories'})
    print('SIACCESS_ENGINEERING_PASS constructors=20 roots=0 launch=0 RHS=0 integral=0')


if __name__=='__main__':main()
