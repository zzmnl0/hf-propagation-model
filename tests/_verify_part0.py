"""Part 0 verification: imports, config sanity, and module stub check."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config as cfg
import utils
import numpy as np

print('BG_X shape:', cfg.BG_X.shape, '  range:', cfg.BG_X[0], '~', cfg.BG_X[-1], 'km')
print('BG_Z shape:', cfg.BG_Z.shape, '  range:', cfg.BG_Z[0], '~', cfg.BG_Z[-1], 'km')
print('k0 @ 10MHz :', round(utils.freq_to_k0(cfg.FREQ_MHZ), 2), 'km^-1')
print('lambda     :', utils.wavelength_m(cfg.FREQ_MHZ), 'm')
print('PE dz      :', cfg.PE['dz_m'], 'm')
lfs = utils.free_space_loss_dB(1169, cfg.FREQ_MHZ)
tau = utils.group_delay_ms(1169 * 2)
print(f'L_free (1169km) : {lfs:.1f} dB')
print(f'tau (2x1169km)  : {tau:.2f} ms  (1-hop approx round-trip)')

# check all modules are importable
from models import ionosphere_model, ray_tracer, point_to_point, es_model, pe_propagator, hybrid_model
print()
print('All 6 model modules imported OK.')

# check implementation status
try:
    ionosphere_model.IonosphereModel().build_Ne_field(cfg.BG_X[:3], cfg.BG_Z[:3])
    print('  ionosphere_model.IonosphereModel -> IMPLEMENTED (Part 1 done)')
except NotImplementedError:
    print('  ionosphere_model.IonosphereModel -> stub (Part 1 pending)')

for name, fn in [
    ('ray_tracer.shoot_rays_fan',       lambda: ray_tracer.shoot_rays_fan(None, None)),
    ('point_to_point.find_all_rays_p2p',lambda: point_to_point.find_all_rays_p2p(None, None, None, None)),
    ('es_model.EsLayerModel.classify',  lambda: es_model.EsLayerModel().classify(10.0)),
    ('pe_propagator.PEPropagator',      lambda: pe_propagator.PEPropagator().propagate(None, None, None, None)),
    ('hybrid_model.HybridPropagationModel', lambda: hybrid_model.HybridPropagationModel({}, {}).compute()),
]:
    try:
        fn()
        print(f'  {name} -> IMPLEMENTED')
    except NotImplementedError:
        print(f'  {name} -> stub (pending)')
    except Exception:
        pass   # other errors (e.g. wrong args) also mean stub is reachable

print()
print('Part 0 scaffolding verified.')
