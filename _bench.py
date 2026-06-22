"""Bench: time scenario 0 P2P with fan-seeded warm start."""
import time, sys, os
import matplotlib; matplotlib.use('Agg')
sys.path.insert(0, os.path.dirname(__file__))
import config as cfg
from models.ionosphere_model import IonosphereModel
from models.ray_tracer_3d import RefractiveIndex3D, shoot_rays_fan_3d
from models.point_to_point_3d import find_all_rays_3d

iono = IonosphereModel()
Ne_2d, _ = iono.build_Ne_field(cfg.BG_X, cfg.BG_Z)
n3d = RefractiveIndex3D(Ne_2d, cfg.BG_X, cfg.BG_Z, cfg.FREQ_MHZ, earth_flat=True)

t1 = time.time()
rays3d = shoot_rays_fan_3d(cfg.TX_POS, n3d, cfg.FREQ_MHZ, cfg.RT_3D)
t2 = time.time()
print('[fan] {} rays: {:.2f} s'.format(len(rays3d), t2-t1))

t3 = time.time()
modes = find_all_rays_3d(cfg.TX_POS, cfg.RX_POS, n3d, cfg.FREQ_MHZ, cfg.P2P_3D,
                          rays_fan=rays3d)
t4 = time.time()
print('[P2P fan-seeded] {} modes: {:.2f} s'.format(len(modes), t4-t3))
print('[P2P cold start] (skipped)')
print('Total: {:.2f} s'.format(t4-t1))
print('OK')
