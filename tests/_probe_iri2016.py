"""Probe iri2016 API"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from datetime import datetime
import numpy as np

import iri2016
print("dir:", [x for x in dir(iri2016) if not x.startswith('_')])

dt  = datetime(2020, 6, 1, 12, 0)
lat, lon = 32.5, 120.0

# iri2016 takes (min, max, step) not an array
altkmrange = (60, 600, 2)  # min, max, step

try:
    res = iri2016.IRI(dt, altkmrange, lat, lon)
    print("IRI() type:", type(res))
    if hasattr(res, 'data_vars'):
        print("data_vars:", list(res.data_vars))
        if 'ne' in res.data_vars:
            ne = res['ne'].values
            print("ne shape:", ne.shape)
            print("ne peak:", ne.max(), "m^-3  at alt index", ne.argmax())
            alts = res.coords['alt_km'].values if 'alt_km' in res.coords else np.arange(60,602,2)
            print("altitude at peak:", alts[ne.argmax()], "km")
        # show all available variables
        for k in list(res.data_vars)[:6]:
            print(f"  [{k}]  shape={res[k].values.shape}  dtype={res[k].values.dtype}")
    print("coords:", list(res.coords))
except Exception as e:
    print("IRI() error:", e)
    import traceback; traceback.print_exc()
