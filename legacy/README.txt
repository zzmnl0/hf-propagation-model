Legacy coordinate source archive

This directory contains pre-C0R local Cartesian, earth-flattening, 2-D,
local_2p5d PE, and flux-tube implementations retained only to reproduce
historical results.

Production rules:
1. main.py, main2.py, and models/ must not import legacy.*.
2. No compatibility copy or fallback is kept in models/.
3. Legacy tests must import legacy.models explicitly.
4. A legacy capability may return to production only after migration to the
   geocentric/ECEF contract and independent acceptance.

Contents:
- legacy/models: local Cartesian, 2-D/2.5-D ray, P2P, PE, tube and radar code.
- legacy/analysis: analyses that depend on the archived coordinate contract.
- legacy/viz: plots tied to the archived 2-D background contract.
- legacy/_bench.py: historical benchmark entry.

Maintenance:
1. Shared current physics utilities may be imported explicitly as models.*;
   production modules must never import legacy.*.
2. Do not copy an archived module back into models/. Migrate its coordinate
   and data contract, add independent acceptance, then move the implementation.
3. Run tests/_verify_legacy_source_separation.py after layout changes.
