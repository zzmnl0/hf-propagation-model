"""A1 check that O/X identities remain continuous during a frequency scan."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from models.dispersion_ah import appleton_hartree_roots, select_root


def run():
    previous = {"O": None, "X": None}
    histories = {"O": [], "X": []}
    for freq_MHz in np.linspace(8.0, 14.0, 121):
        roots = appleton_hartree_roots(
            2.0 * np.pi * freq_MHz * 1e6, 4.0e11, [1.0, 0.2, 0.3],
            [0.45, -0.15, -0.88], 1.2)
        for branch in ("O", "X"):
            root = select_root(roots, branch, previous[branch])
            histories[branch].append(root)
            previous[branch] = root.n2
    ok = all(root.branch_id == branch and root.status != "branch_lost"
             for branch, values in histories.items() for root in values)
    ok = ok and min(root.branch_gap for root in histories["O"]) > 0.0
    print(" frequency_continuity {}".format("OK" if ok else "FAIL"))
    print("A1 frequency continuity: {} fail(s)".format(0 if ok else 1))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(run())
