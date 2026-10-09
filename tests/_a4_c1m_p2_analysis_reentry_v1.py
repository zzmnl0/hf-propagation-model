"""Read-only P2 v2 identity reentry; preserve the frozen core and its CLI failure."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import types
sys.dont_write_bytecode = True


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(run: Path) -> None:
    run=run.resolve();receipt=run/"analysis/p2_reentry_v1"
    reentry=json.loads((receipt/"manifest.json").read_text(encoding="utf-8"))
    own=reentry["source"]
    if sha(Path(__file__).resolve())!=own["sha256"] or sha(receipt/own["snapshot"])!=own["sha256"]:
        raise ValueError("Frozen reentry source changed")
    analysis=run/"analysis/p2_independent_v2"
    binding=json.loads((analysis/"analysis_manifest.json").read_text(encoding="utf-8"))
    source=binding["sources"]["core"];snapshot=analysis/source["snapshot"]
    if sha(snapshot)!=source["sha256"]:raise ValueError("Frozen v2 core changed")
    code=Path(__file__).resolve().parent.parent
    sys.path[:0]=[str(code/"tests"),str(code)]
    module=types.ModuleType("_a4_p2_v2_readonly_reentry");module.__file__=str(snapshot)
    sys.modules[module.__name__]=module
    exec(compile(snapshot.read_text(encoding="utf-8"),str(snapshot),"exec"),module.__dict__)
    # The original read-only CLI fails only because its final print lacks this
    # standard-library dependency. Supply it in this isolated reference module.
    module.json=json
    module.verify_supplement(run)


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--run",required=True)
    args=parser.parse_args();verify(Path(args.run))
