"""Run the five independent report scenes concurrently with bounded P2P work."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output2"
SCENES = ("baseline", "tid", "spreadf", "es", "bubble")


def run() -> int:
    processes = []
    streams = []
    for scene in SCENES:
        log = (OUT / "run_unified_ox_p2p_{}.log".format(scene)).open(
            "w", encoding="utf-8")
        streams.append(log)
        command = [
            sys.executable, "-u", str(ROOT / "main2.py"),
            "--scenario", scene, "--ray-fan-x", "--p2p",
            "--p2p-seeds", "1", "--p2p-profile", "smoke",
            "--p2p-final-profile", "research", "--p2p-max-nfev", "5",
            "--csv",
        ]
        processes.append((scene, subprocess.Popen(
            command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)))
        print("started {} pid={}".format(scene, processes[-1][1].pid))
    failures = []
    for scene, process in processes:
        returncode = process.wait()
        print("finished {} exit={}".format(scene, returncode))
        if returncode:
            failures.append(scene)
    for stream in streams:
        stream.close()
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(run())
