import time
import subprocess
from pathlib import Path
import sys
import re
import threading
from dataclasses import dataclass, field
from typing import Optional

_RUN_RE = re.compile(r"mol(\d+)")
_SUB_RE = re.compile(r"_sub_(\d+)_")


def extract_run_index(dat_file: Path) -> int | None:
    m = _RUN_RE.search(dat_file.name)
    return int(m.group(1)) if m else None


def extract_sub_index(dat_file: Path) -> int | None:
    m = _SUB_RE.search(dat_file.name)
    return int(m.group(1)) if m else None


@dataclass
class WatchConfig:
    watch_dir: Path
    scenario_root: Path
    backmap_script: Path
    python_exe: str = field(default_factory=lambda: sys.executable)
    poll_interval: float = 0.5
    max_backmap: int = 1
    stable_for: float = 1.0
    stable_poll: float = 0.2
    stable_timeout: float = 60.0
    defer_backmap_seconds: float = 0.0
    out_dir: Path | None = None
    overwrite: bool = False
    do_foxs: bool = False
    foxs_py: Optional[str] = None
    saxs_dat: Optional[Path] = None
    max_q: Optional[float] = None
    dat_glob: str = "*.dat"
    ignore_suffixes: tuple[str, ...] = (".tmp", ".part")
    backend: str = "modeller"
    cg2all_exec: Optional[str] = None
    disulfide_file: Optional[Path] = None


def wait_until_stable(path: Path, stable_for: float, poll: float, timeout: float) -> bool:
    start = time.time()
    last_size = None
    stable_start = None
    while True:
        if not path.exists():
            time.sleep(poll)
            if time.time() - start > timeout:
                return False
            continue
        size = path.stat().st_size
        if last_size is None or size != last_size:
            last_size = size
            stable_start = time.time()
        else:
            if stable_start is not None and (time.time() - stable_start) >= stable_for:
                return True
        if time.time() - start > timeout:
            return False
        time.sleep(poll)


def fingerprint_for_dat(cfg: WatchConfig, dat_file: Path) -> Path:
    sub_i = extract_sub_index(dat_file)
    if sub_i is None:
        return cfg.scenario_root / "fingerPrint1.dat"
    return cfg.scenario_root / f"fingerPrint{sub_i + 1}.dat"


def run_backmap(cfg: WatchConfig, dat_file: Path) -> int:
    run_i = extract_run_index(dat_file)
    run_dir = cfg.watch_dir / (f"allAtomRun{run_i}" if run_i is not None else "allAtomRunUnknown")
    run_dir.mkdir(parents=True, exist_ok=True)

    fp = fingerprint_for_dat(cfg, dat_file)
    if not fp.exists():
        print(f"[FAIL] Missing fingerprint for {dat_file.name}: {fp}", flush=True)
        return 3

    name = dat_file.stem
    cmd = [
        cfg.python_exe, str(cfg.backmap_script),
        "--coords", str(dat_file),
        "--fingerprint", str(fp),
        "--outdir", str(run_dir),
        "--name", name,
        "--backend", cfg.backend,
        "--scenario-root", str(cfg.scenario_root),
    ]

    if cfg.backend == "cg2all":
        if cfg.cg2all_exec is None:
            raise ValueError("cg2all backend selected but cg2all_exec not set")
        cmd += ["--cg2all-exec", str(cfg.cg2all_exec)]

    if cfg.disulfide_file is not None:
        cmd += ["--disulfide-file", str(cfg.disulfide_file)]

    if cfg.do_foxs:
        if not (cfg.foxs_py and cfg.saxs_dat and cfg.max_q is not None):
            raise ValueError("FoXS enabled but foxs_py/saxs_dat/max_q not set.")
        foxs_out = run_dir / "foxs_results.txt"
        cmd += [
            "--do-foxs",
            "--foxs-py", str(cfg.foxs_py),
            "--saxs", str(cfg.saxs_dat),
            "--max-q", str(cfg.max_q),
            "--foxs-out", str(foxs_out),
        ]

    print(f"[RUN ] {' '.join(cmd)}", flush=True)
    p = subprocess.run(cmd, capture_output=True, text=True)

    if p.returncode == 0:
        out_pdb = run_dir / f"{name}_AA.pdb"
        print(f"[DONE] {out_pdb}", flush=True)
    else:
        print(f"[FAIL] {dat_file.name} (exit {p.returncode})", flush=True)
        if p.stdout.strip():
            print("  stdout:", p.stdout.strip()[:500], flush=True)
        if p.stderr.strip():
            print("  stderr:", p.stderr.strip()[:500], flush=True)

    return p.returncode


class PollingWatcher:
    def __init__(self, cfg: WatchConfig):
        self.cfg = cfg
        self._stop = threading.Event()
        self._thread = None
        self._last_processed_mtime = {}
        self._inflight = set()
        self._sem = threading.Semaphore(cfg.max_backmap)
        self._start_time = time.time()
        self._activation_time = self._start_time + cfg.defer_backmap_seconds

    def _run_backmap_limited(self, dat_file: Path):
        if self._stop.is_set():
            return
        acquired = self._sem.acquire(timeout=0.5)
        if not acquired:
            self._inflight.discard(dat_file)
            return
        try:
            if self._stop.is_set():
                return
            print(f"[BACKMAP] starting {dat_file.name}", flush=True)
            run_backmap(self.cfg, dat_file)
        finally:
            self._inflight.discard(dat_file)
            self._sem.release()

    def start(self):
        if self._thread and self._thread.is_alive():
            print("[INFO] Watcher already running.")
            return
        self._stop.clear()
        self._thread = threading.Thread(target=self._loop, daemon=True)
        self._thread.start()
        print(f"[INFO] Watching {self.cfg.watch_dir} (poll={self.cfg.poll_interval}s)")

    def stop(self):
        self._stop.set()
        if self._thread:
            self._thread.join(timeout=2)
        print("[INFO] Watcher stopped.")

    def _loop(self):
        self.cfg.watch_dir.mkdir(parents=True, exist_ok=True)
        while not self._stop.is_set():
            for dat in sorted(self.cfg.watch_dir.glob(self.cfg.dat_glob)):
                if not (dat.name.startswith("mol") and dat.name.endswith("_xyz.dat")):
                    continue
                if not dat.is_file():
                    continue
                if dat.suffix.lower() != ".dat":
                    continue
                if any(str(dat).endswith(sfx) for sfx in self.cfg.ignore_suffixes):
                    continue
                if dat in self._inflight:
                    continue

                mtime = dat.stat().st_mtime
                if self._last_processed_mtime.get(dat) == mtime:
                    continue

                # Hard cutoff: only files created/written after the activation
                # threshold are ever eligible for backmapping.
                if mtime < self._activation_time:
                    self._last_processed_mtime[dat] = mtime
                    continue

                ok = wait_until_stable(dat, self.cfg.stable_for, self.cfg.stable_poll, self.cfg.stable_timeout)
                if not ok:
                    print(f"[WARN] Never stabilized: {dat.name}")
                    self._last_processed_mtime[dat] = mtime
                    continue

                if dat.stat().st_size == 0:
                    print(f"[WARN] Empty file, will retry later: {dat.name}", flush=True)
                    continue

                print(f"[WATCHER] running backmap for {dat.name}", flush=True)
                self._inflight.add(dat)
                threading.Thread(target=self._run_backmap_limited, args=(dat,), daemon=True).start()
                self._last_processed_mtime[dat] = dat.stat().st_mtime

            time.sleep(self.cfg.poll_interval)


_WATCHER = None


def start_watcher(watch_dir, scenario_root, backmap_script,
                  out_dir=None, poll_interval=0.5, overwrite=False,
                  max_backmap=1, defer_backmap_seconds=0.0,
                  do_foxs=False, foxs_py=None, saxs_dat=None, max_q=None,
                  backend="modeller", cg2all_exec=None, disulfide_file=None):
    global _WATCHER
    cfg = WatchConfig(
        watch_dir=Path(watch_dir),
        scenario_root=Path(scenario_root),
        backmap_script=Path(backmap_script),
        out_dir=Path(out_dir) if out_dir else None,
        poll_interval=poll_interval,
        overwrite=overwrite,
        max_backmap=max_backmap,
        defer_backmap_seconds=defer_backmap_seconds,
        do_foxs=do_foxs,
        foxs_py=str(foxs_py) if foxs_py else None,
        saxs_dat=Path(saxs_dat).resolve() if saxs_dat else None,
        max_q=max_q,
        backend=backend,
        cg2all_exec=cg2all_exec,
        disulfide_file=Path(disulfide_file).resolve() if disulfide_file else None,
    )
    _WATCHER = PollingWatcher(cfg)
    _WATCHER.start()
    return _WATCHER


def stop_watcher():
    global _WATCHER
    if _WATCHER is not None:
        _WATCHER.stop()
        _WATCHER = None


def main():
    import argparse
    import signal
    ap = argparse.ArgumentParser()
    ap.add_argument("--watch-dir", required=True)
    ap.add_argument("--scenario-root", required=True)
    ap.add_argument("--backmap-script", required=True)
    ap.add_argument("--poll", type=float, default=0.5)
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--max-backmap", type=int, default=1)
    ap.add_argument("--do-foxs", action="store_true")
    ap.add_argument("--foxs-py", default=None)
    ap.add_argument("--saxs", default=None)
    ap.add_argument("--max-q", type=float, default=None)
    ap.add_argument("--defer-backmap-seconds", type=float, default=0.0)
    ap.add_argument("--backend", choices=["modeller", "cg2all"], default="modeller")
    ap.add_argument("--cg2all-exec", default=None)
    ap.add_argument("--disulfide-file", default=None)
    args = ap.parse_args()

    if args.do_foxs and (args.saxs is None or args.max_q is None or args.foxs_py is None):
        ap.error("--do-foxs requires --foxs-py, --saxs, and --max-q")
    if args.backend == "cg2all" and args.cg2all_exec is None:
        ap.error("--backend cg2all requires --cg2all-exec")

    cfg = WatchConfig(
        watch_dir=Path(args.watch_dir).resolve(),
        scenario_root=Path(args.scenario_root).resolve(),
        backmap_script=Path(args.backmap_script).resolve(),
        poll_interval=args.poll,
        overwrite=args.overwrite,
        max_backmap=args.max_backmap,
        do_foxs=args.do_foxs,
        foxs_py=str(args.foxs_py) if args.foxs_py else None,
        saxs_dat=Path(args.saxs).resolve() if args.saxs else None,
        max_q=args.max_q,
        defer_backmap_seconds=args.defer_backmap_seconds,
        backend=args.backend,
        cg2all_exec=args.cg2all_exec,
        disulfide_file=Path(args.disulfide_file).resolve() if args.disulfide_file else None,
    )

    print("[WATCHER] started", flush=True)
    print("watch_dir      :", cfg.watch_dir, flush=True)
    print("scenario_root  :", cfg.scenario_root, flush=True)
    print("backmap_script :", cfg.backmap_script, flush=True)
    print("max_backmap    :", cfg.max_backmap, flush=True)
    print(f"[WATCHER] backmapping activates after {cfg.defer_backmap_seconds} s", flush=True)

    if not cfg.backmap_script.exists():
        raise FileNotFoundError(cfg.backmap_script)

    watcher = PollingWatcher(cfg)
    watcher.start()

    def _handle(sig, frame):
        print(f"[WATCHER] got signal {sig}, stopping...", flush=True)
        watcher.stop()
        raise SystemExit(0)

    signal.signal(signal.SIGINT, _handle)
    signal.signal(signal.SIGTERM, _handle)

    while True:
        time.sleep(1)


if __name__ == "__main__":
    main()
