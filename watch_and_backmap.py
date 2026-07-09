import time
import os
import signal
import subprocess
import shlex
from pathlib import Path
import sys
import re
import threading
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

_RUN_RE = re.compile(r"mol(\d+)")
_SUB_RE = re.compile(r"_sub_(\d+)_")

_STEP_RE = re.compile(r"^mol(\d+)_sub_(\d+)_step_(\d+)_xyz\.dat$")
_INITIAL_RE = re.compile(r"^mol(\d+)_sub_(\d+)_initial_xyz\.dat$")
_END_RE = re.compile(r"^mol(\d+)_sub_(\d+)_end_xyz\.dat$")


def extract_run_index(dat_file: Path) -> int | None:
    m = _RUN_RE.search(dat_file.name)
    return int(m.group(1)) if m else None


def extract_sub_index(dat_file: Path) -> int | None:
    m = _SUB_RE.search(dat_file.name)
    return int(m.group(1)) if m else None


def parse_prediction_key(dat_file: Path):
    """
    Parse Carbonara prediction filenames.

    Returns
    -------
    tuple | None
        (run_i, sub_i, tag), where tag is one of:
          ("step", step_i), ("initial", 0), ("end", large_int)

    Examples
    --------
    mol61_sub_0_step_3_xyz.dat -> (61, 0, ("step", 3))
    mol61_sub_1_step_3_xyz.dat -> (61, 1, ("step", 3))
    mol61_sub_0_initial_xyz.dat -> (61, 0, ("initial", 0))
    mol61_sub_0_end_xyz.dat -> (61, 0, ("end", 10**12))
    """
    name = dat_file.name

    m = _STEP_RE.match(name)
    if m:
        return int(m.group(1)), int(m.group(2)), ("step", int(m.group(3)))

    m = _INITIAL_RE.match(name)
    if m:
        return int(m.group(1)), int(m.group(2)), ("initial", 0)

    m = _END_RE.match(name)
    if m:
        return int(m.group(1)), int(m.group(2)), ("end", 10**12)

    return None


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

    # Number of mixture components / structures expected per prediction state.
    # For ordinary monomer fitting this is 1.
    no_structures: int = 1

    # FoXS
    do_foxs: bool = False
    foxs_py: Optional[str] = None
    saxs_dat: Optional[Path] = None
    max_q: Optional[float] = None

    # File patterns
    dat_glob: str = "*.dat"
    ignore_suffixes: tuple[str, ...] = (".tmp", ".part")

    # Backmapping backend
    backend: str = "modeller"
    cg2all_exec: Optional[str] = None
    disulfide_file: Optional[Path] = None

    # Optional batch-screening mode. Default is maximal exploration:
    # leave every predictStructureQvary process running to maxNoFitSteps.
    terminate_on_foxs: bool = False
    terminate_threshold: Optional[float] = None
    terminate_confirmation_count: int = 1


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


def infer_sub_index_base(sub_indices: list[int], no_structures: int) -> int:
    """
    Infer whether Carbonara sub indices are 0-based or 1-based.

    The current codebase has historically used sub_0 -> fingerPrint1, but some
    mixture examples are naturally described as sub_1..sub_N. This keeps the
    watcher robust to either convention.
    """
    subs = sorted(sub_indices)
    if subs == list(range(no_structures)):
        return 0
    if subs == list(range(1, no_structures + 1)):
        return 1
    # Fallback: preserve the existing historical behaviour.
    return 0


def fingerprint_for_dat(cfg: WatchConfig, dat_file: Path, sub_index_base: int = 0) -> Path:
    sub_i = extract_sub_index(dat_file)
    if sub_i is None:
        return cfg.scenario_root / "fingerPrint1.dat"

    if sub_index_base == 1:
        fp_num = sub_i
    else:
        fp_num = sub_i + 1

    return cfg.scenario_root / f"fingerPrint{fp_num}.dat"


def run_backmap(cfg: WatchConfig, dat_file: Path, sub_index_base: int = 0, do_foxs_override: Optional[bool] = None) -> int:
    run_i = extract_run_index(dat_file)
    run_dir = cfg.watch_dir / (f"allAtomRun{run_i}" if run_i is not None else "allAtomRunUnknown")
    run_dir.mkdir(parents=True, exist_ok=True)

    fp = fingerprint_for_dat(cfg, dat_file, sub_index_base=sub_index_base)
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

    do_foxs = cfg.do_foxs if do_foxs_override is None else bool(do_foxs_override)

    # Safety guard: for mixture/ensemble runs, FoXS must be evaluated once at
    # group level after all components are backmapped. Never let the per-file
    # backmap CLI write individual component entries to foxs_results.txt.
    if cfg.no_structures > 1:
        do_foxs = False

    print(f"[BACKMAP] {dat_file.name}: per_file_foxs={do_foxs}", flush=True)

    if do_foxs:
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
            print("  stdout:", p.stdout.strip()[:1000], flush=True)
        if p.stderr.strip():
            print("  stderr:", p.stderr.strip()[:1000], flush=True)

    return p.returncode


def aa_pdb_for_dat(cfg: WatchConfig, dat_file: Path) -> Path:
    run_i = extract_run_index(dat_file)
    run_dir = cfg.watch_dir / (f"allAtomRun{run_i}" if run_i is not None else "allAtomRunUnknown")
    return run_dir / f"{dat_file.stem}_AA.pdb"


def read_single_foxs_score_for_dat(cfg: WatchConfig, dat_file: Path) -> Optional[float]:
    """
    Read the FoXS chi^2 that backmap_cli.py writes for one ordinary
    single-structure prediction. Mixture/ensemble scores are handled by
    run_foxs_mixture_group().
    """
    run_i = extract_run_index(dat_file)
    if run_i is None:
        return None
    run_dir = cfg.watch_dir / f"allAtomRun{run_i}"
    summary = run_dir / "foxs_results.txt"
    aa = aa_pdb_for_dat(cfg, dat_file)
    if not summary.exists():
        return None

    aa_str = str(aa)
    aa_name = aa.name
    try:
        lines = summary.read_text().splitlines()
    except Exception:
        return None

    for line in reversed(lines):
        parts = line.split()
        if len(parts) < 2:
            continue
        # Current backmap_cli writes: <absolute/relative AA pdb path> <chi2|ERROR>
        if parts[0] == aa_str or parts[0].endswith("/" + aa_name) or aa_name in parts[0]:
            try:
                return float(parts[-1])
            except Exception:
                return None
    return None


def _load_numeric_table(path: Path):
    try:
        data = np.loadtxt(path)
    except Exception:
        return None
    if data.ndim == 1:
        data = data[None, :]
    if data.shape[0] < 2 or data.shape[1] < 2:
        return None
    return data


def _load_experimental_saxs(saxs_dat: Path, max_q: float | None = None):
    data = _load_numeric_table(Path(saxs_dat))
    if data is None:
        raise ValueError(f"Could not read numeric SAXS data from {saxs_dat}")
    q = data[:, 0]
    I = data[:, 1]
    if data.shape[1] >= 3:
        sigma = data[:, 2]
    else:
        sigma = np.ones_like(I)
    sigma = np.where(np.asarray(sigma, dtype=float) <= 0, 1.0, sigma)
    if max_q is not None:
        mask = q <= float(max_q)
        q, I, sigma = q[mask], I[mask], sigma[mask]
    return q.astype(float), I.astype(float), sigma.astype(float)


def _extract_calc_curve(profile_path: Path, q_target: np.ndarray):
    data = _load_numeric_table(profile_path)
    if data is None:
        return None
    q = data[:, 0].astype(float)
    # FoXS-like fitted curves are usually q, I_exp, I_calc; if more columns
    # exist, use the last column as the calculated intensity.
    Icalc = data[:, -1].astype(float)
    order = np.argsort(q)
    q = q[order]
    Icalc = Icalc[order]
    if len(q) != len(q_target) or np.max(np.abs(q - q_target)) > 1e-8:
        Icalc = np.interp(q_target, q, Icalc)
    return Icalc


def _find_foxs_profile(outdir: Path, aa_pdb: Path, started_at: float):
    # pyFoXS may choose slightly different extensions/names depending on version.
    # Prefer newly written numeric files containing the AA stem, excluding our logs/results.
    candidates = []
    for pat in ("*.dat", "*.fit", "*.txt"):
        candidates.extend(outdir.glob(pat))
    filtered = []
    for f in candidates:
        if f.name.endswith("_foxs.log"):
            continue
        if f.name in {"foxs_results.txt", "foxs_mixture_results.txt"}:
            continue
        if f.stat().st_mtime + 1e-6 < started_at:
            continue
        if aa_pdb.stem not in f.stem:
            continue
        if _load_numeric_table(f) is not None:
            filtered.append(f)
    if not filtered:
        return None
    return max(filtered, key=lambda x: x.stat().st_mtime)


def _run_foxs_for_profile(cfg: WatchConfig, aa_pdb: Path, outdir: Path):
    if not (cfg.foxs_py and cfg.saxs_dat and cfg.max_q is not None):
        raise ValueError("FoXS enabled but foxs_py/saxs_dat/max_q not set.")

    started_at = time.time()
    cmd = shlex.split(str(cfg.foxs_py)) + [
        str(aa_pdb),
        str(cfg.saxs_dat),
        "--max_q", str(cfg.max_q),
    ]
    log_path = outdir / f"{aa_pdb.stem}_foxs.log"
    p = subprocess.run(cmd, cwd=str(outdir), capture_output=True, text=True)
    log_path.write_text(
        f"$ {' '.join(cmd)}\n\n=== STDOUT ===\n{p.stdout}\n\n=== STDERR ===\n{p.stderr}\n"
    )
    if p.returncode != 0:
        return None, log_path, p.returncode
    profile = _find_foxs_profile(outdir, aa_pdb, started_at)
    return profile, log_path, p.returncode


def _fit_simplex_weights(component_profiles: np.ndarray, I_exp: np.ndarray, sigma: np.ndarray):
    """
    Fit non-negative weights summing to one, plus a single global scale factor.
    component_profiles has shape (n_components, n_q).
    """
    try:
        from scipy.optimize import minimize
    except Exception as exc:
        raise ImportError("scipy is required for mixture FoXS weight optimisation") from exc

    comps = np.asarray(component_profiles, dtype=float)
    y = np.asarray(I_exp, dtype=float)
    sig = np.asarray(sigma, dtype=float)
    sig = np.where(sig <= 0, 1.0, sig)
    n = comps.shape[0]

    def best_scale(mix):
        wt = 1.0 / (sig * sig)
        denom = np.sum(wt * mix * mix)
        if denom <= 0:
            return 1.0
        return float(np.sum(wt * y * mix) / denom)

    def objective(w):
        mix = np.dot(w, comps)
        c = best_scale(mix)
        r = (c * mix - y) / sig
        return float(np.mean(r * r))

    w0 = np.ones(n, dtype=float) / float(n)
    cons = ({"type": "eq", "fun": lambda w: np.sum(w) - 1.0},)
    bounds = [(0.0, 1.0)] * n
    res = minimize(objective, w0, method="SLSQP", bounds=bounds, constraints=cons)
    if not res.success:
        # Still return the best point SLSQP found; this is usually useful in monitoring.
        print(f"[WARN] mixture weight optimisation did not fully converge: {res.message}", flush=True)
    w = np.clip(np.asarray(res.x, dtype=float), 0.0, 1.0)
    s = float(np.sum(w))
    w = w / s if s > 0 else w0
    mix = np.dot(w, comps)
    scale = best_scale(mix)
    chi2 = objective(w)
    return w, scale, chi2


def _group_label(group_key):
    run_i, tag = group_key
    kind, val = tag
    if kind == "step":
        return f"mol{run_i}_step_{val}"
    return f"mol{run_i}_{kind}"


def _update_mixture_summary(summary_file: Path, label: str, chi2: float, weights, scale: float, pdbs):
    summary_file.parent.mkdir(parents=True, exist_ok=True)
    old = summary_file.read_text().splitlines() if summary_file.exists() else []
    old = [line for line in old if not line.startswith(label + " ")]
    wtxt = ",".join(f"{float(w):.6g}" for w in weights)
    pdbtxt = ",".join(str(p) for p in pdbs)
    old.append(f"{label} chi2={chi2:.8g} scale={scale:.8g} weights={wtxt} pdbs={pdbtxt}")
    summary_file.write_text("\n".join(old) + "\n")


def run_foxs_mixture_group(cfg: WatchConfig, group_key, files: list[Path]) -> tuple[int, Optional[float]]:
    """
    Score a grouped mixture state by running FoXS on each all-atom component,
    then fitting best non-negative mixture weights that sum to one.
    """
    if not cfg.do_foxs:
        return 0, None
    if len(files) <= 1:
        return 0, None

    run_i = group_key[0]
    outdir = cfg.watch_dir / f"allAtomRun{run_i}"
    aa_pdbs = [aa_pdb_for_dat(cfg, f) for f in files]
    missing = [p for p in aa_pdbs if not p.exists()]
    if missing:
        print(f"[FOXS-MIX] missing AA PDBs: {', '.join(str(p) for p in missing)}", flush=True)
        return 4, None

    q_exp, I_exp, sigma = _load_experimental_saxs(Path(cfg.saxs_dat), cfg.max_q)
    curves = []
    profiles = []
    for aa in aa_pdbs:
        profile, log_path, rc = _run_foxs_for_profile(cfg, aa, outdir)
        if rc != 0 or profile is None:
            print(f"[FOXS-MIX] FoXS failed/no profile for {aa.name}; see {log_path}", flush=True)
            return 5, None
        Icalc = _extract_calc_curve(profile, q_exp)
        if Icalc is None:
            print(f"[FOXS-MIX] Could not parse FoXS curve: {profile}", flush=True)
            return 6, None
        curves.append(Icalc)
        profiles.append(profile)

    weights, scale, chi2 = _fit_simplex_weights(np.vstack(curves), I_exp, sigma)
    label = _group_label(group_key)

    # Write mixture curve for inspection.
    mix = scale * np.dot(weights, np.vstack(curves))
    curve_file = outdir / f"{label}_foxs_mixture_fit.dat"
    np.savetxt(curve_file, np.column_stack([q_exp, I_exp, sigma, mix]),
               header="q I_exp sigma I_foxs_mixture_fit")

    summary_file = outdir / "foxs_mixture_results.txt"
    _update_mixture_summary(summary_file, label, chi2, weights, scale, aa_pdbs)
    print(f"[FOXS-MIX] {label} chi2={chi2:.6g} weights={','.join(f'{w:.4g}' for w in weights)}", flush=True)
    return 0, float(chi2)


class PollingWatcher:
    def __init__(self, cfg: WatchConfig):
        self.cfg = cfg
        self._stop = threading.Event()
        self._thread = None
        self._processed_groups = {}   # group_key -> tuple[(path, mtime), ...]
        self._completed_groups = set()
        self._inflight_groups = set()
        self._sem = threading.Semaphore(cfg.max_backmap)
        self._start_time = time.time()
        self._activation_time = self._start_time + cfg.defer_backmap_seconds
        self._good_foxs_counts = {}  # run_i -> number of qualifying FoXS scores
        self._termination_lock = threading.Lock()

    def _pid_file_for_run(self, run_i: int) -> Path:
        return self.cfg.watch_dir / f"run{run_i}.pid"

    def _stopped_file_for_run(self, run_i: int) -> Path:
        return self.cfg.watch_dir / f"run{run_i}.stopped"

    def _stop_failed_file_for_run(self, run_i: int) -> Path:
        return self.cfg.watch_dir / f"run{run_i}.stop_failed"

    def _run_is_stopped(self, run_i: int) -> bool:
        return self._stopped_file_for_run(run_i).exists()

    @staticmethod
    def _pid_is_alive(pid: int) -> bool:
        try:
            os.kill(pid, 0)
            return True
        except ProcessLookupError:
            return False
        except PermissionError:
            # It exists, but we cannot signal it. Treat as alive so we do not
            # accidentally report success.
            return True

    @staticmethod
    def _pid_cmdline(pid: int) -> str:
        proc_cmd = Path(f"/proc/{pid}/cmdline")
        try:
            raw = proc_cmd.read_bytes()
        except Exception:
            return ""
        return raw.replace(b"\x00", b" ").decode(errors="ignore").strip()

    def _write_stop_record(self, run_i: int, text: str) -> None:
        path = self._stopped_file_for_run(run_i)
        try:
            path.write_text(text.rstrip() + "\n")
        except Exception as exc:
            print(f"[STOP] Could not write {path}: {exc}", flush=True)

    def _write_stop_failed_record(self, run_i: int, text: str) -> None:
        path = self._stop_failed_file_for_run(run_i)
        try:
            path.write_text(text.rstrip() + "\n")
        except Exception as exc:
            print(f"[STOP] Could not write {path}: {exc}", flush=True)

    def _terminate_predictor_for_run(self, run_i: int, label: str, chi2: float, count: int) -> bool:
        """Terminate exactly one predictStructureQvary process using run<i>.pid."""
        pid_file = self._pid_file_for_run(run_i)
        if not pid_file.exists():
            msg = (
                f"No PID file for run {run_i}: {pid_file}\n"
                f"Wanted to stop after {label}, FoXS chi2={chi2:.8g}, qualifying_count={count}\n"
            )
            print(f"[STOP] {msg.strip()}", flush=True)
            self._write_stop_failed_record(run_i, msg)
            return False

        try:
            pid = int(pid_file.read_text().strip())
        except Exception as exc:
            msg = f"Could not read PID file {pid_file}: {exc}"
            print(f"[STOP] {msg}", flush=True)
            self._write_stop_failed_record(run_i, msg)
            return False

        cmdline = self._pid_cmdline(pid)
        if cmdline and "predictStructureQvary" not in cmdline:
            msg = (
                f"Refusing to kill PID {pid} for run {run_i}: command line does not look like predictStructureQvary.\n"
                f"cmdline={cmdline}\n"
            )
            print(f"[STOP] {msg.strip()}", flush=True)
            self._write_stop_failed_record(run_i, msg)
            return False

        if not self._pid_is_alive(pid):
            msg = (
                f"Run {run_i} already finished before termination signal.\n"
                f"Trigger: {label}, FoXS chi2={chi2:.8g}, qualifying_count={count}\n"
            )
            print(f"[STOP] {msg.strip()}", flush=True)
            self._write_stop_record(run_i, msg)
            return True

        try:
            print(f"[STOP] Terminating run {run_i} PID={pid}: {label}, FoXS chi2={chi2:.6g}", flush=True)
            os.kill(pid, signal.SIGTERM)
            time.sleep(1.0)
            if self._pid_is_alive(pid):
                print(f"[STOP] PID={pid} still alive; sending SIGKILL", flush=True)
                os.kill(pid, signal.SIGKILL)
        except ProcessLookupError:
            # Finished between checks. This is still a successful stop decision.
            pass
        except Exception as exc:
            msg = f"Failed to terminate run {run_i} PID={pid}: {exc}"
            print(f"[STOP] {msg}", flush=True)
            self._write_stop_failed_record(run_i, msg)
            return False

        msg = (
            f"Stopped run {run_i}\n"
            f"PID: {pid}\n"
            f"Trigger: {label}\n"
            f"FoXS chi2: {chi2:.8g}\n"
            f"Threshold: {self.cfg.terminate_threshold}\n"
            f"Qualifying count: {count}\n"
            f"Required count: {self.cfg.terminate_confirmation_count}\n"
        )
        self._write_stop_record(run_i, msg)
        return True

    def _consider_foxs_termination(self, run_i: int, label: str, chi2: Optional[float]) -> None:
        if not self.cfg.terminate_on_foxs:
            return
        if chi2 is None:
            return
        if self.cfg.terminate_threshold is None:
            return

        try:
            chi2_f = float(chi2)
            threshold = float(self.cfg.terminate_threshold)
        except Exception:
            return

        if chi2_f > threshold:
            return

        with self._termination_lock:
            if self._run_is_stopped(run_i):
                return
            count = int(self._good_foxs_counts.get(run_i, 0)) + 1
            self._good_foxs_counts[run_i] = count
            need = max(1, int(self.cfg.terminate_confirmation_count))
            print(
                f"[STOP-CHECK] run={run_i} {label} FoXS chi2={chi2_f:.6g} <= {threshold:.6g} "
                f"({count}/{need})",
                flush=True,
            )
            if count >= need:
                self._terminate_predictor_for_run(run_i, label, chi2_f, count)

    def _is_candidate_dat(self, dat: Path) -> bool:
        # Initial structures are produced at run launch and should not be part of
        # delayed monitoring/backmapping by default.
        if "_initial_xyz.dat" in dat.name:
            return False
        if not (dat.name.startswith("mol") and dat.name.endswith("_xyz.dat")):
            return False
        if not dat.is_file():
            return False
        if dat.suffix.lower() != ".dat":
            return False
        if any(str(dat).endswith(sfx) for sfx in self.cfg.ignore_suffixes):
            return False
        return parse_prediction_key(dat) is not None

    def _discover_groups(self):
        groups = {}
        for dat in sorted(self.cfg.watch_dir.glob(self.cfg.dat_glob)):
            if not self._is_candidate_dat(dat):
                continue
            parsed = parse_prediction_key(dat)
            if parsed is None:
                continue
            run_i, sub_i, tag = parsed
            if self._run_is_stopped(run_i):
                continue
            group_key = (run_i, tag)
            groups.setdefault(group_key, []).append((sub_i, dat))
        return groups

    def _group_signature(self, files: list[Path]):
        return tuple((str(f), f.stat().st_mtime if f.exists() else None) for f in files)

    def _group_is_complete(self, sub_files: list[tuple[int, Path]]) -> bool:
        if self.cfg.no_structures <= 1:
            return len(sub_files) >= 1
        return len({sub for sub, _ in sub_files}) >= self.cfg.no_structures

    def _group_is_ready(self, files: list[Path]) -> bool:
        for f in files:
            if not f.exists() or not f.is_file():
                return False
            if f.stat().st_size == 0:
                print(f"[WARN] Empty file, will retry later: {f.name}", flush=True)
                return False
            ok = wait_until_stable(
                f,
                self.cfg.stable_for,
                self.cfg.stable_poll,
                self.cfg.stable_timeout,
            )
            if not ok:
                print(f"[WARN] Never stabilized: {f.name}", flush=True)
                return False
        return True

    def _run_backmap_group_limited(self, group_key, sub_files: list[tuple[int, Path]], group_signature):
        if self._stop.is_set():
            return

        with self._sem:
            try:
                if self._stop.is_set():
                    return

                sub_indices = [sub for sub, _ in sub_files]
                sub_index_base = infer_sub_index_base(sub_indices, self.cfg.no_structures)
                files = [dat for _, dat in sorted(sub_files, key=lambda x: x[0])]
                names = ", ".join(f.name for f in files)
                print(f"[BACKMAP-GROUP] starting group={group_key}: {names}", flush=True)

                ok_all = True
                # For mixtures, FoXS must be evaluated at group level after all
                # components are backmapped. For ordinary single-structure runs,
                # keep the existing per-file FoXS path in backmap_cli.py.
                per_file_foxs = self.cfg.do_foxs and self.cfg.no_structures <= 1
                for dat_file in files:
                    if self._stop.is_set():
                        return
                    rc = run_backmap(self.cfg, dat_file, sub_index_base=sub_index_base,
                                     do_foxs_override=per_file_foxs)
                    if rc != 0:
                        ok_all = False
                    elif per_file_foxs:
                        chi2 = read_single_foxs_score_for_dat(self.cfg, dat_file)
                        if chi2 is not None:
                            self._consider_foxs_termination(group_key[0], dat_file.stem, chi2)

                if ok_all and self.cfg.do_foxs and self.cfg.no_structures > 1:
                    rc, chi2 = run_foxs_mixture_group(self.cfg, group_key, files)
                    if rc != 0:
                        ok_all = False
                    elif chi2 is not None:
                        self._consider_foxs_termination(group_key[0], _group_label(group_key), chi2)

                if ok_all:
                    self._completed_groups.add(group_key)
                    self._processed_groups[group_key] = group_signature
            finally:
                self._inflight_groups.discard(group_key)

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
            groups = self._discover_groups()

            for group_key in sorted(groups.keys(), key=lambda x: (x[0], x[1][0], x[1][1])):
                if group_key in self._inflight_groups:
                    continue
                if (not self.cfg.overwrite) and (group_key in self._completed_groups):
                    continue

                sub_files = sorted(groups[group_key], key=lambda x: x[0])

                if not self._group_is_complete(sub_files):
                    continue

                # If more files are present than expected, take the first no_structures
                # by sub-index. This avoids one malformed duplicate blocking the group.
                if self.cfg.no_structures > 0 and len(sub_files) > self.cfg.no_structures:
                    sub_files = sub_files[: self.cfg.no_structures]

                files = [dat for _, dat in sub_files]
                group_signature = self._group_signature(files)

                if self._processed_groups.get(group_key) == group_signature:
                    continue

                # Defer/ignore old groups. If all files in a group were written before
                # activation, mark that exact group signature as processed so it does
                # not backlog after activation.
                mtimes = [f.stat().st_mtime for f in files if f.exists()]
                if mtimes and max(mtimes) < self._activation_time:
                    self._processed_groups[group_key] = group_signature
                    continue

                if not self._group_is_ready(files):
                    continue

                print(
                    f"[WATCHER] running backmap group run={group_key[0]}, tag={group_key[1]}, n={len(files)}",
                    flush=True,
                )
                self._inflight_groups.add(group_key)
                threading.Thread(
                    target=self._run_backmap_group_limited,
                    args=(group_key, sub_files, group_signature),
                    daemon=True,
                ).start()

            time.sleep(self.cfg.poll_interval)


_WATCHER = None


def start_watcher(watch_dir, scenario_root, backmap_script,
                  out_dir=None, poll_interval=0.5, overwrite=False,
                  max_backmap=1, defer_backmap_seconds=0.0,
                  no_structures=1,
                  do_foxs=False, foxs_py=None, saxs_dat=None, max_q=None,
                  backend="modeller", cg2all_exec=None, disulfide_file=None,
                  terminate_on_foxs=False, terminate_threshold=None,
                  terminate_confirmation_count=1):
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
        no_structures=int(no_structures),
        do_foxs=do_foxs,
        foxs_py=str(foxs_py) if foxs_py else None,
        saxs_dat=Path(saxs_dat).resolve() if saxs_dat else None,
        max_q=max_q,
        backend=backend,
        cg2all_exec=cg2all_exec,
        disulfide_file=Path(disulfide_file).resolve() if disulfide_file else None,
        terminate_on_foxs=bool(terminate_on_foxs),
        terminate_threshold=float(terminate_threshold) if terminate_threshold is not None else None,
        terminate_confirmation_count=max(1, int(terminate_confirmation_count)),
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
    ap.add_argument("--defer-backmap-seconds", type=float, default=0.0)
    ap.add_argument("--no-structures", type=int, default=1)

    ap.add_argument("--do-foxs", action="store_true")
    ap.add_argument("--foxs-py", default=None)
    ap.add_argument("--saxs", default=None)
    ap.add_argument("--max-q", type=float, default=None)

    ap.add_argument("--backend", choices=["modeller", "cg2all"], default="modeller")
    ap.add_argument("--cg2all-exec", default=None)
    ap.add_argument("--disulfide-file", default=None)

    ap.add_argument("--terminate-on-foxs", action="store_true",
                    help="Opt-in mode: stop individual predictor runs once FoXS chi^2 is good enough")
    ap.add_argument("--terminate-threshold", type=float, default=2.5,
                    help="FoXS chi^2 threshold used with --terminate-on-foxs")
    ap.add_argument("--terminate-confirmation-count", type=int, default=1,
                    help="Number of qualifying FoXS scores required before stopping a run")

    args = ap.parse_args()

    if args.no_structures < 1:
        ap.error("--no-structures must be >= 1")
    if args.do_foxs and (args.saxs is None or args.max_q is None or args.foxs_py is None):
        ap.error("--do-foxs requires --foxs-py, --saxs, and --max-q")
    if args.backend == "cg2all" and args.cg2all_exec is None:
        ap.error("--backend cg2all requires --cg2all-exec")
    if args.terminate_confirmation_count < 1:
        ap.error("--terminate-confirmation-count must be >= 1")
    if args.terminate_threshold <= 0:
        ap.error("--terminate-threshold must be > 0")
    if args.terminate_on_foxs and not args.do_foxs:
        ap.error("--terminate-on-foxs requires --do-foxs")

    cfg = WatchConfig(
        watch_dir=Path(args.watch_dir).resolve(),
        scenario_root=Path(args.scenario_root).resolve(),
        backmap_script=Path(args.backmap_script).resolve(),
        poll_interval=args.poll,
        overwrite=args.overwrite,
        max_backmap=args.max_backmap,
        defer_backmap_seconds=args.defer_backmap_seconds,
        no_structures=args.no_structures,
        do_foxs=args.do_foxs,
        foxs_py=str(args.foxs_py) if args.foxs_py else None,
        saxs_dat=Path(args.saxs).resolve() if args.saxs else None,
        max_q=args.max_q,
        backend=args.backend,
        cg2all_exec=args.cg2all_exec,
        disulfide_file=Path(args.disulfide_file).resolve() if args.disulfide_file else None,
        terminate_on_foxs=args.terminate_on_foxs,
        terminate_threshold=float(args.terminate_threshold) if args.terminate_on_foxs else None,
        terminate_confirmation_count=max(1, int(args.terminate_confirmation_count)),
    )

    print("[WATCHER] started", flush=True)
    print("watch_dir      :", cfg.watch_dir, flush=True)
    print("scenario_root  :", cfg.scenario_root, flush=True)
    print("backmap_script :", cfg.backmap_script, flush=True)
    print("max_backmap    :", cfg.max_backmap, flush=True)
    print("no_structures  :", cfg.no_structures, flush=True)
    if cfg.terminate_on_foxs:
        print("mode           : terminate-on-FoXS", flush=True)
        print("term_threshold :", cfg.terminate_threshold, flush=True)
        print("term_confirm   :", cfg.terminate_confirmation_count, flush=True)
    else:
        print("mode           : maximal exploration", flush=True)
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
