import os
import sys
import re
import subprocess
import signal
import threading
import time
from pathlib import Path
from IPython.display import display, Markdown, update_display

from monitor_utils import sweep_quality


def check_modeller_available() -> bool:
    """
    Return True only if MODELLER is importable and can initialise an environ().

    MODELLER is an optional external dependency because it requires a separate
    licence. It should therefore be checked only when the user requests the
    Modeller backmapping backend.
    """
    try:
        from modeller import environ  # type: ignore
        _ = environ()
        return True
    except Exception:
        return False


def modeller_error_message() -> str:
    return (
        "MODELLER backmapping was requested, but MODELLER is not installed "
        "or is not licensed/configured in this Python environment.\n\n"
        "MODELLER is optional and is not installed by setupPython.sh because "
        "it requires a separate licence. Either install/configure MODELLER on "
        "this machine, or choose a non-Modeller backmapping method such as "
        "CG2ALL if available.\n\n"
        "Check manually with:\n"
        "    python -c \"from modeller import environ; env=environ(); print('MODELLER OK')\""
    )

def check_biopython_available() -> bool:
    """Return True if Biopython/Bio.PDB is importable in the current kernel Python."""
    try:
        import Bio  # noqa: F401
        from Bio.PDB import PDBParser, PDBIO  # noqa: F401
        return True
    except Exception:
        return False


def biopython_error_message() -> str:
    return (
        "Biopython is required by the Carbonara backmapping code, but Bio/Bio.PDB "
        "is not importable in this Python environment.\n\n"
        f"Current Python executable: {sys.executable}\n\n"
        "Install it in this environment with:\n"
        "    python -m pip install biopython"
    )



def detect_backmap_method_from_run_script(script_path: str | Path) -> str | None:
    """
    Try to detect a backmapping method from a RunMe_*.sh script.

    This deliberately supports several common spellings, since notebook/run
    scripts often evolve:
      BACKMAP_METHOD=modeller
      method=modeller
      --method modeller
      --backmap-method modeller
      method='modeller'

    Returns 'modeller', 'cg2all', or None if no method is identifiable.
    """
    script_path = Path(script_path)
    if not script_path.exists():
        return None

    try:
        text = script_path.read_text(errors="ignore")
    except Exception:
        return None

    patterns = [
        r"(?im)^\s*(?:BACKMAP_METHOD|BACKMAPPING_METHOD|METHOD|method)\s*=\s*[\"']?(modeller|cg2all)[\"']?\b",
        r"(?i)--(?:backmap-method|method)\s+[\"']?(modeller|cg2all)[\"']?\b",
        r"(?i)\bmethod\s*=\s*[\"']?(modeller|cg2all)[\"']?\b",
    ]

    for pat in patterns:
        match = re.search(pat, text)
        if match:
            return match.group(1).lower()

    return None


def detect_no_structures_from_run_script(script_path: str | Path) -> int:
    """
    Read the generated RunMe_*.sh script and extract the Carbonara
    noStructures value. This is the authoritative mixture/component count
    used both by predictStructureQvary and the watcher.
    """
    script_path = Path(script_path)
    if not script_path.exists():
        return 1

    try:
        text = script_path.read_text(errors="ignore")
    except Exception:
        return 1

    # Standard generated form:
    #   noStructures=2
    m = re.search(r"(?im)^\s*noStructures\s*=\s*[\"']?(\d+)[\"']?\s*(?:#.*)?$", text)
    if m:
        try:
            return max(1, int(m.group(1)))
        except Exception:
            return 1

    # Fallback: watcher arg form:
    #   --no-structures "$noStructures"  or --no-structures 2
    m = re.search(r"(?i)--no-structures\s+[\"']?(\d+)[\"']?", text)
    if m:
        try:
            return max(1, int(m.group(1)))
        except Exception:
            return 1

    return 1


def _read_run_script_text(script_path: str | Path) -> str:
    script_path = Path(script_path)
    if not script_path.exists():
        return ""
    try:
        return script_path.read_text(errors="ignore")
    except Exception:
        return ""


def _detect_bool_assignment(text: str, name: str, default: bool = False) -> bool:
    m = re.search(rf"(?im)^\s*{re.escape(name)}\s*=\s*[\"']?(True|False|true|false|1|0|yes|no)[\"']?", text)
    if not m:
        return default
    return m.group(1).strip().lower() in {"true", "1", "yes"}


def _detect_numeric_assignment(text: str, name: str, default):
    m = re.search(rf"(?im)^\s*{re.escape(name)}\s*=\s*[\"']?([0-9.eE+-]+)[\"']?", text)
    if not m:
        return default
    try:
        if isinstance(default, int):
            return int(float(m.group(1)))
        return float(m.group(1))
    except Exception:
        return default


def detect_run_mode_from_run_script(script_path: str | Path) -> dict:
    """
    Detect whether a generated RunMe_*.sh is in maximal-exploration mode or
    terminate-on-FoXS mode. Missing/old scripts are treated as maximal exploration.
    """
    text = _read_run_script_text(script_path)
    terminate_on_foxs = _detect_bool_assignment(text, "TERMINATE_ON_FOXS", default=False)

    # Fallback for scripts that pass the watcher arg directly.
    if not terminate_on_foxs and re.search(r"(?i)--terminate-on-foxs\b", text):
        terminate_on_foxs = True

    threshold = _detect_numeric_assignment(text, "TERMINATE_FOXS_THRESHOLD", 2.5)
    confirmation = _detect_numeric_assignment(text, "TERMINATE_CONFIRMATION_COUNT", 1)
    confirmation = max(1, int(confirmation))

    return {
        "terminate_on_foxs": bool(terminate_on_foxs),
        "mode": "terminate-on-FoXS" if terminate_on_foxs else "maximal exploration",
        "terminate_threshold": float(threshold),
        "terminate_confirmation_count": confirmation,
    }


class CarbonaraRunner:

    def __init__(
        self,
        project_name,
        base_dir="carbonara_runs",
        foxs_cmd="pyfoxs",
        backmap_method="auto",
        require_modeller_check=True,
    ):
        self.project_name = project_name
        self.base_dir = Path(base_dir)
        self.script = f"RunMe_{project_name}.sh"
        self.fitdata = self.base_dir / project_name / "fitdata"
        self.foxs_cmd = foxs_cmd
        self.backmap_method = backmap_method
        self.require_modeller_check = require_modeller_check

        self.proc = None
        self._monitor_thread = None
        self._stop_event = threading.Event()
        self._display_id = None
        self._monitor_state = {
            "good": 0,
            "total": 0,
            "errors": 0,
            "best": None,
            "threshold": 2.5,
            "last_update": None,
            "start_time": None,
            "defer_backmap_seconds": 600,
            "mixture_n": 1,
            "mode": "maximal exploration",
            "terminate_on_foxs": False,
            "terminate_threshold": None,
            "terminate_confirmation_count": 1,
            "stopped_runs": 0,
        }

    # ------------------------
    # Run control
    # ------------------------
    def _resolve_backmap_method(self) -> str | None:
        """
        Resolve requested backmapping method.

        If self.backmap_method is 'auto', inspect the generated RunMe script.
        If no method can be detected, return None and do not block startup.
        """
        if self.backmap_method is None:
            return None

        method = str(self.backmap_method).strip().lower()
        if method in {"", "none", "false", "off", "no"}:
            return None

        if method == "auto":
            return detect_backmap_method_from_run_script(self.script)

        return method

    def _preflight_checks(self):
        script_path = Path(self.script)
        if not script_path.exists():
            raise FileNotFoundError(
                f"Run script not found: {script_path}\n"
                f"Expected to find RunMe_{self.project_name}.sh in the current directory."
            )

        method = self._resolve_backmap_method()

        # The generated shell script receives this kernel Python via PYTHON_EXE,
        # so check core Python-side backmapping dependencies here before launching.
        if not check_biopython_available():
            raise RuntimeError(biopython_error_message())
        print("✅ Biopython/Bio.PDB check passed")

        if method == "modeller" and self.require_modeller_check:
            if not check_modeller_available():
                raise RuntimeError(modeller_error_message())
            print("✅ MODELLER check passed")
        elif method == "cg2all":
            print("ℹ️ Backmapping method detected: CG2ALL")
        elif method is None:
            print("ℹ️ No backmapping method detected in run script; skipping MODELLER preflight check")
        else:
            print(f"ℹ️ Backmapping method set to {method!r}; no MODELLER preflight needed")

    def start(self):
        if self.proc is not None:
            print("Already running.")
            return

        self._preflight_checks()

        env = os.environ.copy()
        env["PYTHON_EXE"] = sys.executable

        # Put the kernel Python directory first so any plain `python` calls inside
        # child shell scripts also resolve to the same environment. The generated
        # RunMe script should use PYTHON_EXE explicitly, but this makes older
        # generated scripts safer too.
        python_bin = str(Path(sys.executable).resolve().parent)
        env["PATH"] = python_bin + os.pathsep + env.get("PATH", "")

        self.proc = subprocess.Popen(
            ["bash", self.script, self.foxs_cmd],
            preexec_fn=os.setsid,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            env=env,
        )

        print(f"🚀 Carbonara started (PID={self.proc.pid})")
        print(f"🐍 Watcher/backmapping Python: {sys.executable}")

    def stop(self):
        if self.proc is None:
            print("Nothing running.")
            return

        pgid = os.getpgid(self.proc.pid)
        print("🛑 Stopping Carbonara...")

        try:
            os.killpg(pgid, signal.SIGINT)
            self.proc.wait(timeout=5)
        except Exception:
            try:
                os.killpg(pgid, signal.SIGTERM)
                self.proc.wait(timeout=5)
            except Exception:
                os.killpg(pgid, signal.SIGKILL)
                self.proc.wait()

        self.proc = None
        print("✅ Stopped")

    # ------------------------
    # Monitoring
    # ------------------------
    def _render_monitor_text(self):
        stats = self._monitor_state
        lines = []
        lines.append("## 🔬 Carbonara Live Monitor")
        lines.append("")

        if stats["start_time"] is not None:
            elapsed = time.time() - stats["start_time"]
            remaining = stats["defer_backmap_seconds"] - elapsed
            if remaining > 0:
                lines.append(f"**Backmapping starts in:** {int(remaining)} s")
                lines.append("**Status:** Waiting period active. Early predictions are being ignored.")
            else:
                lines.append("**Status:** Backmapping window active. New eligible predictions should now be picked up.")

        lines.append(f"**Monitor threshold (χ²):** {stats['threshold']}")
        lines.append(f"**Mixture components:** {stats.get('mixture_n', 1)}")
        lines.append(f"**Run mode:** {stats.get('mode', 'maximal exploration')}")
        if stats.get("terminate_on_foxs"):
            lines.append(
                f"**Termination rule:** FoXS χ² ≤ {stats.get('terminate_threshold')} "
                f"for {stats.get('terminate_confirmation_count', 1)} qualifying score(s)"
            )
            lines.append(f"**Early-stopped runs:** {stats.get('stopped_runs', 0)}")
        lines.append(f"**Good models:** {stats['good']} / {stats['total']}")
        lines.append(f"**Errors:** {stats['errors']}")
        if stats["best"] is not None:
            lines.append(f"**Best χ²:** {stats['best']:.4f}")
        if stats["last_update"] is not None:
            lines.append(f"**Last sweep:** {stats['last_update']}")

        return "\n".join(lines)

    def start_monitor(self, threshold=2.5, every_s=10, defer_backmap_seconds=600):
        if self._monitor_thread is not None and self._monitor_thread.is_alive():
            print("Monitoring already running.")
            return

        self._stop_event.clear()
        self._monitor_state["threshold"] = threshold
        self._monitor_state["start_time"] = time.time()
        self._monitor_state["defer_backmap_seconds"] = defer_backmap_seconds
        mixture_n = detect_no_structures_from_run_script(self.script)
        self._monitor_state["mixture_n"] = mixture_n
        mode_info = detect_run_mode_from_run_script(self.script)
        self._monitor_state.update(mode_info)
        self._monitor_state["stopped_runs"] = len(list(self.fitdata.glob("run*.stopped")))

        initial = Markdown(self._render_monitor_text())
        handle = display(initial, display_id=True)
        self._display_id = handle.display_id

        def loop():
            while not self._stop_event.is_set():
                try:
                    stats = sweep_quality(self.fitdata, threshold, mixture_n=mixture_n)
                    self._monitor_state.update(stats)
                except Exception:
                    self._monitor_state["errors"] += 1

                self._monitor_state["stopped_runs"] = len(list(self.fitdata.glob("run*.stopped")))
                self._monitor_state["last_update"] = time.strftime("%H:%M:%S")
                update_display(Markdown(self._render_monitor_text()), display_id=self._display_id)
                time.sleep(every_s)

        self._monitor_thread = threading.Thread(target=loop, daemon=True)
        self._monitor_thread.start()
        print(f"📊 Monitoring started (mixture_n={mixture_n}, mode={mode_info['mode']})")

    def stop_monitor(self):
        self._stop_event.set()
        print("📊 Monitoring stopped")
