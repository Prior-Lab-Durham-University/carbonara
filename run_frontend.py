import os
import subprocess
import signal
import threading
import time
from pathlib import Path
from IPython.display import display, Markdown, update_display

from monitor_utils import sweep_quality


class CarbonaraRunner:

    def __init__(self, project_name, base_dir="carbonara_runs", foxs_cmd="pyfoxs"):
        self.project_name = project_name
        self.base_dir = Path(base_dir)
        self.script = f"RunMe_{project_name}.sh"
        self.fitdata = self.base_dir / project_name / "fitdata"
        self.foxs_cmd = foxs_cmd

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
        }

    # ------------------------
    # Run control
    # ------------------------
    def start(self):
        if self.proc is not None:
            print("Already running.")
            return

        self.proc = subprocess.Popen(
            ["bash", self.script, self.foxs_cmd],
            preexec_fn=os.setsid,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        print(f"🚀 Carbonara started (PID={self.proc.pid})")

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

        lines.append(f"**Threshold (χ²):** {stats['threshold']}")
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

        initial = Markdown(self._render_monitor_text())
        handle = display(initial, display_id=True)
        self._display_id = handle.display_id

        def loop():
            while not self._stop_event.is_set():
                try:
                    stats = sweep_quality(self.fitdata, threshold)
                    self._monitor_state.update(stats)
                except Exception:
                    self._monitor_state["errors"] += 1

                self._monitor_state["last_update"] = time.strftime("%H:%M:%S")
                update_display(Markdown(self._render_monitor_text()), display_id=self._display_id)
                time.sleep(every_s)

        self._monitor_thread = threading.Thread(target=loop, daemon=True)
        self._monitor_thread.start()
        print("📊 Monitoring started")

    def stop_monitor(self):
        self._stop_event.set()
        print("📊 Monitoring stopped")
