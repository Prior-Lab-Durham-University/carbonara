import os
import subprocess
import signal
import threading
import time
from pathlib import Path
from IPython.display import display, Markdown

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
        self._display_handle = None

    # ------------------------
    # Run control
    # ------------------------
    def start(self):
      if self.proc is not None:
          print("Already running.")
          return

      self.proc = subprocess.Popen(
          ["bash", self.script, self.foxs_cmd],
          preexec_fn=os.setsid
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
        except:
            try:
                os.killpg(pgid, signal.SIGTERM)
                self.proc.wait(timeout=5)
            except:
                os.killpg(pgid, signal.SIGKILL)
                self.proc.wait()

        self.proc = None
        print("✅ Stopped")

    # ------------------------
    # Monitoring
    # ------------------------
    def start_monitor(self, threshold=2.5, every_s=10):
        self._stop_event.clear()

        # Create a persistent display handle once
        if self._display_handle is None:
            self._display_handle = display(Markdown("Starting monitor..."), display_id=True)

        def loop():
            while not self._stop_event.is_set():
                stats = sweep_quality(self.fitdata, threshold)

                lines = []
                lines.append("## 🔬 Carbonara Live Monitor")
                lines.append("")
                lines.append(f"**Threshold (χ²):** {threshold}")
                lines.append(f"**Good models:** {stats['good']} / {stats['total']}")
                lines.append(f"**Errors:** {stats['errors']}")
                if stats["best"] is not None:
                    lines.append(f"**Best χ²:** {stats['best']:.4f}")

                text = "\n".join(lines)

        print("📊 Monitoring started")

    def stop_monitor(self):
        self._stop_event.set()
        print("📊 Monitoring stopped")
