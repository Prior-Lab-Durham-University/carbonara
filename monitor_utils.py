from pathlib import Path
import re

def parse_foxs_results_file(path: Path):
    vals = []
    for line in path.read_text().splitlines():
        parts = line.strip().split()
        if len(parts) < 2:
            continue
        if parts[1] == "ERROR":
            vals.append(None)
        else:
            try:
                vals.append(float(parts[1]))
            except:
                vals.append(None)
    return vals

_MIX_CHI_RE = re.compile(r"\bchi2=([0-9.eE+-]+)")

def read_foxs_scores(fitdata_dir, mixture_n=1):
    fitdata_dir = Path(fitdata_dir)

    rows = []

    if mixture_n <= 1:
        for f in fitdata_dir.glob("allAtomRun*/foxs_results.txt"):
            for line in f.read_text().splitlines():
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        rows.append({
                            "source": str(f),
                            "model": parts[0],
                            "chi2": float(parts[1]),
                        })
                    except ValueError:
                        pass
    else:
        for f in fitdata_dir.glob("allAtomRun*/foxs_mixture_results.txt"):
            for line in f.read_text().splitlines():
                m = _MIX_CHI_RE.search(line)
                if m:
                    rows.append({
                        "source": str(f),
                        "model": line.split()[0],
                        "chi2": float(m.group(1)),
                        "line": line,
                    })

    return rows

def sweep_quality(fitdata_dir: Path, threshold: float, mixture_n: int = 1):
    """
    Count good FoXS-scored predictions.

    For ordinary single-structure runs, read allAtomRun*/foxs_results.txt.
    For mixture/ensemble runs, read allAtomRun*/foxs_mixture_results.txt,
    where one line corresponds to one complete mixture state.
    """
    rows = read_foxs_scores(fitdata_dir, mixture_n=mixture_n)

    vals = [r.get("chi2") for r in rows]
    valid = [v for v in vals if v is not None]
    good = [v for v in valid if v <= threshold]

    return {
        "good": len(good),
        "total": len(rows),
        "errors": sum(v is None for v in vals),
        "best": min(valid) if valid else None,
    }
