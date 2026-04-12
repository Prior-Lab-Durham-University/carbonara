from pathlib import Path

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


def sweep_quality(fitdata_dir: Path, threshold: float):
    run_dirs = sorted([p for p in fitdata_dir.glob("allAtomRun*") if p.is_dir()])

    total_good = total_total = total_err = 0
    best_global = None

    for rd in run_dirs:
        res = rd / "foxs_results.txt"
        if not res.exists():
            continue

        chis = parse_foxs_results_file(res)

        valid = [c for c in chis if c is not None]
        good = [c for c in valid if c <= threshold]

        total_good += len(good)
        total_total += len(chis)
        total_err += sum(c is None for c in chis)

        if valid:
            m = min(valid)
            best_global = m if best_global is None else min(best_global, m)

    return {
        "good": total_good,
        "total": total_total,
        "errors": total_err,
        "best": best_global,
    }
