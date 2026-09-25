"""Week 1 experiment (part 2): pi approximations, error versus work.

Extends week-01-pi-approximation.py (imported, not modified) with:
  - a high-precision Leibniz partial sum, so the error can be studied
    far below float resolution;
  - Machin's formula, pi = 16 arctan(1/5) - 4 arctan(1/239), as a fast
    series to compare against Leibniz on the same error-vs-work axes.

Work is measured as number of series terms summed.

Outputs (in experiments/output/):
  week-01-error-vs-work.csv   method, terms, error at working precision
  week-01-error-vs-work.png   log-log plot of error vs terms

from_the_book() is still a stub in the baseline file. Once the book's
digits-of-pi method is read, implement it in a new week-01 file and add it
to the comparison here.
"""
import csv
import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpmath import mp, mpf

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
DPS = 60  # working precision (decimal digits); verification reruns at 2*DPS


def _load_baseline():
    spec = importlib.util.spec_from_file_location(
        "week01_baseline", HERE / "week-01-pi-approximation.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


baseline = _load_baseline()


# --- Extensions -----------------------------------------------------------

def leibniz_partial_sums(checkpoints):
    """Yield (n, L(n)) at each n in checkpoints, summing once at mp.dps."""
    checkpoints = sorted(checkpoints)
    total = mpf(0)
    k = 0
    for n in checkpoints:
        while k < n:
            term = mpf(1) / (2 * k + 1)
            total += term if k % 2 == 0 else -term
            k += 1
        yield n, 4 * total


def arctan_inv(x, n_terms):
    """arctan(1/x) truncated to n_terms of its Taylor series."""
    total = mpf(0)
    power = mpf(1) / x
    x2 = x * x
    for k in range(n_terms):
        term = power / (2 * k + 1)
        total += term if k % 2 == 0 else -term
        power /= x2
    return total


def machin_pi(n_terms):
    """Machin: pi = 16 arctan(1/5) - 4 arctan(1/239), n_terms per arctan."""
    return 16 * arctan_inv(5, n_terms) - 4 * arctan_inv(239, n_terms)


def error_vs_work():
    """Rows of (method, terms, |pi - approx|) at the current mp.dps."""
    rows = []
    leib_n = [10**j for j in range(1, 6)]
    for n, approx in leibniz_partial_sums(leib_n):
        rows.append(("leibniz", n, abs(mp.pi - approx)))
    for n in range(1, 41):
        rows.append(("machin", n, abs(mp.pi - machin_pi(n))))
    return rows


# --- Driver ---------------------------------------------------------------

def write_csv(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow([mp.nstr(v, 20) if isinstance(v, type(mpf(0))) else v for v in r])


def plot(rows, path):
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    for method, color in (("leibniz", "#1f6feb"), ("machin", "#d9480f")):
        pts = [(n, float(err)) for m, n, err in rows if m == method and err > 0]
        ax.loglog(*zip(*pts), marker="o", ms=3, lw=1.2, color=color, label=method)
    ax.set_xlabel("terms summed (work)")
    ax.set_ylabel("|pi - approximation|")
    ax.set_title(f"Error vs work (mp.dps = {DPS})")
    ax.grid(True, which="major", alpha=0.3)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)

    # Baseline, unchanged (float arithmetic).
    baseline.report("Leibniz series (baseline, float)", baseline.leibniz_pi,
                    [10, 100, 1000, 10000])

    mp.dps = DPS
    rows = error_vs_work()
    write_csv(OUT / "week-01-error-vs-work.csv", ["method", "terms", "abs_error"], rows)
    plot(rows, OUT / "week-01-error-vs-work.png")

    print(f"--- Machin, mp.dps={DPS} ---")
    for m, n, err in rows:
        if m == "machin" and n in (1, 2, 5, 10, 20, 40):
            print(f"  terms={n:>3}  error={mp.nstr(err, 3)}")
    print(f"Wrote CSV and plot to {OUT}")


if __name__ == "__main__":
    main()
