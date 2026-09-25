"""Week 1 experiment (part 2): pi approximations, error versus work.

Extends week-01-pi-approximation.py (imported, not modified) with:
  - a high-precision Leibniz partial sum, so the error can be studied
    far below float resolution;
  - Machin's formula, pi = 16 arctan(1/5) - 4 arctan(1/239), as a fast
    series to compare against Leibniz on the same error-vs-work axes;
  - a look at the *shape* of the Leibniz error, not just its size.

Work is measured as number of series terms summed.

Outputs (in experiments/output/):
  week-01-error-vs-work.csv   method, terms, error at working precision
  week-01-leibniz-tail.csv    n, scaled error n*(-1)^n*(pi - L(n)) and corrections
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


def leibniz_tail(ns):
    """Scaled signed error and its successive corrections.

    s1 = n * (-1)^n * (pi - L(n))              -> 1 ?
    s3 = n^3 * ((-1)^n (pi - L(n)) - 1/n)      -> -1/4 ?
    s5 = n^5 * ((-1)^n (pi - L(n)) - 1/n + 1/(4n^3))  -> 5/16 ?
    """
    rows = []
    for n, approx in leibniz_partial_sums(ns):
        e = (-1) ** n * (mp.pi - approx)
        n_ = mpf(n)
        s1 = n_ * e
        s3 = n_**3 * (e - 1 / n_)
        s5 = n_**5 * (e - 1 / n_ + 1 / (4 * n_**3))
        rows.append((n, s1, s3, s5))
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
    print()

    ns = [10, 11, 100, 101, 1000, 1001, 10000, 10001]
    tail = leibniz_tail(ns)
    write_csv(OUT / "week-01-leibniz-tail.csv", ["n", "s1", "s3", "s5"], tail)

    # Verify at 2x precision: the scaled constants must not move.
    mp.dps = 2 * DPS
    tail2 = leibniz_tail(ns)
    mp.dps = DPS
    print(f"--- Leibniz error shape (dps={DPS}; max drift vs dps={2 * DPS}) ---")
    print("  n        s1 (->1?)         s3 (->-1/4?)        s5 (->5/16?)")
    worst = mpf(0)
    for (n, s1, s3, s5), (_, t1, t3, t5) in zip(tail, tail2):
        worst = max(worst, abs(s1 - t1), abs(s3 - t3), abs(s5 - t5))
        print(f"  {n:<6} {mp.nstr(s1, 15):>18} {mp.nstr(s3, 15):>18} {mp.nstr(s5, 15):>18}")
    print(f"  max |dps vs 2*dps| drift: {mp.nstr(worst, 3)}")
    print(f"\nWrote CSVs and plot to {OUT}")


if __name__ == "__main__":
    main()

# --- Experimental loop (filled in from the run on 2026-09-25) -------------
#
# Compute:  see output/ -- Leibniz needs ~10^5 terms for 5 digits; Machin
#           gets ~1.45 digits per term, ~58 digits at 40 terms.
# Observe:  the Leibniz error is not just "small, about 1/n" -- it is
#           *exactly structured*: n*(-1)^n*(pi - L(n)) -> 1, the next
#           correction -> -1/4, the next -> 5/16, stable at 2x precision.
# Conjecture: (-1)^n (pi - L(n)) = sum_m E_{2m} * 2 / (2n)^(2m+1), with
#           E_0, E_2, E_4, ... = 1, -1, 5, -61, ... the Euler numbers.
# Attack:   predict the next coefficient (-61 * 2 / 2^7 = -61/64) and test
#           s7 = n^7 * (remainder after the 5/16 term) at n = 10^4 and
#           10^5 with dps = 120 then 240; any drift away from -61/64 kills it.
