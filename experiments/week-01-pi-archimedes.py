"""Week 1 experiment: Archimedes' polygon method for pi, error versus work.

The Week 1 target: reproduce one pi approximation from the text and record
error versus work.

Work is measured as the number of side-doublings, starting from hexagons
(k doublings -> polygons with 6 * 2^k sides).

Outputs (in experiments/output/):
  week-01-archimedes.csv   doublings, sides, lower, upper, errors, log10 errors
  week-01-archimedes.png   log-scale error vs doublings, both bounds

Log: fill in week-01-pi-archimedes-log.md (from LOG-TEMPLATE.md).
"""
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpmath import log10, mp, mpf, sqrt

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
DPS = 60  # working precision (decimal digits); verification reruns at 2*DPS
DOUBLINGS = 30


def from_the_book(doublings):
    """Archimedes: bracket pi between inscribed and circumscribed polygons.

    CHECK AGAINST THE BOOK. This is the standard harmonic/geometric-mean
    form of the method, not copied from the text. If the book uses side
    lengths, angles, or a different recurrence, replace this body with its
    version.

    For a unit-diameter circle, let a = circumscribed perimeter and
    b = inscribed perimeter, starting from hexagons (a = 2*sqrt(3), b = 3).
    Doubling the number of sides gives
        a' = 2ab / (a + b)      (harmonic mean)
        b' = sqrt(a' * b)       (geometric mean)
    and b < pi < a at every step.

    Returns a list of (k, sides, lower, upper) for k = 0..doublings.
    """
    a = 2 * sqrt(3)
    b = mpf(3)
    rows = [(0, 6, b, a)]
    for k in range(1, doublings + 1):
        a = 2 * a * b / (a + b)
        b = sqrt(a * b)
        rows.append((k, 6 * 2**k, b, a))
    return rows


def with_errors(rows):
    return [(k, n, lo, hi, mp.pi - lo, hi - mp.pi) for k, n, lo, hi in rows]


def write_csv(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["doublings", "sides", "lower", "upper", "err_lower", "err_upper",
                    "log10_err_lower", "log10_err_upper"])
        for k, n, lo, hi, el, eu in rows:
            w.writerow([k, n] + [mp.nstr(v, 20) for v in (lo, hi, el, eu)]
                       + [mp.nstr(log10(v), 10) for v in (el, eu)])


def plot(rows, path):
    ks = [r[0] for r in rows]
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.semilogy(ks, [float(r[4]) for r in rows], marker="o", ms=3, lw=1.2,
                color="#1f6feb", label="pi - inscribed (lower)")
    ax.semilogy(ks, [float(r[5]) for r in rows], marker="s", ms=3, lw=1.2,
                color="#d9480f", label="circumscribed - pi (upper)")
    ax.set_xlabel("side doublings k  (6 * 2^k sides)")
    ax.set_ylabel("error")
    ax.set_title(f"Archimedes' polygons, error vs work (mp.dps = {DPS})")
    ax.grid(True, which="major", alpha=0.3)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)

    mp.dps = DPS
    rows = with_errors(from_the_book(DOUBLINGS))
    write_csv(OUT / "week-01-archimedes.csv", rows)
    plot(rows, OUT / "week-01-archimedes.png")

    # Sanity: the bounds must bracket pi at every step.
    assert all(lo < mp.pi < hi for _, _, lo, hi, _, _ in rows)

    # Verify at 2x precision: the errors must not move.
    mp.dps = 2 * DPS
    rows2 = with_errors(from_the_book(DOUBLINGS))
    drift = max(abs(r[i] - s[i]) for r, s in zip(rows, rows2) for i in (4, 5))
    mp.dps = DPS

    print(f"--- Archimedes, mp.dps={DPS} ---")
    print("  k      sides        lower bound            upper bound        err_lo     err_hi"
          "   log10_lo   log10_hi")
    for k, n, lo, hi, el, eu in rows:
        print(f"  {k:<3} {n:>12}  {mp.nstr(lo, 18):<20} {mp.nstr(hi, 18):<20} "
              f"{mp.nstr(el, 3):>9}  {mp.nstr(eu, 3):>9}  "
              f"{float(log10(el)):>9.3f}  {float(log10(eu)):>9.3f}")
    print(f"  max |dps vs 2*dps| drift in errors: {mp.nstr(drift, 3)}")
    print(f"Wrote CSV and plot to {OUT}")


if __name__ == "__main__":
    main()
