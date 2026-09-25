# Week 1 spoilers — read after doing the loop yourself

This is Claude's experimental loop for `week-01-pi-error-vs-work.py`. It was
written on 2026-09-25, before the Week 1 reading. It was moved out of the
script so you can do your own loop first and then compare.

## Code

The error-shape analysis, moved out of `week-01-pi-error-vs-work.py`.
It prints its table; the old `output/week-01-leibniz-tail.csv` is no longer kept.

```python
# Run from experiments/. Reuses leibniz_partial_sums from the week-01 script.
import importlib.util
from pathlib import Path
from mpmath import mp, mpf

spec = importlib.util.spec_from_file_location("w01", Path("week-01-pi-error-vs-work.py"))
w01 = importlib.util.module_from_spec(spec); spec.loader.exec_module(w01)
leibniz_partial_sums = w01.leibniz_partial_sums
DPS = 60


def leibniz_tail(ns):
    """Scaled signed error and its successive corrections.

    s1 = n * (-1)^n * (pi - L(n))
    s3 = n^3 * ((-1)^n (pi - L(n)) - 1/n)
    s5 = n^5 * ((-1)^n (pi - L(n)) - 1/n + 1/(4n^3))
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


def main():
    mp.dps = DPS
    ns = [10, 11, 100, 101, 1000, 1001, 10000, 10001]
    tail = leibniz_tail(ns)

    # Verify at 2x precision: the scaled constants must not move.
    mp.dps = 2 * DPS
    tail2 = leibniz_tail(ns)
    mp.dps = DPS
    print(f"--- Leibniz error shape (dps={DPS}; max drift vs dps={2 * DPS}) ---")
    print("  n                      s1                 s3                 s5")
    worst = mpf(0)
    for (n, s1, s3, s5), (_, t1, t3, t5) in zip(tail, tail2):
        worst = max(worst, abs(s1 - t1), abs(s3 - t3), abs(s5 - t5))
        print(f"  {n:<6} {mp.nstr(s1, 15):>18} {mp.nstr(s3, 15):>18} {mp.nstr(s5, 15):>18}")
    print(f"  max |dps vs 2*dps| drift: {mp.nstr(worst, 3)}")


main()
```

## Compute

See `output/`. Leibniz needs ~10^5 terms for 5 digits. Machin gets ~1.45
digits per term, ~58 digits at 40 terms.

## Observe

The Leibniz error is not just "small, about 1/n"; it has an exact
structure. As n grows:

- n·(−1)^n·(π − L(n)) → 1
- the next correction → −1/4
- the one after that → 5/16

These are stable at 2x precision (60 vs 120 digits).

## Conjecture

(−1)^n (π − L(n)) = Σ_m 2·E_{2m} / (2n)^(2m+1), where
E_0, E_2, E_4, E_6, … = 1, −1, 5, −61, … are the Euler numbers.

## Attack

Predict the next coefficient: 2·(−61)/2^7 = −61/64 = −0.953125. Measure
s7 = n^7 · (remainder after the 5/16 term) at n = 10^4 and 10^5, at 120
then 240 digits. Any drift away from −61/64 kills the conjecture.

Result (scratch check, n = 10^4 only): s7 = −0.9531249459 at both 120 and
240 digits, so it survived. n = 10^5 has not been tested yet.

This is a known theorem (Borwein, Borwein & Dilcher, 1989), so it is a
rediscovery, not a new result.
