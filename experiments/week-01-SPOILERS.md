# Week 1 spoilers — read after doing the loop yourself

This is Claude's experimental loop for `week-01-pi-error-vs-work.py`. It was
written on 2026-09-25, before the Week 1 reading. It was moved out of the
script so you can do your own loop first and then compare.

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
