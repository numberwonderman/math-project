"""Week 1 experiment: reproduce one pi approximation from the text.

Reading plan, Week 1 (Sep 26 - Oct 2, 2026):
  Experiment: reproduce one pi approximation from the text and record
  error versus work.

This script starts you with the Leibniz series so it runs on day one.
Your job: replace it with (or add alongside it) the approximation method
from the digits-of-pi section, then compare error versus work for each.
"""
import math


def leibniz_pi(n_terms):
    """pi/4 = 1 - 1/3 + 1/5 - 1/7 + ...  (slow, but simple and exact)."""
    total = 0.0
    for k in range(n_terms):
        total += ((-1) ** k) / (2 * k + 1)
    return 4 * total


def from_the_book(n_terms):
    """TODO: implement the approximation from the digits-of-pi section.

    Replace this body with the method you read about, then compare its
    error-versus-work against leibniz_pi above.
    """
    raise NotImplementedError("Fill this in after reading the digits-of-pi section")


def report(name, func, work_sizes):
    print(f"--- {name} ---")
    for n in work_sizes:
        approx = func(n)
        error = abs(approx - math.pi)
        print(f"  work={n:>8}  approx={approx:.10f}  error={error:.2e}")
    print()


def main():
    sizes = [10, 100, 1000, 10000]
    report("Leibniz series (starter)", leibniz_pi, sizes)
    # Uncomment once you implement from_the_book():
    # report("From the book", from_the_book, sizes)


if __name__ == "__main__":
    main()
