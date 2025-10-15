#!/usr/bin/env python

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional

import numpy as np


@dataclass(frozen=True)
class RateMap:
    """Binned recombination map with precomputed cumulative intensity (Λ)."""

    length: int                 # chromosome length (bp)
    bins: int                   # number of bins
    dx: float                   # bp per bin
    rates: List[float]          # per-bp rate in each bin (len=bins)
    cum_intensity: List[float]  # Λ at bin edges (len=bins+1)
    s_max: float                # total expected crossovers (= Λ[-1])


def make_sine_rate_map(length: int, bins: int, mean_rate: float) -> RateMap:
    """Return a simple sine-parameterized RateMap over [0, length).

    Uses c(x) = mean_rate + (mean_rate/2) * sin(2π * 2 * x / length).
    Edit inside if you want amplitude/phase/cycles as arguments.
    """
    if length <= 0 or bins <= 0:
        raise ValueError("length and bins must be positive.")
    amplitude = mean_rate * 0.9
    if mean_rate <= 0 or mean_rate < abs(amplitude):
        raise ValueError("need mean_rate > 0 and mean_rate >= |amplitude|")

    cycles = 4.0
    phase = 0.0
    tau = getattr(math, "tau", 2.0 * math.pi)
    dx = length / bins

    rates: List[float] = []
    for i in range(bins):
        x_mid = (i + 0.5) * dx
        rate = mean_rate + amplitude * math.sin(tau * cycles * (x_mid / length) + phase)
        rates.append(max(rate, 0.0))  # guard tiny negatives from FP error

    cum_intensity = [0.0]
    acc = 0.0
    for i in range(bins):
        acc += rates[i] * dx
        cum_intensity.append(acc)

    return RateMap(
        length=int(length),
        bins=int(bins),
        dx=float(dx),
        rates=rates,
        cum_intensity=cum_intensity,
        s_max=cum_intensity[-1],
    )


def make_uniform_rate_map(length: int, bins: int, rate_per_bp: float) -> RateMap:
    """Return a uniform-rate RateMap over [0, length)."""
    if length <= 0 or bins <= 0:
        raise ValueError("length and bins must be positive.")
    if rate_per_bp <= 0:
        raise ValueError("rate_per_bp must be positive.")
    dx = length / bins
    rates = [rate_per_bp] * bins
    cum_intensity = [0.0]
    for i in range(bins):
        cum_intensity.append(cum_intensity[-1] + rate_per_bp * dx)
    return RateMap(
        length=int(length),
        bins=int(bins),
        dx=float(dx),
        rates=rates,
        cum_intensity=cum_intensity,
        s_max=cum_intensity[-1],
    )


class RateMapInterferenceSampler:
    """Gamma-renewal crossover sampler on a binned RateMap.

    Simulates in map units (Λ) with a stationary start, then inverts to bp.
    Interference strength `nu`: 1 => Poisson; larger => stronger interference.
    """

    def __init__(self, rate_map: RateMap) -> None:
        if rate_map.s_max <= 0:
            raise ValueError("rate_map has zero total rate (s_max <= 0)")
        self.rate_map = rate_map

    def _invert_monotone(self, s_vals: List[float]) -> List[int]:
        """Map increasing s-values (Λ) to bp via one forward scan + linear interp."""
        out: List[int] = []
        j = 1  # cum_intensity[0] = 0
        cum_intensity = self.rate_map.cum_intensity
        dx = self.rate_map.dx

        for s in s_vals:
            while j < len(cum_intensity) and cum_intensity[j] < s:
                j += 1
            if j >= len(cum_intensity):
                out.append(self.rate_map.length - 1)
                continue
            s0, s1 = cum_intensity[j - 1], cum_intensity[j]
            x0 = (j - 1) * dx
            frac = 0.0 if s1 == s0 else (s - s0) / (s1 - s0)
            x = x0 + frac * dx
            out.append(int(x))

        return out

    def sample_meiosis_gamma(self, nu: float = 5.0, rng: Optional[np.random.Generator] = None) -> List[int]:
        """Draw one meiosis of crossovers; return sorted integer bp positions."""
        if nu <= 0:
            return []
        r = rng or np.random.default_rng()
        s_max = self.rate_map.s_max

        # Stationary start in map space (mean gap = 1):
        # first event is a random fraction of a gamma-distributed gap.
        gap = r.gamma(shape=nu, scale=1.0 / nu)  # mean 1
        t = gap * r.random()

        s_points: List[float] = []
        while t < s_max:
            s_points.append(t)
            t += r.gamma(shape=nu, scale=1.0 / nu)

        return self._invert_monotone(s_points)


if __name__ == "__main__":
    rate_map = make_sine_rate_map(length=100_000_000, bins=500, mean_rate=2e-8)
    sampler = RateMapInterferenceSampler(rate_map)

    rng = np.random.default_rng(123)

    for i in range(5):
        positions = sampler.sample_meiosis_gamma(nu=5.0, rng=rng)
        print(i, positions)
