from __future__ import annotations

import math

import matplotlib.pyplot as plt
import numpy as np


class AngleHistogram:
    def __init__(self, n_bins: int = 80):
        self.bin_edges = np.linspace(0.0, math.pi, n_bins + 1)
        self.fixed_order = np.zeros(n_bins)
        self.showered = np.zeros(n_bins)

    def fill_fixed_order(self, value: float | None, weight: float) -> None:
        self._fill(self.fixed_order, value, weight)

    def fill_showered(self, value: float | None, weight: float) -> None:
        self._fill(self.showered, value, weight)

    def _fill(self, histogram, value: float | None, weight: float) -> None:
        if value is None or not math.isfinite(value):
            return
        index = np.searchsorted(self.bin_edges, value, side="right") - 1
        if 0 <= index < len(histogram):
            histogram[index] += weight

    def save(self, prefix: str, title: str) -> None:
        centers = 0.5 * (self.bin_edges[:-1] + self.bin_edges[1:])
        widths = self.bin_edges[1:] - self.bin_edges[:-1]

        np.savez(
            f"{prefix}.npz",
            bin_edges=self.bin_edges,
            fixed_order=self.fixed_order,
            showered=self.showered,
        )

        plt.figure(figsize=(7.0, 4.8))
        plt.step(centers, self.fixed_order / widths, where="mid", label="fixed order")
        plt.step(centers, self.showered / widths, where="mid", label="showered")
        plt.xlabel(r"$\theta(j_2,j_3)$")
        plt.ylabel("weighted events / rad")
        plt.title(title)
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{prefix}.png", dpi=160)
        plt.close()

