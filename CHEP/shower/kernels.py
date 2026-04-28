from __future__ import annotations

import math
import random

from CHEP.shower.qcd import CA, CF, TR


class Kernel:
    def __init__(self, flavours: list[int]):
        self.flavours = flavours


class MEKernel(Kernel):
    @staticmethod
    def cs_to_sp(z: float, y: float, m2: float) -> tuple[float, float, float]:
        return y * m2, (z - 1.0) * (y - 1.0) * m2, z * (1.0 - y) * m2


class PqqME(MEKernel):
    def __init__(self, flavours: list[int], power_correction: float = 1.0):
        super().__init__(flavours)
        self.power_correction = power_correction

    def value(self, z: float, y: float, *, m2: float) -> float:
        sij, sjk, sik = self.cs_to_sp(z, y, m2)
        me_kernel = (
            2.0 * sij * (sik + sjk)
            + sij**2
            + 2.0 * sik * sjk
            + 2.0 * sik**2
            + sjk**2
        ) / (sij * sjk)
        propagator = m2 / sij
        partial_fractioner = sjk / (sij + sjk)
        value = CF * ((me_kernel * partial_fractioner) - 1.0 + self.power_correction)
        return value / propagator

    def estimate(self, z: float) -> float:
        return CF * 2.0 / (1.0 - z) + CF * self.power_correction

    def integral(self, zm: float, zp: float) -> float:
        return CF * 2.0 * math.log((1.0 - zm) / (1.0 - zp)) + self.power_correction * CF * (zp - zm)

    def generate_z(self, zm: float, zp: float, rng: random.Random) -> float:
        return 1.0 + (zp - 1.0) * ((1.0 - zm) / (1.0 - zp)) ** rng.random()


class Pgg(Kernel):
    def value(self, z: float, y: float, **opts) -> float:
        return CA / 2.0 * (2.0 / (1.0 - z * (1.0 - y)) - 2.0 + z * (1.0 - z))

    def estimate(self, z: float) -> float:
        return CA / (1.0 - z)

    def integral(self, zm: float, zp: float) -> float:
        return CA * math.log((1.0 - zm) / (1.0 - zp))

    def generate_z(self, zm: float, zp: float, rng: random.Random) -> float:
        return 1.0 + (zp - 1.0) * ((1.0 - zm) / (1.0 - zp)) ** rng.random()


class Pgq(Kernel):
    def value(self, z: float, y: float, **opts) -> float:
        return TR / 2.0 * (1.0 - 2.0 * z * (1.0 - z))

    def estimate(self, z: float) -> float:
        return TR / 2.0

    def integral(self, zm: float, zp: float) -> float:
        return TR / 2.0 * (zp - zm)

    def generate_z(self, zm: float, zp: float, rng: random.Random) -> float:
        return zm + (zp - zm) * rng.random()


def default_kernels() -> list[Kernel]:
    kernels: list[Kernel] = [
        PqqME([flavour, flavour, 21])
        for flavour in [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]
    ]
    kernels += [Pgq([21, flavour, -flavour]) for flavour in [1, 2, 3, 4, 5]]
    kernels += [Pgg([21, 21, 21])]
    return kernels

