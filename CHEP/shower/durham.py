from __future__ import annotations

import math

from CHEP.shower.kinematics import safe_angle
from CHEP.utils.vectors import LorentzVector


def yij(p: LorentzVector, q: LorentzVector, ecm2: float) -> float:
    denominator = abs(p.space()) * abs(q.space())
    if denominator == 0.0:
        return math.inf
    cosine = p.space().dot(q.space()) / denominator
    cosine = min(max(cosine, -1.0), 1.0)
    return 2.0 * min(p[0], q[0]) ** 2 * (1.0 - cosine) / ecm2


def cluster_to_n_jets(momenta: list[LorentzVector], n_jets: int, ecm2: float) -> list[LorentzVector]:
    jets = [LorentzVector(momentum) for momentum in momenta]
    if len(jets) <= n_jets:
        return jets

    while len(jets) > n_jets:
        best = None
        best_value = math.inf
        for i in range(len(jets)):
            for j in range(i):
                value = yij(jets[i], jets[j], ecm2)
                if value < best_value:
                    best_value = value
                    best = (i, j)
        if best is None:
            break
        i, j = best
        jets[j] = LorentzVector(jets[j] + jets[i])
        del jets[i]
    return jets


def leading_jet_angle_23(momenta: list[LorentzVector], ecm2: float) -> float | None:
    if len(momenta) < 3:
        return None
    jets = cluster_to_n_jets(momenta, 3, ecm2)
    if len(jets) < 3:
        return None
    jets.sort(key=lambda jet: jet[0], reverse=True)
    angle = safe_angle(jets[1], jets[2])
    if math.isnan(angle):
        return None
    return angle

