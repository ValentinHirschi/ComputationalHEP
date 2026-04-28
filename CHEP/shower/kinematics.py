from __future__ import annotations

import math

from CHEP.utils.vectors import LorentzVector


def spatial_cross(a: LorentzVector, b: LorentzVector) -> LorentzVector:
    cross = a.space().cross(b.space())
    return LorentzVector([0.0, cross[0], cross[1], cross[2]])


def boost_to_rest(reference: LorentzVector, vector: LorentzVector) -> LorentzVector:
    boosted = LorentzVector(vector)
    boosted.boost(-reference.boostVector())
    return boosted


def boost_from_rest(reference: LorentzVector, vector: LorentzVector) -> LorentzVector:
    boosted = LorentzVector(vector)
    boosted.boost(reference.boostVector())
    return boosted


def safe_angle(a: LorentzVector, b: LorentzVector) -> float:
    denominator = abs(a.space()) * abs(b.space())
    if denominator == 0.0:
        return math.nan
    cosine = a.space().dot(b.space()) / denominator
    cosine = min(max(cosine, -1.0), 1.0)
    return math.acos(cosine)

