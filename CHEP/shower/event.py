from __future__ import annotations

from dataclasses import dataclass

from CHEP.utils.lhe_parser import CHEPEvent, CHEPParticle, LegState
from CHEP.utils.vectors import LorentzVector, LorentzVectorList


@dataclass
class ShowerParticle:
    pid: int
    momentum: LorentzVector
    colors: list[int]

    def color_connected(self, other: "ShowerParticle") -> bool:
        return (
            self.colors[0] > 0 and self.colors[0] == other.colors[1]
        ) or (
            self.colors[1] > 0 and self.colors[1] == other.colors[0]
        )

    def set(self, pid: int, momentum: LorentzVector, colors: list[int]) -> None:
        self.pid = pid
        self.momentum = LorentzVector(momentum)
        self.colors = list(colors)


@dataclass
class ShowerEvent:
    incoming: list[ShowerParticle]
    final: list[ShowerParticle]
    weight: float
    aqed: float
    aqcd: float
    scale: float = 0.0
    ievent: int = 0

    def all_particles(self) -> list[ShowerParticle]:
        return self.incoming + self.final

    def final_momenta(self) -> LorentzVectorList:
        return LorentzVectorList([particle.momentum for particle in self.final])

    def hard_scale2(self) -> float:
        total = LorentzVector()
        for particle in self.incoming:
            total += particle.momentum
        return total.square()

    def max_color(self) -> int:
        colors = [c for particle in self.final for c in particle.colors if c > 0]
        return max(colors, default=0)

    def to_chep_event(self) -> CHEPEvent:
        event = CHEPEvent()
        event.nexternal = len(self.incoming) + len(self.final)
        event.ievent = self.ievent
        event.wgt = self.weight
        event.scale = self.scale
        event.aqed = self.aqed
        event.aqcd = self.aqcd

        for particle in self.incoming:
            event.append(
                CHEPParticle(
                    event,
                    LegState.INITIAL,
                    particle.pid,
                    particle.momentum[1],
                    particle.momentum[2],
                    particle.momentum[3],
                    particle.momentum[0],
                    mass=0.0,
                    color1=particle.colors[0],
                    color2=particle.colors[1],
                )
            )
        for particle in self.final:
            event.append(
                CHEPParticle(
                    event,
                    LegState.FINAL,
                    particle.pid,
                    particle.momentum[1],
                    particle.momentum[2],
                    particle.momentum[3],
                    particle.momentum[0],
                    mass=0.0,
                    color1=particle.colors[0],
                    color2=particle.colors[1],
                )
            )
        return event


def momentum_from_lhe_particle(particle) -> LorentzVector:
    return LorentzVector([particle.E, particle.px, particle.py, particle.pz])


def shower_event_from_chep_event(event: CHEPEvent) -> ShowerEvent:
    incoming = []
    final = []
    for particle in event:
        shower_particle = ShowerParticle(
            pid=particle.pid,
            momentum=momentum_from_lhe_particle(particle),
            colors=[int(particle.color1 or 0), int(particle.color2 or 0)],
        )
        if particle.status == LegState.INITIAL.to_lhe():
            incoming.append(shower_particle)
        elif particle.status == LegState.FINAL.to_lhe():
            final.append(shower_particle)

    return ShowerEvent(
        incoming=incoming,
        final=final,
        weight=event.wgt,
        aqed=event.aqed,
        aqcd=event.aqcd,
        scale=event.scale,
        ievent=event.ievent,
    )


def make_epem_event(
    momenta: list[LorentzVector],
    final_pids: list[int],
    final_colors: list[list[int]],
    weight: float,
    aqed: float,
    aqcd: float,
    scale: float,
) -> CHEPEvent:
    event = CHEPEvent()
    event.nexternal = 2 + len(final_pids)
    event.wgt = weight
    event.scale = scale
    event.aqed = aqed
    event.aqcd = aqcd

    incoming = [
        (-11, momenta[0]),
        (11, momenta[1]),
    ]
    for pid, momentum in incoming:
        event.append(
            CHEPParticle(
                event,
                LegState.INITIAL,
                pid,
                momentum[1],
                momentum[2],
                momentum[3],
                momentum[0],
                mass=0.0,
            )
        )

    for pid, colors, momentum in zip(final_pids, final_colors, momenta[2:]):
        event.append(
            CHEPParticle(
                event,
                LegState.FINAL,
                pid,
                momentum[1],
                momentum[2],
                momentum[3],
                momentum[0],
                mass=0.0,
                color1=colors[0],
                color2=colors[1],
            )
        )
    return event

