from __future__ import annotations

from dataclasses import dataclass
import math
import random

from CHEP.shower.event import ShowerEvent, ShowerParticle
from CHEP.shower.kernels import Kernel, default_kernels
from CHEP.shower.kinematics import boost_from_rest, boost_to_rest, spatial_cross
from CHEP.utils.vectors import LorentzVector


def clone_particles(particles: list[ShowerParticle]) -> list[ShowerParticle]:
    return [
        ShowerParticle(particle.pid, LorentzVector(particle.momentum), list(particle.colors))
        for particle in particles
    ]


def invariant_mass(momentum: LorentzVector) -> float:
    return math.sqrt(abs(momentum.square()))


def format_particle_table(
    title: str,
    particles: list[ShowerParticle],
    *,
    n_incoming: int = 0,
    subtract_final: bool = False,
) -> str:
    columns = [5, 4, 8, 8, 8, 23, 23, 23, 23, 23]
    template = " ".join(f"%-{width}s" for width in columns)
    line = "-" * (sum(columns) + len(columns) - 1)
    out = [title]
    out.append(template % ("#", "st", "pid", "col1", "col2", "E", "p_x", "p_y", "p_z", "M"))
    out.append(line)

    residual = LorentzVector()
    for i_particle, particle in enumerate(particles, start=1):
        momentum = particle.momentum
        is_incoming = i_particle <= n_incoming
        status = "I" if is_incoming else "F"
        if subtract_final and not is_incoming:
            residual -= momentum
        else:
            residual += momentum
        out.append(
            template
            % (
                i_particle,
                status,
                particle.pid,
                particle.colors[0],
                particle.colors[1],
                f"{momentum[0]: .16e}",
                f"{momentum[1]: .16e}",
                f"{momentum[2]: .16e}",
                f"{momentum[3]: .16e}",
                f"{invariant_mass(momentum): .16e}",
            )
        )

    out.append(line)
    out.append(
        template
        % (
            "Sum",
            "",
            "",
            "",
            "",
            f"{residual[0]: .16e}",
            f"{residual[1]: .16e}",
            f"{residual[2]: .16e}",
            f"{residual[3]: .16e}",
            f"{invariant_mass(residual): .16e}",
        )
    )
    return "\n".join(out)


@dataclass
class EmissionInfo:
    t: float
    z: float
    y: float
    phi: float
    splitter_before: ShowerParticle
    spectator_before: ShowerParticle
    emitted: ShowerParticle
    splitter_after: ShowerParticle
    spectator_after: ShowerParticle
    full_event_before: list[ShowerParticle]
    full_event_after: list[ShowerParticle]
    n_incoming: int

    def format(self) -> str:
        local_before = [self.splitter_before, self.spectator_before]
        local_after = [self.splitter_after, self.emitted, self.spectator_after]
        n_final_before = len(self.full_event_before) - self.n_incoming
        n_final_after = len(self.full_event_after) - self.n_incoming
        return (
            f"t={self.t:.10e}, z={self.z:.10e}, y={self.y:.10e}, phi={self.phi:.10e}\n"
            f"event final-state multiplicity: {n_final_before} -> {n_final_after}\n"
            f"splitter {self.splitter_before.pid} -> {self.splitter_after.pid}, emitted {self.emitted.pid}, "
            f"spectator {self.spectator_before.pid}\n"
            f"new colors: splitter={self.splitter_after.colors}, emitted={self.emitted.colors}, "
            f"spectator={self.spectator_after.colors}\n"
            f"{format_particle_table('before local kinematics:', local_before)}\n"
            f"{format_particle_table('after local kinematics:', local_after)}\n"
            f"{format_particle_table('full event before emission:', self.full_event_before, n_incoming=self.n_incoming, subtract_final=True)}\n"
            f"{format_particle_table('full event after emission:', self.full_event_after, n_incoming=self.n_incoming, subtract_final=True)}"
        )


class Shower:
    def __init__(self, alpha, t0: float, seed: int = 0, kernels: list[Kernel] | None = None):
        if t0 <= 0.0:
            raise ValueError(f"Shower cutoff t0 must be positive, got {t0}")
        self.t0 = t0
        self.alpha = alpha
        self.alphamax = alpha(t0)
        self.kernels = kernels if kernels is not None else default_kernels()
        self.rng = random.Random(seed)
        self.current_t = t0
        self.next_color = 0

    def make_kinematics(
        self,
        z: float,
        y: float,
        phi: float,
        pijt: LorentzVector,
        pkt: LorentzVector,
    ) -> list[LorentzVector]:
        q = pijt + pkt
        q2 = q.square()
        rkt2 = q2 * y * z * (1.0 - z)
        rkt = math.sqrt(max(rkt2, 0.0))

        kt1 = spatial_cross(pijt, pkt)
        if kt1.rho() < 1.0e-12:
            kt1 = spatial_cross(pijt, LorentzVector([0.0, 1.0, 0.0, 0.0]))
        kt1 *= rkt * math.cos(phi) / kt1.rho()

        pijt_cms = boost_to_rest(q, pijt)
        kt2cms = spatial_cross(pijt_cms, kt1)
        kt2cms *= rkt * math.sin(phi) / kt2cms.rho()
        kt2 = boost_from_rest(q, kt2cms)

        pi = LorentzVector(z * pijt + (1.0 - z) * y * pkt + kt1 + kt2)
        pj = LorentzVector((1.0 - z) * pijt + z * y * pkt - kt1 - kt2)
        pk = LorentzVector((1.0 - y) * pkt)
        return [pi, pj, pk]

    def make_colors(
        self,
        flavours: list[int],
        splitter_colors: list[int],
        spectator_colors: list[int],
    ) -> list[list[int]]:
        self.next_color += 1
        new_color = self.next_color

        if flavours[0] != 21:
            if flavours[0] > 0:
                return [[new_color, 0], [splitter_colors[0], new_color]]
            return [[0, new_color], [new_color, splitter_colors[1]]]

        if flavours[1] == 21:
            if splitter_colors[0] == spectator_colors[1]:
                if splitter_colors[1] == spectator_colors[0] and self.rng.random() > 0.5:
                    return [[splitter_colors[0], new_color], [new_color, splitter_colors[1]]]
                return [[new_color, splitter_colors[1]], [splitter_colors[0], new_color]]
            return [[splitter_colors[0], new_color], [new_color, splitter_colors[1]]]

        if flavours[1] > 0:
            return [[splitter_colors[0], 0], [0, splitter_colors[1]]]
        return [[0, splitter_colors[1]], [splitter_colors[0], 0]]

    def generate_emission(self, event: ShowerEvent) -> EmissionInfo | None:
        while self.current_t > self.t0:
            candidate_t = self.t0
            candidate = None

            for splitter in event.final:
                for spectator in event.final:
                    if spectator is splitter or not splitter.color_connected(spectator):
                        continue
                    for kernel in self.kernels:
                        if kernel.flavours[0] != splitter.pid:
                            continue
                        m2 = (splitter.momentum + spectator.momentum).square()
                        if m2 < 4.0 * self.t0:
                            continue
                        zp = 0.5 * (1.0 + math.sqrt(1.0 - 4.0 * self.t0 / m2))
                        integral = kernel.integral(1.0 - zp, zp)
                        if integral <= 0.0:
                            continue
                        g = self.alphamax / (2.0 * math.pi) * integral
                        trial_t = self.current_t * self.rng.random() ** (1.0 / g)
                        if trial_t > candidate_t:
                            candidate_t = trial_t
                            candidate = (splitter, spectator, kernel, m2, zp)

            self.current_t = candidate_t
            if candidate is None or candidate_t <= self.t0:
                return None

            splitter, spectator, kernel, m2, zp = candidate
            z = kernel.generate_z(1.0 - zp, zp, self.rng)
            y = candidate_t / m2 / z / (1.0 - z)
            if y >= 1.0:
                continue

            weight = (1.0 - y) * self.alpha(candidate_t) / self.alphamax
            weight *= kernel.value(z, y, m2=m2) / kernel.estimate(z)
            if weight <= self.rng.random():
                continue

            phi = 2.0 * math.pi * self.rng.random()
            full_event_before = clone_particles(event.incoming + event.final)
            splitter_before = ShowerParticle(splitter.pid, LorentzVector(splitter.momentum), list(splitter.colors))
            spectator_before = ShowerParticle(spectator.pid, LorentzVector(spectator.momentum), list(spectator.colors))

            moms = self.make_kinematics(z, y, phi, splitter.momentum, spectator.momentum)
            colors = self.make_colors(kernel.flavours, splitter.colors, spectator.colors)
            emitted = ShowerParticle(kernel.flavours[2], moms[1], colors[1])
            event.final.append(emitted)
            splitter.set(kernel.flavours[1], moms[0], colors[0])
            spectator.momentum = moms[2]
            full_event_after = clone_particles(event.incoming + event.final)

            return EmissionInfo(
                t=candidate_t,
                z=z,
                y=y,
                phi=phi,
                splitter_before=splitter_before,
                spectator_before=spectator_before,
                emitted=emitted,
                splitter_after=ShowerParticle(splitter.pid, LorentzVector(splitter.momentum), list(splitter.colors)),
                spectator_after=ShowerParticle(spectator.pid, LorentzVector(spectator.momentum), list(spectator.colors)),
                full_event_before=full_event_before,
                full_event_after=full_event_after,
                n_incoming=len(event.incoming),
            )
        return None

    def run(self, event: ShowerEvent, hard_scale2: float | None = None, step: bool = False) -> list[EmissionInfo]:
        self.next_color = event.max_color()
        self.current_t = hard_scale2 if hard_scale2 is not None else event.hard_scale2()
        emissions = []
        while self.current_t > self.t0:
            if step:
                try:
                    input(f"Press Enter to evolve from t={self.current_t:.10e}...")
                except EOFError:
                    print("Input closed; continuing shower without further pauses.")
                    step = False
            emission = self.generate_emission(event)
            if emission is None:
                if step:
                    print(f"No further emission above t0={self.t0:.10e}.")
                break
            emissions.append(emission)
            if step:
                print(emission.format())
        return emissions
