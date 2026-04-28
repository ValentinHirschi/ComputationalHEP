import math

NC = 3.0
TR = 1.0 / 2.0
CA = NC
CF = (NC * NC - 1.0) / (2.0 * NC)


class AlphaS:
    def __init__(self, mz: float, asmz: float, order: int = 1, mb: float = 4.75, mc: float = 1.3):
        self.order = order
        self.mc2 = mc * mc
        self.mb2 = mb * mb
        self.mz2 = mz * mz
        self.asmz = asmz
        self.asmb = self(self.mb2)
        self.asmc = self(self.mc2)

    def beta0(self, nf: int) -> float:
        return 11.0 / 6.0 * CA - 2.0 / 3.0 * TR * nf

    def beta1(self, nf: int) -> float:
        return 17.0 / 6.0 * CA * CA - (5.0 / 3.0 * CA + CF) * TR * nf

    def as0(self, t: float) -> float:
        if t <= 0.0:
            raise ValueError(f"alpha_s scale must be positive, got {t}")

        if t >= self.mb2:
            tref = self.mz2
            asref = self.asmz
            b0 = self.beta0(5) / (2.0 * math.pi)
        elif t >= self.mc2:
            tref = self.mb2
            asref = self.asmb
            b0 = self.beta0(4) / (2.0 * math.pi)
        else:
            tref = self.mc2
            asref = self.asmc
            b0 = self.beta0(3) / (2.0 * math.pi)
        return 1.0 / (1.0 / asref + b0 * math.log(t / tref))

    def as1(self, t: float) -> float:
        if t <= 0.0:
            raise ValueError(f"alpha_s scale must be positive, got {t}")

        if t >= self.mb2:
            tref = self.mz2
            asref = self.asmz
            b0 = self.beta0(5) / (2.0 * math.pi)
            b1 = self.beta1(5) / (2.0 * math.pi) ** 2
        elif t >= self.mc2:
            tref = self.mb2
            asref = self.asmb
            b0 = self.beta0(4) / (2.0 * math.pi)
            b1 = self.beta1(4) / (2.0 * math.pi) ** 2
        else:
            tref = self.mc2
            asref = self.asmc
            b0 = self.beta0(3) / (2.0 * math.pi)
            b1 = self.beta1(3) / (2.0 * math.pi) ** 2
        w = 1.0 + b0 * asref * math.log(t / tref)
        return asref / w * (1.0 - b1 / b0 * asref * math.log(w) / w)

    def __call__(self, t: float) -> float:
        if self.order == 0:
            return self.as0(t)
        return self.as1(t)

