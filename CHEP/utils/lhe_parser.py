import os
import sys
from enum import Enum
try:
    if os.getenv('MG_ROOT_PATH') is not None:
        sys.path.insert(0, str(os.getenv('MG_ROOT_PATH')))
    else:
        sys.path.append('../')
    from madgraph.various.lhe_parser import EventFile, Event, Particle  # type: ignore
except Exception as e:
    print(f"Could not load Magraph lhe parser: {e}. Specify madgraph root path with env. variable MG_ROOT_PATH.")  # nopep8
    sys.exit(1)
from CHEP.utils import CHEP_TEMPLATES


class CHEPEventFile(EventFile):
    def __init__(self, file_path, *args, mode='r', **kwargs):
        super().__init__(file_path, *args, mode=mode, **kwargs)
        if mode == 'w':
            self.banner = EventFile(os.path.join(CHEP_TEMPLATES, "template_events_file.lhe"), mode='r').get_banner()  # nopep8
            self.banner.write(self, close_tag=False)
            self.write("\n"*10)
        self.banner = self.get_banner()


class CHEPEvent(Event):
    pass


class LegState(Enum):
    INITIAL = -1
    FINAL = 1

    def to_lhe(self) -> float:
        return float(self.value)


class CHEPParticle(Particle):
    def __init__(self, event: CHEPEvent, status: LegState, pid: int, px: float, py: float, pz: float, E: float, mass: float, color1: int = 0, color2: int | None = 0):
        super().__init__(event)
        self.pid = pid
        self.color1 = color1
        self.color2 = color2
        self.px = px
        self.py = py
        self.pz = pz
        self.E = E
        self.mass = mass
        self.status = status.to_lhe()
        self.mass = 0
        self.vtim = 0
        self.helicity = 9
        self.rwgt = 0
        self.comment = ''
