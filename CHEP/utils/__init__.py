from enum import StrEnum
import logging
import numpy as np
import random
from CHEP.matrix_elements.epem_lplm.model.parameters import ModelParameters


class CHEPException(BaseException):
    pass


class Parameter(object):
    def __init__(self, name: str, value: float):
        self.name = name
        self.value = value
        self.real = value

    def lower(self) -> str:
        return self.name.lower()


class Particle(object):
    def __init__(self, name: str, mass: str = "zero", width: str = "zero"):
        self.name = name
        self.mass = mass
        self.width = width

    def get(self, characteristics: str) -> str:
        match characteristics:
            case 'mass': return self.mass
            case 'width': return self.width
            case _: raise CHEPException(f"Particle characteristics {characteristics} not implemented")


class ModelInformation(object):

    def __init__(self, model_parameter: ModelParameters):
        self.parameters = model_parameter
        self.parameters_dict = {k: Parameter(k, v)
                                for k, v in self.parameters.__dict__.items()}

    def get_particle(self, pdg: int):

        match abs(pdg):
            case 1: p = Particle('d')
            case 2: p = Particle('u')
            case 3: p = Particle('s')
            case 4: p = Particle('c')
            case 5: p = Particle('b', "mdl_MB")
            case 6: p = Particle('t', "mdl_MT")
            case 11: p = Particle('e-')
            case 12: p = Particle('ve')
            case 13: p = Particle('mu-')
            case 14: p = Particle('vm')
            case 15: p = Particle('ta-', "mdl_MTA")
            case 16: p = Particle('vt')
            case 21: p = Particle('g')
            case 22: p = Particle('a')
            case 23: p = Particle('z', "mdl_MZ", "mdl_WZ")
            case 24: p = Particle('w+', "mdl_MW", "mdl_WW")
            case 25: p = Particle('h', "mdl_MH", "mdl_WH")
            case _: raise CHEPException(f"Particle {pdg} not implemented")

        if pdg < 0:
            match p.name[-1]:
                case '-': p.name = p.name[:-1] + '+'
                case '+': p.name = p.name[:-1] + '-'
                case _: p.name = p.name + '~'

        return p

    def get(self, characteristics: str):
        match characteristics:
            case 'parameter_dict': return self.parameters_dict
            case _: raise CHEPException(f"Model characteristics {characteristics} not implemented")


class Colour(StrEnum):
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'


logger = logging.getLogger('CHEP')


def setup_logging(level: int = logging.INFO) -> None:
    logging.basicConfig(
        format=f'{Colour.GREEN}%(levelname)s{Colour.END} {Colour.BLUE}%(funcName)s l.%(lineno)d{
            Colour.END} {Colour.CYAN}t=%(asctime)s.%(msecs)03d{Colour.END} > %(message)s',
        datefmt='%Y-%m-%d,%H:%M:%S'
    )
    logger.setLevel(level)


class Dimension(object):
    """ A dimension object specifying a specific integration dimension."""

    def __init__(self, name, folded=False):
        self.name = name
        self.folded = folded

    def length(self):
        raise NotImplemented

    def random_sample(self):
        raise NotImplemented


class DiscreteDimension(Dimension):
    """ A dimension object specifying a specific discrete integration dimension."""

    def __init__(self, name, values, **opts):
        try:
            self.normalized = opts.pop('normalized')
        except:
            self.normalized = False
        super(DiscreteDimension, self).__init__(name, **opts)
        assert (isinstance(values, list))
        self.values = values

    def length(self):
        if self.normalized:
            return 1.0/float(len(self.values))
        else:
            return 1.0

    def random_sample(self):
        return np.int64(random.choice(self.values))


class ContinuousDimension(Dimension):
    """ A dimension object specifying a specific discrete integration dimension."""

    def __init__(self, name, lower_bound=0.0, upper_bound=1.0, **opts):
        super(ContinuousDimension, self).__init__(name, **opts)
        assert (upper_bound > lower_bound)
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound

    def length(self):
        return (self.upper_bound-self.lower_bound)

    def random_sample(self):
        return np.float64(self.lower_bound+random.random()*(self.upper_bound-self.lower_bound))


class DimensionList(list):
    """A DimensionList."""

    def __init__(self, *args, **opts):
        super(DimensionList, self).__init__(*args, **opts)

    def volume(self):
        """ Returns the volue of the complete list of dimensions."""
        vol = 1.0
        for d in self:
            vol *= d.length()
        return vol

    def append(self, arg, **opts):
        """ Type-checking. """
        assert (isinstance(arg, Dimension))
        super(DimensionList, self).append(arg, **opts)

    def get_discrete_dimensions(self):
        """ Access all discrete dimensions. """
        return DimensionList(d for d in self if isinstance(d, DiscreteDimension))

    def get_continuous_dimensions(self):
        """ Access all discrete dimensions. """
        return DimensionList(d for d in self if isinstance(d, ContinuousDimension))

    def random_sample(self):
        return np.array([d.random_sample() for d in self])
