from enum import StrEnum
import logging


class CHEPException(BaseException):
    pass


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
