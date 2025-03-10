import os
import sys
try:
    if os.getenv('MG_ROOT_PATH') is not None:
        sys.path.insert(0, str(os.getenv('MG_ROOT_PATH')))
    else:
        sys.path.append('../')
    from madgraph.various.lhe_parser import EventFile  # type: ignore
except Exception as e:
    print(f"Could not load Magraph lhe parser: {e}. Specify madgraph root path with env. variable MG_ROOT_PATH.")  # nopep8
    sys.exit(1)


class CHEPEventFile(EventFile):
    pass
