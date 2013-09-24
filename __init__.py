import sys, os
cur_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(cur_dir, 'python'))
import pySUmb
import pySUmb_C
from pySUmb import SUMB
from pySUmb_C import SUMB_C
__all__ = ['pySUmb', 'pySUmb_C', 'SUMB', 'SUMB_C']
