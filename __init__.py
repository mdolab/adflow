import sys, os
cur_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(cur_dir + '/python')
import pySUmb, pySUmb_C

__all__ = ['pySUmb','pySUMB_C']
