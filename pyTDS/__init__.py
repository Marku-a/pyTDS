"""
pyTDS — Time Delay Stability (TDS) library
==========================================
General-purpose algorithm for detecting stable coupling between time-series signals.

Based on:
  Bashan et al., Nature Communications 3:702 (2012)
  Bartsch et al., PLoS ONE (2015)

Original MATLAB implementation by Ronny Bartsch.
Python port and generalization by Asaf Markuza.

Quick start
-----------
>>> import pyTDS
>>> result = pyTDS.tds(s1, s2)
>>> print(result['score'])   # TDS % in [0, 100]
"""

from .params import TDSParams
from .core import tds, tds_score, stable_label, time_delay_interaction
from .network import tds_matrix, fix_symmetry, apply_threshold, to_networkx
from .surrogate import surrogate_tds, significance_threshold

__all__ = [
    "TDSParams",
    "tds",
    "tds_score",
    "stable_label",
    "time_delay_interaction",
    "tds_matrix",
    "fix_symmetry",
    "apply_threshold",
    "to_networkx",
    "surrogate_tds",
    "significance_threshold",
]
