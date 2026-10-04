import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

from adaptive import adaptive_params, decorrelation_time  # noqa: E402
from generators import ar1  # noqa: E402


def test_decorrelation_time_matches_tc():
    for tc in (2, 10, 40):
        x = ar1(30000, tc, np.random.default_rng(tc))
        assert abs(decorrelation_time(x) - tc) <= 0.2 * tc


def test_window_scales_with_tc():
    rng = np.random.default_rng(1)
    w = {}
    for tc in (10, 20):
        sig = [ar1(30000, tc, rng) for _ in range(3)]
        w[tc] = adaptive_params(sig)[0].window
    assert 1.6 <= w[20] / w[10] <= 2.4


def test_calibrated_null_score_low():
    rng = np.random.default_rng(2)
    sig = [ar1(30000, 10, rng) for _ in range(4)]
    p, info = adaptive_params(sig, tolerance="calibrated", n_shifts=5, max_pairs=6)
    assert info["null_score_by_tolerance"][info["tolerance"]] <= 5.0
    # independent pair, not used in calibration
    from pyTDS.core import tds
    a, b = ar1(30000, 10, rng), ar1(30000, 10, rng)
    assert tds(a, b, p)["score"] <= 5.0
