"""
TDSParams — central configuration for all TDS computations.

All parameters have defaults matching the original NCOM paper
(Bashan et al., Nature Communications 2012).
"""

from dataclasses import dataclass

VALID_WINDOW_ANCHORS = ("center", "end")
VALID_STABILITY_ANCHORS = ("any", "end")


@dataclass
class TDSParams:
    """
    Configuration for the TDS algorithm. Pass a single instance through
    all functions so settings are consistent end-to-end.

    Parameters
    ----------
    window : int
        Sliding window size in samples. Default 60 (= 60 s at 1 Hz).
    overlap : int
        Step size between windows in samples. Default 30 (50% overlap).
    max_lag : int
        Maximum cross-correlation lag to search, in samples (±max_lag).
        Default 30.
    stability_window : int
        Number of consecutive τ₀ points assessed for stability. Default 5.
    stability_min : int
        Minimum number of points within tolerance to label as stable.
        Default 4 (i.e. 4 out of 5).
    tolerance : int
        Allowed deviation in τ₀ (samples) to still be considered stable.
        Default ±1.
    window_anchor : str
        Where on the window to stamp each t_vec value. ``"center"`` (default)
        uses ``start + L//2``; ``"end"`` uses ``start + L`` (the last sample
        of the window — the earliest moment the result could be known).
    stability_anchor : str
        Which points inside a stable window get labeled. ``"any"`` (default)
        reproduces the Bashan et al. 2012 offline definition and **uses
        look-ahead** — a label at index ``j`` can be set by a window starting
        at ``j``, so it depends on tau up to index ``j + stability_window - 1``.
        ``"end"`` is causal: within each stability window only the last index
        may be labeled, and only if that last point is itself within
        ``tolerance`` of the winning candidate delay; guarantees
        ``stable_label(tau)[:k] == stable_label(tau[:k])``. Note the first
        ``stability_window - 1`` labels are always 0 in ``"end"`` mode
        (warm-up) and that TDS scores in ``"end"`` mode are never higher
        than in ``"any"`` mode.
    n_surrogates : int
        Number of surrogate subjects for null distribution. Default 1000.
    alpha : float
        Significance level for threshold (1-alpha percentile). Default 0.05.
    """

    window: int = 60
    overlap: int = 30
    max_lag: int = 30
    stability_window: int = 5
    stability_min: int = 4
    tolerance: int = 1
    n_surrogates: int = 1000
    alpha: float = 0.05
    window_anchor: str = "center"  # "center" → start + L//2; "end" → start + L
    stability_anchor: str = "any"  # "any" → any point in a stable window; "end" → causal, last point only

    def __post_init__(self) -> None:
        if self.window_anchor not in VALID_WINDOW_ANCHORS:
            raise ValueError(
                f"window_anchor must be one of {VALID_WINDOW_ANCHORS}, "
                f"got {self.window_anchor!r}"
            )
        if self.stability_anchor not in VALID_STABILITY_ANCHORS:
            raise ValueError(
                f"stability_anchor must be one of {VALID_STABILITY_ANCHORS}, "
                f"got {self.stability_anchor!r}"
            )
