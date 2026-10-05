"""Faithful Python port of the thesis MATLAB TDS (msc-project, '01 repreduce NCOM').

Time_Delay_Interaction.m : L=60, overlap=30, NL=floor(2N/L-1), segments zscored (MATLAB zscore, ddof=1),
max_xc_time_delay.m      : xcorr(s1,s2,30,'coeff'), tau = lag of max |xc| (first on ties),
StableLabel.m            : 4 of 5 consecutive taus within +-1 of one of the window's delays.
GetLinksStrength.m       : per stage, % stable windows whose time (1-based sample (k-1)*30+30) lies in the stage.
"""
import numpy as np

STAGES = ("DS", "LS", "REM", "awake")          # column order of GetLinksStrength
SIGNALS = ("HR", "Resp", "Chin", "Leg", "Eye", "delta", "theta", "alpha", "sigma", "beta")


def _good():
    """triu(ones(10))-eye(10) with the EEG-EEG block (signals 6..10) removed: 35 links."""
    g = np.triu(np.ones((10, 10)), 1).astype(bool)
    g[5:, 5:] = False
    return g


def _order(g):
    """MATLAB maps(good_links') order: column-major over good_links', i.e. element maps(r, c) where g[c, r]."""
    gt = g.T
    return [(r, c) for c in range(10) for r in range(10) if gt[r, c]]


def zs(x):
    sd = x.std(ddof=1)
    return (x - x.mean()) / sd if sd > 0 else np.full_like(x, np.nan)


def xcorr_coeff(a, b, maxlag):
    """MATLAB xcorr(a,b,maxlag,'coeff'): R(m) = sum_n a(n+m) b(n) / sqrt(sum a^2 sum b^2), m=-maxlag..maxlag."""
    full = np.correlate(a, b, mode="full")       # index k <-> lag k-(N-1), c[k] = sum a[n+lag] b[n]
    N = len(a)
    mid = N - 1
    r = full[mid - maxlag: mid + maxlag + 1]
    den = np.sqrt(np.sum(a * a) * np.sum(b * b))
    return r / den, np.arange(-maxlag, maxlag + 1)


def time_delay_interaction(s1, s2, L=60, overlap=30, maxlag=30):
    N = len(s1)
    NL = int(np.floor(2 * N / L - 1))
    tau = np.full(NL, np.nan)
    t_vec = np.zeros(NL, int)
    for k in range(NL):
        a = k * overlap
        seg1, seg2 = zs(s1[a:a + L]), zs(s2[a:a + L])
        t_vec[k] = a + L // 2            # 1-based seg_inds(end/2) = a+30
        if np.isnan(seg1).any() or np.isnan(seg2).any():
            xc = np.full(2 * maxlag + 1, np.nan)
        else:
            xc, lags = xcorr_coeff(seg1, seg2, maxlag)
        ax = np.abs(xc)
        tau[k] = (np.nanargmax(ax) - maxlag) if np.isfinite(ax).any() else -maxlag  # MATLAB max of all-NaN -> index 1
    return tau, t_vec


def stable_label(tau, tol=1):
    N = len(tau)
    lbl = np.zeros(N, int)
    for i in range(N - 4):
        seg = tau[i:i + 5]
        for d in np.unique(seg):
            inr = (seg >= d - tol) & (seg <= d + tol)
            if inr.sum() >= 4:
                lbl[i:i + 5][inr] = 1
                break
    return lbl


def links_strength(X, stages):
    """X: T x 10, stages: array of str per sample. Returns 35 x 4 array (rows good-link order, cols STAGES)."""
    stages = np.asarray(stages)
    out = np.zeros((35, 4))
    for li, (i, j) in enumerate(_order(_good())):
        tau, t = time_delay_interaction(X[:, i], X[:, j])
        lbl = stable_label(tau)
        st_at = stages[t - 1]  # MATLAB time t (1-based) -> 0-based t-1
        for k, s in enumerate(STAGES):
            m = st_at == s
            out[li, k] = 100 * lbl[m].mean() if m.any() else np.nan
    return out



def links_counts(X, stages):
    """Like links_strength but returns 35 x 4 x 2: (sum of stable labels, n windows) per link x stage, for pooling
    over subjects as makeDataForReconFig2.m ('Divide to sleep stages') does (centre-sample stage rule)."""
    stages = np.asarray(stages)
    out = np.zeros((35, 4, 2))
    for li, (i, j) in enumerate(_order(_good())):
        tau, t = time_delay_interaction(X[:, i], X[:, j])
        lbl = stable_label(tau)
        st_at = stages[t - 1]
        for k, s in enumerate(STAGES):
            m = st_at == s
            out[li, k] = (lbl[m].sum(), m.sum())
    return out
