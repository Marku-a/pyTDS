"""R1.2: thesis surrogate subjects, ported from msc-project '01 repreduce NCOM/surrogate test/':
GetSurPatient.m (pass 1 + 'add more Data' pass 2 with arr_data(1:end*0.5,:)) and randomizeSurData.m (5 rounds),
then GetLinksStrength on the surrogate (thesis_port.links_strength) -> 140 values (flatten order F).

Usage: python r1_surrogates.py DATA.npz N SEED OUT.npy   (4 processes, per-surrogate seeds from SeedSequence(SEED).spawn(N))

MATLAB semantics decisions: randperm(35,10) -> rng.permutation(35)[:10]; randi(n,[2,1]) -> rng.integers(1,n+1,2);
1:end*0.5 on odd length -> floor(n/2) rows; stage labels = the 10th signal's subject (stage_dum uses jnd==10).
Seeds differ from MATLAB, so only distributions are comparable.
"""
import os, sys
from multiprocessing import Pool
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from thesis_port import links_strength

ORDER = ("LS", "DS", "REM", "awake")        # GetSurPatient stage order (differs from GetLinksStrength's)


def load_subjects(path):
    z = np.load(path)
    return [(z["X_" + n], z["S_" + n]) for n in z["names"]]


def _one_pass(subjects, rng):
    """One 'paitent = randperm(35,10)' pass -> (data n x 10, stages n)."""
    pat = rng.permutation(len(subjects))[:10]
    cols, sts = [], []
    for sens in range(10):
        X, S = subjects[pat[sens]]
        d, s = X[:, sens], S
        if sens > 0:
            n = len(cols[0])
            if len(d) > n:
                d, s = d[:n], s[:n]
            else:
                cols = [c[:len(d)] for c in cols]
                sts = [c[:len(d)] for c in sts]
        cols.append(d); sts.append(s)
    sur = np.stack(cols, 1)
    stg = np.stack(sts, 1)
    data, labs = [], []
    for st in ORDER:
        inst = stg == st
        parts = [sur[inst[:, j], j] for j in range(10)]
        shortest = min(len(p) for p in parts)
        data.append(np.stack([p[:shortest] for p in parts], 1))
        labs.append(np.full(shortest, st, dtype=stg.dtype))
    return np.concatenate(data), np.concatenate(labs)


def randomize(data, stages, rng):
    """randomizeSurData.m: 5 rounds of random-block shuffling (labels move with the data)."""
    or_size = len(data)
    idx = np.arange(or_size)
    for _ in range(5):
        rest, new = idx, []
        while len(rest):
            a, b = rng.integers(1, len(rest) + 1, 2)         # randi(n,[2,1]), 1-based inclusive
            lo, hi = min(a, b), max(a, b)
            new.append(rest[lo - 1:hi])
            rest = np.concatenate([rest[:lo - 1], rest[hi:]])
            if len(rest) < or_size / 20:
                new.append(rest)
                rest = rest[:0]
        idx = np.concatenate(new)
    return data[idx], stages[idx]


def make_surrogate(subjects, rng):
    """GetSurPatient.m. subjects: list of 35 (X T x 10, stages T) from load_subjects. Returns (data, stages)."""
    d1, s1 = _one_pass(subjects, rng)
    d2, s2 = _one_pass(subjects, rng)
    h = len(d2) // 2                                         # 1:end*0.5 -> floor
    data = np.concatenate([d1, d2[:h]])
    stages = np.concatenate([s1, s2[:h]])
    return randomize(data, stages, rng)


_SUBJ = None


def _work(args):
    global _SUBJ
    path, seed = args
    if _SUBJ is None:
        _SUBJ = load_subjects(path)
    data, stages = make_surrogate(_SUBJ, np.random.default_rng(seed))
    return links_strength(data, stages).flatten(order="F")


if __name__ == "__main__":
    path, N, seed, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    seeds = np.random.SeedSequence(seed).spawn(N)
    with Pool(4) as p:
        res = np.array(p.map(_work, [(path, s) for s in seeds], chunksize=4))
    np.save(out, res)
    print(res.shape, "nan:", int(np.isnan(res).sum()), "mean:", np.nanmean(res))
